#!/usr/bin/env Rscript

## Script for processing paied reads (Illumina - AVITI) into ASV tables

# Load argparse library to parse arguments
suppressPackageStartupMessages({
  library(argparse)
  })

## Main program starts here
###############################################################################

# Arguments parsing
parser <- ArgumentParser()
parser$add_argument("-i", "--input", help = "input directory path where raw reads are stored")
parser$add_argument("-m", "--intermediates", help = "intermediate files directory path")
parser$add_argument("-o", "--output", help = "output directory path")
parser$add_argument("--run_name", help = "run name (prefix for output files)")
parser$add_argument("--forward_reads_regex", default = ".*R1(_001)?.fastq.gz", help = "regex to find forward reads files")
parser$add_argument("--reverse_reads_regex", default = ".*R2(_001)?.fastq.gz", help = "regex to find reverse reads files")
parser$add_argument("--run_regex", help = "regex to extract run code from raw reads filenames")
parser$add_argument("--sample_regex", help = "regex to extract sample code from raw reads filenames")
parser$add_argument("--database_file", help = "database for classification")
parser$add_argument("--forward_primer", help = "Forward primer")
parser$add_argument("--reverse_primer", help = "Reverse primer")
parser$add_argument("--min_amplicon_size", help = "Minimum amplicon size to keep")
parser$add_argument("--max_amplicon_size", help = "Maximum amplicon size to keep")
parser$add_argument("--processors", default = 12, help = "Number of processors for ASV classification and alignment")
parser$add_argument("--clustering_treshold", default = 0.03, help = "Number of processors for ASV classification and alignment")
parser$add_argument("--amplicon_size_step", default = 3, help = "Amplicon size step")
parser$add_argument("--cluster_by_size",  help = "Cluster sequence after dividing them by sequence length")
parser$add_argument("--translation_check", default = FALSE, help = "Do a translation amplicon error check, TRUE or FALSE")
parser$add_argument("--translation_frame", default = 1, help = "Translation frame for amplicon error check, only active if translation_check = TRUE")


arguments <- parser$parse_args()

INPUT_FOLDER = arguments$input
INTERMEDIATE_FOLDER = arguments$intermediates
OUTPUT_FOLDER = arguments$output
RUN_NAME = arguments$run_name
FORWARD_READS_REGEX = arguments$forward_reads_regex
REVERSE_READS_REGEX = arguments$reverse_reads_regex
RUN_REGEX = arguments$run_regex
SAMPLE_REGEX = arguments$sample_regex
DATABASE_PATH = arguments$database_file
FWD = arguments$forward_primer
REV = arguments$reverse_primer
min_amplicon_size = arguments$min_amplicon_size
max_amplicon_size = arguments$max_amplicon_size
amplicon_size_step = arguments$amplicon_size_step
PROCESSORS = arguments$processors
CLUSTERING_TRESHOLD = arguments$clustering_treshold
CLUSTER_BY_SIZE = arguments$cluster_by_size
TRANSLATION_CHECK = arguments$translation_check
TRANSLATION_FRAME = arguments$translation_frame


# Save workspace image with all the parameters
save.image("./intermediates/workspace.RData")


# Load libraryes and custom scripts
suppressPackageStartupMessages({
  options(tidyverse.quiet = TRUE)
  library(tidyverse, quietly = TRUE)
  library(magrittr, quietly = TRUE)
  library(dada2, quietly = TRUE)
  library(ShortRead, quietly = TRUE)
  library(Biostrings, quietly = TRUE)
  library(DECIPHER, quietly = TRUE)
  library(ape, quietly = TRUE)
  source("./code/999_utils.R")
  })

#### Get the samples organized ####

# define paths
raw_path = file.path(INPUT_FOLDER)
trim_path = file.path(INTERMEDIATE_FOLDER, "01_trim")
filt_path = file.path(INTERMEDIATE_FOLDER, "02_filt")
backup_path = file.path(INTERMEDIATE_FOLDER, "backups")
results_path = file.path(OUTPUT_FOLDER)

# create paths if missing
if (!dir.exists(filt_path)) dir.create(filt_path, recursive = TRUE)
if (!dir.exists(trim_path)) dir.create(trim_path, recursive = TRUE)
if (!dir.exists(backup_path)) dir.create(backup_path, recursive = TRUE)
if (!dir.exists(results_path)) dir.create(results_path, recursive = TRUE)


# find files
sample_table <- tibble::tibble(
  fastq_R1 = sort(list.files(raw_path, FORWARD_READS_REGEX, 
                             full.names = TRUE)),
  fastq_R2 = sort(list.files(raw_path, REVERSE_READS_REGEX, 
                             full.names = TRUE))
) %>% dplyr::mutate(
  runcode = stringr::str_extract(fastq_R1, RUN_REGEX), # Extract run code
  sample = stringr::str_extract(fastq_R1, SAMPLE_REGEX) # Extract sample code
)

# generate filenames for trimmed and filtered reads
sample_table <- sample_table %>% dplyr::mutate(
  trim_R1 = file.path(trim_path, paste(runcode, 
                                       sample, 
                                       "R1_trim.fastq.gz", 
                                       sep = "_")),
  trim_R2 = file.path(trim_path, paste(runcode, 
                                       sample, 
                                       "R2_trim.fastq.gz", 
                                       sep = "_")),
  filt_R1 = file.path(filt_path, paste(runcode, 
                                       sample, 
                                       "R1_filt.fastq.gz", 
                                       sep = "_")),
  filt_R2 = file.path(filt_path, paste(runcode, 
                                       sample, 
                                       "R2_filt.fastq.gz", 
                                       sep = "_"))
)

# Check that there are no duplicates
assertthat::assert_that(
  !any(duplicated(sample_table$fastq_R1)),
  !any(duplicated(sample_table$fastq_R2)),
  !any(duplicated(sample_table$trim_R1)),
  !any(duplicated(sample_table$trim_R2)),
  !any(duplicated(sample_table$filt_R1)),
  !any(duplicated(sample_table$filt_R2))
)

write.table(sample_table, 
            file = paste0(backup_path,"/filenames_table.txt"), 
            sep = "\t", 
            quote = FALSE, 
            col.names = TRUE)

cat("Found", nrow(sample_table), "samples.\n")


# Flags for Cutadapt and trimming primers, filtering reads without primers
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

CUTADAPT_FLAGS <- paste0("cutadapt --cores ",PROCESSORS," -a ^", 
                         FWD, 
                         "...", 
                         REV.RC, 
                         " -A ^", 
                         REV, 
                         "...", 
                         FWD.RC, 
                         " --discard-untrimmed -o")

for (i in seq_along(sample_table$fastq_R1)) {
  system(paste(CUTADAPT_FLAGS, 
               sample_table$trim_R1[i],
               "-p", 
               sample_table$trim_R2[i], 
               sample_table$fastq_R1[i], 
               sample_table$fastq_R2[i]))
}

# Quality-filtering of reads
out <- filterAndTrim(sample_table$trim_R1, 
                     sample_table$filt_R1,
                     sample_table$trim_R2, 
                     sample_table$filt_R2,
  maxN = 0, truncQ = 2, rm.phix = TRUE,
  compress = TRUE, multithread = PROCESSORS, minLen = 190
)

# Error-rate estimation, merging and chimera filtering done separ. for each run
error_results <- list()
backups <- list()
for (i in 1:length(unique(sample_table$runcode))) {
  sample_table_temp <- subset(sample_table, 
                              runcode == unique(sample_table$runcode)[i])

  ## calculating errors - on a subset of data to save time 
  ## see(https://benjjneb.github.io/dada2/bigdata.html)
  set.seed(100)
  errF <- learnErrors(sample_table_temp$filt_R1, 
                      multithread = PROCESSORS, nbases = 1e+8, randomize = TRUE)
  errR <- learnErrors(sample_table_temp$filt_R2, 
                      multithread = PROCESSORS, nbases = 1e+8, randomize = TRUE)
  # Sample inference
  dadaFs <- dada(sample_table_temp$filt_R1, err = errF, multithread = PROCESSORS)
  dadaRs <- dada(sample_table_temp$filt_R2, err = errR, multithread = PROCESSORS)
  # Merging paired reads
  mergers <- mergePairs(dadaFs, 
                        sample_table_temp$filt_R1, 
                        dadaRs, 
                        sample_table_temp$filt_R2, 
                        minOverlap = 8)
  # Creating sequence table
  seqtab <- makeSequenceTable(mergers)
  # Removing chimeras
  seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                      method = "consensus", 
                                      multithread = PROCESSORS, 
                                      verbose = TRUE)
  error_results[[i]] <- seqtab.nochim
  backups[[i]] <- list(errF, errR, dadaFs, dadaRs, mergers, seqtab)

  rm(sample_table_temp)
}

# Merge the results from different runs together
if (length(error_results) == 1) {
  seqtab.nochim <- error_results[[1]]
}
if (length(error_results) != 1) {
  seqtab.nochim <- dada2::mergeSequenceTables(tables = error_results)}


# Save the sequences table
saveRDS(seqtab.nochim, file.path(backup_path, "savepoint1_seqtable.rds"))
# Read the sequence table (only a part of the script can be run from here)
seqtab.nochim <- readRDS(file.path(backup_path, "savepoint1_seqtable.rds"))


# Abundance filtering
seqtab.nochim_min_abundance <- 
  seqtab.nochim[, apply(seqtab.nochim, 2, sum) >= 10] # filtering reads with counts < 10

if(CLUSTER_BY_SIZE == TRUE){
# Run translation check and clustering separately for different lengths to save time and memory
  results_clustering = list()
  amplicon_sizes = seq(as.numeric(min_amplicon_size), as.numeric(max_amplicon_size), by = as.numeric(amplicon_size_step))

  for(i in 1:length(amplicon_sizes)){
    seqtab_temp = seqtab.nochim_min_abundance[ ,nchar(colnames(seqtab.nochim_min_abundance)) == amplicon_sizes[i]]

    if(TRANSLATION_CHECK == TRUE){seqtab_temp_trans = translation_filtering(seqtab_temp, frame = TRANSLATION_FRAME)}
    if(TRANSLATION_CHECK == FALSE){seqtab_temp_trans = seqtab_temp}
  
    seqtab_temp_clust = cluster_seqs(seqtab_temp_trans, nproc = PROCESSORS, cutoff = CLUSTERING_TRESHOLD, align = F)
    seqtab_temp_clust[[1]]$otu = paste0(seqtab_temp_clust[[1]]$otu,"_",amplicon_sizes[i])
    seqtab_temp_clust[[2]]$otu = paste0(seqtab_temp_clust[[2]]$otu,"_",amplicon_sizes[i])
    results_clustering[[i]] = seqtab_temp_clust
    rm(seqtab_temp)
    rm(seqtab_temp_trans)
    rm(seqtab_temp_clust)
    }
}

if(CLUSTER_BY_SIZE == FALSE){
  results_clustering = list()
  amplicon_sizes = seq(as.numeric(min_amplicon_size), as.numeric(max_amplicon_size), by = as.numeric(amplicon_size_step))

  seqtab_temp = seqtab.nochim_min_abundance[ ,nchar(colnames(seqtab.nochim_min_abundance)) %in% amplicon_sizes]
  seqtab_temp_trans = seqtab_temp

  seqtab_temp_clust = cluster_seqs(seqtab_temp_trans, nproc = PROCESSORS, cutoff = CLUSTERING_TRESHOLD, align = TRUE)
  seqtab_temp_clust[[1]]$otu = paste0(seqtab_temp_clust[[1]]$otu)
  seqtab_temp_clust[[2]]$otu = paste0(seqtab_temp_clust[[2]]$otu)
  results_clustering[[1]] = seqtab_temp_clust
  rm(seqtab_temp)
  rm(seqtab_temp_trans)
  rm(seqtab_temp_clust)
}


# Extract OTU table
results_clustering_OTUs = lapply(results_clustering, FUN = function(x){x[[1]]}) 
otu_seqtab = do.call(rbind,results_clustering_OTUs)
otu_seqtab = otu_seqtab[order(rowSums(otu_seqtab[,2:(ncol(otu_seqtab)-1)]), decreasing = T),]

otu_newnames = paste0("OTU_",1:nrow(otu_seqtab))
otu_oldnames = otu_seqtab$otu
otu_seqtab$otu = otu_newnames

map_table_otu = data.frame(otu_newnames = otu_newnames, otu_oldnames = otu_oldnames)

# Extract esv table
results_clustering_asv = lapply(results_clustering, FUN = function(x){x[[2]]})   
esv_seqtab = do.call(rbind,results_clustering_asv)
esv_seqtab = esv_seqtab[order(rowSums(esv_seqtab[,2:(ncol(esv_seqtab)-2)]), decreasing = T),]
esv_seqtab = merge(esv_seqtab, map_table_otu, by.x = "otu", by.y = "otu_oldnames")

esv_seqtab_new = esv_seqtab[,colnames(esv_seqtab) %in% basename(sample_table$filt_R1)]

colnames(esv_seqtab_new) <- gsub("-", 
                             ".", 
                             sample_table$sample[match(colnames(esv_seqtab_new), 
                                                       basename(sample_table$filt_R1))], 
                             fixed = TRUE)

esv_seqtab_new$otu = esv_seqtab$otu
esv_seqtab_new$sequence = esv_seqtab$sequence

# Write ASV table as file
write.table(esv_seqtab,
            file.path(OUTPUT_FOLDER,paste0(RUN_NAME,"_ASV.tsv")), 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

# Save the OTU table as backup
saveRDS(otu_seqtab, file.path(backup_path, "savepoint2_seqtable.rds"))
# Read the sequence table (only a part of the script can be run from here)
otu_seqtab <- readRDS(file.path(backup_path, "savepoint2_seqtable.rds"))



# Saving OTU sequences to fasta file
seqs_otu<-otu_seqtab$seq
names_otu<-otu_seqtab$otu
dna <- DNAStringSet(seqs_otu)
names(dna) = names_otu
Biostrings::writeXStringSet(dna, file.path(OUTPUT_FOLDER,paste0(RUN_NAME,"_OTU.fasta")))

# Classify OTUs
otus_classification = classify_blast(ape::read.FASTA(file.path(OUTPUT_FOLDER,paste0(RUN_NAME,"_OTU.fasta"))), DATABASE_PATH)

classified_otu_table = cbind(otu_seqtab, otus_classification)

col_newnames =   gsub("-", 
                             ".", 
                             sample_table$sample[match(colnames(classified_otu_table), 
                                                       basename(sample_table$filt_R1))], 
                             fixed = TRUE)

colnames(classified_otu_table)[!is.na(col_newnames )] = col_newnames[!is.na(col_newnames )]

write.table(classified_otu_table,
            file.path(OUTPUT_FOLDER,paste0(RUN_NAME,"_classified_otutab.tsv")), 
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")
