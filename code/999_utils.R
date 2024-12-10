#!/usr/bin/env Rscript

# Filter sequences based on absence of stop codons (set on genetic code 5)
translation_filtering = function(seqtable, frame = 2){
  require(Biostrings)
  
  GENETIC_CODE <- Biostrings::getGeneticCode("5")  # invertebrate mitochondrial
  
  if(frame == 1){  AMINO<-lapply(lapply( sub('', '', colnames(seqtable)),  DNAString  ), function(x){Biostrings::translate(x, GENETIC_CODE)})}
  if(frame == 2){  AMINO<-lapply(lapply( sub('.', '', colnames(seqtable)),  DNAString  ), function(x){Biostrings::translate(x, GENETIC_CODE)})}
  if(frame == 3){  AMINO<-lapply(lapply( sub('..', '', colnames(seqtable)),  DNAString  ), function(x){Biostrings::translate(x, GENETIC_CODE)})}
  
  AMINO_char<-lapply(AMINO, as.character)
  with_stop_codons<- unlist(lapply(AMINO_char, function(x){grepl("\\*",x)}))
  
  # Frequency of sequence variants with stop codons:
  #mean(with_stop_codons)
  
  # Filtering
  seqtable_filtered<- seqtable[,!with_stop_codons]
  
  return(seqtable_filtered)}

# Cluster the sequences
cluster_seqs = function(seqtable, nproc, cutoff, align = F){
  asv_sequences <- colnames(seqtable)
  sample_names <- rownames(seqtable)
  
  aln <- Biostrings::DNAStringSet(asv_sequences)
  
  if(align == T){aln = DECIPHER::AlignSeqs(aln, processors = nproc)}

  d <- DECIPHER::DistanceMatrix(aln, processors = nproc,type="dist")  # type="dist"  saves memory

  
  clusters<- DECIPHER::TreeLine(
    myDistMatrix = d, 
    method = "complete",
    cutoff = cutoff,
    type = "clusters",
    processors = nproc,
    showPlot = FALSE,
    verbose = TRUE)
  
  clusters <- clusters %>% add_column(sequence = asv_sequences)
  rm(d)  # free up memory
  
  ### First (most common) sequence as representative for OTU
  otu_sequences<-clusters %>% group_by(cluster) %>% summarise(seq = sequence[1])
  
  merged_seqtab <- seqtable %>%
    # setup: turn seqtab into a tibble with rows = ASVs and columns = samples
    t %>%
    as_tibble(rownames = "sequence") %>%
    # add the cluster information
    left_join(clusters, by = "sequence") %>%
    # merge ESVs in the same cluster, summing abundances within samples
    group_by(cluster) %>%
    summarize_at(vars(-sequence), sum) %>%
    left_join(otu_sequences, by = "cluster")%>%
    # Set new taxa names to OTU<cluster #> 
    mutate(cluster = paste0("OTU", cluster)) %>%
    column_to_rownames("cluster")%>%   
    # move seq column to beginning
    select(seq, everything())

  
  merged_seqtab<-merged_seqtab[order(-rowSums(merged_seqtab[-c(1)])),]  	
  old_otu = rownames(merged_seqtab)
  merged_seqtab$otu<-paste0("OTU",1:nrow(merged_seqtab))
  rownames(merged_seqtab)<-merged_seqtab$otu
  otu_conv_table = data.frame(old_otu = old_otu, new_otu =  merged_seqtab$otu)
  clusters = clusters %>% mutate(otu = paste0("OTU", cluster))
  
  asv_otu_table = seqtable %>%
    # setup: turn seqtab into a tibble with rows = ASVs and columns = samples
    t %>%
    as_tibble(rownames = "sequence") %>%
    # add the cluster information
    left_join(clusters, by = "sequence")
  
  asv_otu_table = merge(asv_otu_table, otu_conv_table, by.x = "otu", by.y = "old_otu")
  asv_otu_table = asv_otu_table[,2:ncol(asv_otu_table)]
  colnames(asv_otu_table)[ncol(asv_otu_table)] = "otu"
  
  results = list(otu_tab = merged_seqtab, esv_tab = asv_otu_table)
  
  return(results)}


######## CLASSIFY COI OTUS WITH BLAST ################
# Database needs to be in fasta file with the names of the sequences as:
# Kingdom;Phylum;Class;Order;Family;Genus;Species

classify_blast = function(query, db_path, remove_temp_dir = T, ...){
  # create directory for db and results
  dir.create("temp_blast")
  query_path = "./temp_blast/query.fas"
  ape::write.FASTA(query, query_path)
  
  # make blast database
  system(paste0("makeblastdb -in ", db_path, " -dbtype nucl -out ./temp_blast/mydb -title mydb"))
  
  # Run blastn
  system(paste0("blastn -db ./temp_blast/mydb -query ",
                query_path,
                " -perc_identity 75 -out ./temp_blast/out.txt -evalue 1e-6 -outfmt \"6 qseqid sseqid pident length evalue bitscore\" -qcov_hsp_perc 50 -max_target_seqs 25"))
  
  
  # Load blast results
  blast_results = read.table("./temp_blast/out.txt", header = F, sep = "\t")
  colnames(blast_results) = c("qseqid","sseqid","pident","length","evalue","bitscore")
  
  # classify OTUs based on BLAST file
  otus_id = names(query)
  class = lapply(otus_id, FUN = classify_otu_blast, blast_results = blast_results)
  class = do.call(rbind, class)
  class = cbind(otus_id, class)
  class = data.frame(class)
  colnames(class) = c("OTU","n_BLAST_matches","avg_pident","Kingdom","Phylum","Class","Order","Family","Genus","Species",
                      "sc_Kingdom","sc_Phylum","sc_Class","sc_Order","sc_Family","sc_Genus","sc_Species","tax_level")
  
  if(remove_temp_dir == TRUE){unlink("temp_blast", recursive  = T)}
  
  return(class)}

classify_vector_tax = function(x){
  if(length(unique(x))==1){results = c(x[1],100)}
  if(length(unique(x))!=1){
    tabtaxa = table(x)
    tabtaxa = tabtaxa/sum(tabtaxa)*100
    if(max(tabtaxa)<70){results = c("","")}
    if(max(tabtaxa)>=70){results = c(names(tabtaxa)[which.max(tabtaxa)], max(tabtaxa))}
  }
  return(results)
}

classify_otu_blast = function(OTU, blast_results, diff_high_match = 1, min_dist = c(96,94,90,85)){
  
  
  blas_results_temp = subset(blast_results, qseqid == OTU)
  
  if(nrow(blas_results_temp) == 0){
    empty_results = c(n_matches = 0, avg_pident=NA,
                      rep("",14),tax_level = 0)
    return(empty_results)
  }
  
  if(nrow(blas_results_temp) != 0){
    
    blas_results_temp = blas_results_temp[order(blas_results_temp$pident, decreasing =T),]
    max_pid = blas_results_temp[1,"pident"]
    blas_results_temp = subset(blas_results_temp, pident >= (max_pid-diff_high_match))
    
    avg_pident = mean(blas_results_temp$pident)
    n_matches = nrow(blas_results_temp)
    
    ###################
    
    
    tax_table = blas_results_temp$sseqid
    
    tax_table = do.call(rbind,strsplit(as.character(tax_table), ";")) # Format taxonomy in names as table
    tax_table = data.frame(tax_table)
    colnames(tax_table) = c("kingdom","phylum","class","order","family","genus","species")
    
    tax_class = apply(tax_table, MARGIN = 2, FUN = classify_vector_tax)
    results_string = c(n_matches = n_matches, avg_pident=avg_pident,tax_class[1,],tax_class[2,])
    
    if(avg_pident < min_dist[1]){results_string[c(9,16)] = c("","")} # relaxed from 97 to 96 as tardigrades have generally higher COI intraspec. diversity
    if(avg_pident < min_dist[2]){results_string[c(8,15)] = c("","")} # relaxed from 95 as above
    if(avg_pident < min_dist[3]){results_string[c(7,14)] = c("","")}
    if(avg_pident < min_dist[4]){results_string[c(6,13)] = c("","")}
    
    tax_level = sum(results_string[3:9] != "")
    results_string = c(results_string,tax_level = tax_level)
    
    return(results_string)
    
  }
}


##### Subtract reads in blank from all the samples
# Following https://doi.org/10.1002/edn3.372
subtract_blank <- function(phyloseq_object, blank_name) {
  require(phyloseq)


  otu_tab_temp <- data.frame(otu_table(phyloseq_object))
  otu_tab_temp <- sweep(otu_tab_temp, 
                        MARGIN = 1, 
                        STATS = otu_tab_temp[, blank_name], 
                        FUN = "-")
  otu_tab_temp[otu_tab_temp < 0] <- 0
  new_otu_table <- phyloseq::otu_table(otu_tab_temp, 
                                       taxa_are_rows = TRUE)

  phyloseq_object_new <- phyloseq(new_otu_table, 
                                  sample_data(phyloseq_object), 
                                  tax_table(phyloseq_object), 
                                  refseq(phyloseq_object))

  return(phyloseq_object_new)
}

##### Save RDS and return object (to use with magrittr %>%)
saveRDS_pipe <- function(object_to_save, filename) {
  saveRDS(object_to_save, filename)
  return(object_to_save)
}
