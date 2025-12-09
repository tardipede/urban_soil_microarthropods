# This script calculates the alpha diversity of Arthropods based on the OTU table from metabarcoding
# Before running it move into this folder the "COI_soil_classified_otutab.tsv" output file from the metabarcoding pipeline

library(tidyverse)
library(iNEXT)
library(writexl)

### Load data and metadata
data = read.table("COI_soil_classified_otutab.tsv",
                  sep = "\t",
                  header = TRUE)


### Prepare data
# Remove negs function to remove reads with negative count after subtracting blank
remove_negs = function(x){x[x<0] = 1}

# Extract OTU table and remove blank
otutable = data %>% 
  dplyr::select(colnames(.)[grepl(pattern = "SIN",
                      x = colnames(.),
                      fixed = TRUE)]) %>% # Keep only sample columns
  sweep(MARGIN = 1, STATS = .$SIN000, FUN = "-")  # Subtract blank
otutable[otutable < 0] = 0
rownames(otutable) = data$OTU


# Select Arthropod OTUS with Blast avg similarity > 80%
otutable = otutable[data$avg_pident > 80 & grepl("Arthropoda",data$Phylum),]

# Remove samples with 0 reads
otutable = otutable[,colSums(otutable) != 0]

### Community analysis - alpha diversity

# Rarefaction/Extrapolation: calculate asymptotic richness
alpha_inext = iNEXT(otutable, 
                    q=c(0), 
                    knots = 50,
                    datatype="abundance", 
                    endpoint=min(colSums(otutable)))

# Check the difference between observed and predicted richness
diversity_summary = alpha_inext$AsyEst %>%
  subset(Diversity == "Species richness") %>%
  mutate(tot_reads = colSums(otutable))


writexl::write_xlsx(diversity_summary, "diversity_summary.xlsx")


# Get Arthropoda OTUs taxonomy
a_otus = rownames(otutable)
data_ar = subset(data, data$OTU %in% a_otus)
writexl::write_xlsx(data_ar, "Arthropod_OTUs.xlsx")
table(data_ar$Class)
