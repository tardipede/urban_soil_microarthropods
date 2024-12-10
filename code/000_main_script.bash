#!/usr/bin/bash

# Process raw reads into ASV tables
Rscript code/001_coi_reads_process.R       --input "./data"  \
                                           --intermediates "./intermediates" \
                                           --output "./results" \
                                           --run_name "COI_soil" \
                                           --run_regex "(ID[0-9]{4})" \
                                           --sample_regex "((SIN)[0-9]{3})" \
                                           --database_file "./databases/db_all.fas" \
                                           --forward_primer "GGWACWGGWTGAACWGTWTAYCCYCC" \
                                           --reverse_primer "TANACYTCNGGRTGNCCRAARAAYCA" \
                                           --min_amplicon_size 290 \
                                           --max_amplicon_size 350 \
                                           --amplicon_size_step 1 \
                                           --cluster_by_size FALSE

# Primers ref: https://doi.org/10.1186/1742-9994-10-34
# Inosines in rev primer (TAIACYTCIGGRTGICCRAARAAYCA) were replaced with Ns


