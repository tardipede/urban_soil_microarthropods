# Urban soil microarthropods study scripts

Data analysis scripts associated with the manuscript:
Lami F., Zavatta L., Morelli A., Ciurli A., Bazzocchi G.G., Cavani G., Vecchi M. (2026). Complex interactions between local habitat features, landscape factors and seasonality shape soil microarthropod communities in urban green areas. *Urban Ecosystems*.


## 1) Create a folder name data and place there the raw reads files
Download the raw reads file from the NCBI project PRJNA1369507
The runs names should be downloaded with their original names, for example:  
    * Forward reads: 170514_ID3543_378-SIN000-COI-MET-plate4-H2_S289_L001_R1_001.fastq.gz  
    * Reverse reads: 170514_ID3543_378-SIN000-COI-MET-plate4-H2_S289_L001_R2_001.fastq.gz  
  
Otherwise change accordingly the *--run_regex* and  *--sample_regex* options in the *./code/000_main_script.bash* file.  
If you use runs with different filenames, check also if changing the *--forward_reads_regex* and *--reverse_reads_regex* arguments from their default.  

## 2) Create conda env
```
conda env create -f ./conda_envs/environment.yml
```
## 3) Activate conda env
```
conda activate metabarcoding_env
```
## 4) Download and trim database based on primers
```
cd databases
wget https://www.reference-midori.info/download/Databases/GenBank264_2024-12-14/RAW_sp/uniq/MIDORI2_UNIQ_SP_NUC_GB264_CO1_RAW.fasta.gz

cutadapt -g GGWACWGGWTGAACWGTWTAYCCYCC...TGRTTYTTYGGNCAYCCNGARGTNTA \
         --discard-untrimmed \
         --revcomp \
         --report minimal \
         --cores 12 \
         -m 200 \
         -M 400 \
         -e 0.25 \
         -o ./MIDORI2_DB_trimmed.fasta \
         ./MIDORI2_UNIQ_SP_NUC_GB264_CO1_RAW.fasta.gz
```
## 5) Dereplicate based on sequence
```
seqkit rmdup ./MIDORI2_DB_trimmed.fasta --by-seq -o ./MIDORI2_DB_trimmed_derep.fasta
```
## 6) Format taxonomy
```
chmod +x ./*
Rscript db_format.R
cd ..
```

## 7) Run main script
```
chmod +x code/*
bash ./code/000_main_script.bash
```

## 8) Calculate alpha diversity
Move the "COI_soil_classified_otutab.tsv" output file to the alpha_diversity_calculations folder and run the R script in it




[![DOI](https://zenodo.org/badge/901294411.svg)](https://doi.org/10.5281/zenodo.17877587)

