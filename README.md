# Urban_soil_microarthropods


## 1) Create a folder name data and place there the raw reads files
Download the raw reads file from the NCBI project PRJNA1369507
The runs names should be downloaded with their original names, for example:  
⋅⋅⋅⋅* 145055_ID3302_557-IT-210-concrete-pl6-E10-COIA_S1_L001_R1_001.fastq.gz  
⋅⋅⋅⋅* 145055_ID3302_557-IT-210-concrete-pl6-E10-COIA_S1_L001_R2_001.fastq.gz  
  
Otherwise change accordingly the *--run_regex* and  *--sample_regex* options in the *./code/000_main_script.bash* file. 

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
Rscript db_format.R
cd ..
```

## 7) Run main script
```
bash ./code/000_main_script.bash
```

## 8) Calculate alpha diversity
Move the "COI_soil_classified_otutab.tsv" output file to the alpha_diversity_calculations folder and run the R script in it
