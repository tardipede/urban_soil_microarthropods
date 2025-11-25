# Urban_soil_microarthropods


## 1) Place the raw reads files in the data folder
Download the raw reads file from the NCBI project 

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
