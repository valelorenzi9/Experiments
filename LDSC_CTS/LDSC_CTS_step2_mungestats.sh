#!/usr/bin/env bash

############################################################
# Convert GWAS summary statistics file to .sumstats format #
############################################################

while getopts w:g: flag
do
    case "${flag}" in
        w) workdir=${OPTARG};;
        g) gwas=${OPTARG};;
    esac
done
echo -e "Workdir: $workdir\n";
echo -e "GWAS of interest: $gwas\n"

filename=$gwas
name="${filename%%.*}"

# Activate ldsc virtual environment 
#set -e 
#source /home/jovyan/my-conda-envs/ldsc/bin/activate 

# Convert GWAS summary statistics to .sumstats format by using munge_sumstats.py script and keeping only SNPs that are listed in the HapMap3 project
# ! Though not mentioned in the tutorial, it's important to reduce the chunk size otherwise it takes ages to run: by default chunksize = 5000000
python /home/jovyan/ldsc/munge_sumstats.py \
--sumstats $workdir/$gwas \
--merge-alleles $workdir/w_hm3.snplist \
--out $name \
--chunksize 500000
