#!/usr/bin/env bash

############################################################
# Convert GWAS summary statistics file to .sumstats format #
############################################################

# Assign script argument to variable 
while getopts g:o: flag
do
    case "${flag}" in
        g) gwas=${OPTARG};;
        o) outdir=${OPTARG};;
    esac
done
echo "GWAS name: $gwas";
echo "Outdir: $outdir";

# Download SNPs that are listed in the HapMap3 project, if they haven't been downloaded already 

if [ ! -f $outdir/w_hm3.snplist.bz2]
then
    echo "List of SNPs in HapMap3 project not found in specified outdir, it will be downloaded"
    wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
    bunzip2 w_hm3.snplist.bz2
else
    echo "List of SNPs in HapMap3 project found in specified outdir, it will NOT be re-downloaded"
fi

# Activate ldsc virtual environment 
set -e 
source /home/jovyan/my-conda-envs/ldsc/bin/activate 

# Convert GWAS summary statistics to .sumstats format by using munge_sumstats.py script and keeping only SNPs that are listed in the HapMap3 project
python /home/jovyan/ldsc/munge_sumstats.py \
--sumstats $outdir/$gwas.sumstats.gz \
--merge-alleles $outdir/w_hm3.snplist \
--out $gwas
