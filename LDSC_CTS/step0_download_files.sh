#!/usr/bin/env bash

################################################################
# Download necessary files for running cell type-specific LDSC #
################################################################

# Assign script argument to variable 

while getopts w: flag
do
    case "${flag}" in
        w) workdir=${OPTARG};;
    esac
done
echo -e "Workdir: $workdir\n";

# Download SNPs that are listed in the HapMap3 project, if they haven't been downloaded already 

if [ ! -f $workdir/w_hm3.snplist.bz2 ]
then
    echo -e "w_hm3.snplist.bz2 not found in specified directory, it will be downloaded\n"
    wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2 -P $workdir
    bunzip2 w_hm3.snplist.bz2
else
    echo -e "w_hm3.snplist.bz2 found in specified directory, it will NOT be re-downloaded\n"
fi

# If liftOver software is not installed, download it and give it permissions to be executable

if [ ! -f $workdir/liftOver ]
then
    echo -e "LiftOver script not found in specified directory, it will be downloaded\n"
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver -P $workdir
    chmod +x $workdir/liftOver
else
    echo -e "LiftOver script found in specified directory, it will NOT be re-downloaded\n"
fi

if [ ! -f $workdir/hg38ToHg19.over.chain.gz ]
then
    echo -e "hg38ToHg19.over.chain.gz not found in specified directory, it will be downloaded\n"
    wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz -P $workdir
else
    echo -e "hg38ToHg19.over.chain.gz found in specified directory, it will NOT be re-downloaded\n" 
fi 

if [ ! -f $workdir/1000G_Phase3_plinkfiles.tgz ]
then
    echo -e "1000G_Phase3_plinkfiles.tgz not found in specified directory, it will be downloaded\n"
    wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_plinkfiles.tgz -P $workdir
else
    echo -e "1000G_Phase3_plinkfiles.tgz found in specified directory, it will NOT be re-downloaded\n"
fi

if [ ! -f $workdir/hapmap3_snps.tgz ]
then
    echo -e "hapmap3_snps.tgz not found in specified directory, it will be downloaded\n"
    wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/hapmap3_snps.tgz -P $workdir
else
    echo -e "hapmap3_snps.tgz found in specified directory, it will NOT be re-downloaded\n"
fi

if [ ! -f $workdir/1000G_Phase3_baseline_ldscores.tgz ]
then
    echo -e "1000G_Phase3_baseline_ldscores.tgz not found in specified directory, it will be downloaded\n"
    wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_ldscores.tgz -P $workdir
    tar -xvzf 1000G_Phase3_baseline_ldscores.tgz
else
    echo -e "1000G_Phase3_baseline_ldscores.tgz found in specified directory, it will NOT be re-downloaded\n"
fi

if [ ! -f $workdir/weights_hm3_no_hla.tgz ]
then
    echo -e "weights_hm3_no_hla.tgz not found in specified directory, it will be downloaded\n"
    wget https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz -P $workdir
    tar -xvzf weights_hm3_no_hla.tgz
else
    echo -e "weights_hm3_no_hla.tgz found in specified directory, it will NOT be re-downloaded\n"
fi
