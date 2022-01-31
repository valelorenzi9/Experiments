#!/usr/bin/env bash

####################################
# Liftover peaks from hg38 to hg19 #
####################################

# Assign script argument to variable
while getopts o: flag
do
    case "${flag}" in
        o) outdir=${OPTARG};;
    esac
done
echo "Outdir: $outdir";

# If liftOver software is not installed in outdir, download it and give it permissions to be executable
if [ ! -f $outdir/liftOver]
then
    echo "LiftOver script not found in specified outdir, it will be downloaded"
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver -P $outdir
    chmod +x $outdir/liftOver
else
    echo "LiftOver script found in specified outdir, it will NOT be re-downloaded"
fi

if [ ! -f $outdir/hg38ToHg19.over.chain.gz]
then
    echo "hg38ToHg19.over.chain.gz not found in specified outdir, it will be downloaded"
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver -P $outdir
else
    echo "hg38ToHg19.over.chain.gz found in specified outdir, it will NOT be re-downloaded" 
fi

# SNPs files are in hg19 genome assembly, while our peaks are in hg38

for f in $outdir/*.bed; do
    x = basename f .bed	
    echo "Processing ${x}"
    $outdir/liftOver $outdir/$x.bed  \
        $outdir/hg38ToHg19.over.chain.gz \
        $outdir/$x.hg19.bed \
        $outdir/$x.unmapp.bed
done

echo "All bed files coordinates have been lifted from hg38 to hg19"

