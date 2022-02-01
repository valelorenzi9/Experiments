#!/usr/bin/env bash

######################################
# Compute annot files from bed files #
######################################

# Assign script argument to variable 
while getopts u:b:o: flag
do
    case "${flag}" in
	u) ldsc_utils_dir=${OPTARG};;
        b) bed_dir=${OPTARG};;
        o) outdir=${OPTARG};;
    esac
done
echo "LDSC utils dir: $ldsc_utils_dir";
echo "Bed files dir: $bed_dir";
echo "Outdir: $outdir";

# Download required files if not already downloaded in the specified ldsc_utils_dir
if [ ! -f $ldsc_utils_dir/1000G_Phase3_plinkfiles.tgz]
then
    echo "1000G_Phase3_plinkfiles.tgz not found in specified ldsc_utils_dir, it will be downloaded"
    wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/1000G_Phase3_plinkfiles.tgz -P $ldsc_utils_dir
else
    echo "1000G_Phase3_plinkfiles.tgz found in specified ldsc_utils_dir, it will NOT be re-downloaded"
fi

if [ ! -f $ldsc_utils_dir/1000G_Phase3_plinkfiles.tgz]
then
    echo "hapmap3_snps.tgz not found in specified ldsc_utils_dir, it will be downloaded"
    wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/hapmap3_snps.tgz -P $ldsc_utils_dir
else
    echo "hapmap3_snps.tgz found in specified ldsc_utils_dir, it will NOT be re-downloaded"
fi

# Run make_annot.py per cell type and per chromosome 

for x in GermCells coelEpi ovarianSurf preGC_I_OSR1 preGC_II_KITLG preGC_III_GJA1 sPAX8 Oi Gi PV Immune Endothelial ; do
        echo "Processing ${x}"
        for chr in {1..22}
        do
                python /home/jovyan/ldsc/make_annot.py \
                        --bed-file $bed_dir/${x}.hg19.bed \
                        --bimfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
                        --annot-file $outdir/${x}.${chr}.annot.gz
        done
done
