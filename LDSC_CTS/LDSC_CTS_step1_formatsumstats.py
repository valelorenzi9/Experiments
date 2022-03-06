#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()

# workdir
parser.add_argument('-w', '--workdir', nargs='+',  type=str, help='Path to working directory.', required = True)

# gwas
parser.add_argument('-g', '--gwas', type=str, nargs='+', help='Filename of the GWAS summary statistics of interest.',required=True)

args = parser.parse_args()

print("Launching formatting summary statistics script with arguments:/n {}".format(args))

import pandas as pd
sumstats = pd.read_csv(args.workdir[0] + args.gwas[0], compression = 'gzip', 
                      sep = '\t')
sumstats = sumstats.rename(columns={'variant_id': 'rsid', 'effect_allele_frequency' : 'eaf'})
sumstats.to_csv(args.workdir[0] + args.gwas[0], compression = 'gzip', 
                      sep = '\t')
