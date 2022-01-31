#!/usr/bin/env python3

####################
# Import libraries #
####################

# If libraries are not yet installed, install them 

import pip

def import_or_install(package):
    try:
        __import__(package)
    except ImportError:
        pip.main(['install', package])   

required_packages = ['numpy', 'pandas', 'random', 'sys', 'os', 'datetime', 'argparse', 'warnings', 'pybedtools']
for p in required_packages:
    import_or_install(p) 

##################################
# Define command line parameters #
##################################

parser = argparse.ArgumentParser()

# outdir
parser.add_argument('-o', '--outdir', nargs='+',  type=str, help='Path where the bed files containing differentially accessible peaks are saved.', required = True)

# from
parser.add_argument('-p', '--da-peaks', type=str, nargs='+', help='Path to the .csv file containing the differentially accessible peaks between cell types.',required=True)

# to
parser.add_argument('-c','--cutoff', type=float, nargs='+', help='Adjusted p-value cutoff that used to determine how many peaks per cell type are kept.',required=True)

args = parser.parse_args()

####################
# Print parameters #
####################

print("Launching script with arguments:/n {}".format(args))

##################
# Load dataframe #
##################

# Load .csv file containing the differentially accessible peaks per cell type computed with Seurat 
da_peaks = pd.read_csv(args.da_peaks[0], index_col = 0)
print("Column names: ".format(da_peaks.columns))

# Manipulate dataframe to add columns which are needed to convert to bed files
def make_beds(da_peaks):
    peaks = da_peaks.index.to_list()
    chrom = [i.split["-"][0] for i in peaks]
    chromStart = [int(i.split("-")[1]) for i in peaks] 
    chromEnd = [int(i.split("-")[2]) for i in peaks] 
    bed_dict = {'chrom' : chrom, 'chromStart': chromStart, 'chromEnd' : chromEnd}
    bed_df = pd.DataFrame(bed_dict)
    return bed_df

bed_df_list = []
cell_types = np.unique(da_peaks['cluster'])

# Select peaks per cell type that exceed threshold adjusted p-value 
def select_peaks(da_peaks, celltype):
    # Subset to peaks per cell type 
    subset_peaks = da_peaks[da_peaks['cluster'] == celltype]
    # Thresholding 
    subset_peaks = subset_peaks[subset_peaks['p_val_adj'] < args.cutoff]
    return subset_peaks
    
for c in cell_types:
    peaks = subset_peaks(da_peaks, c)
    bed = make_beds(peaks)
    bed_df_list.append(bed)

# Save bed files per cell type to outdir 
for idx in range(len(cell_types)):
    cell_type = cell_types[idx]
    print("Saving bed files for cell type: {}".format(cell_type))
    bed_df = bed_df_list[idx]
    pybedtools.BedTool.from_dataframe(bed_df, outfile = outdir + cell_type + ".bed", index = False, sep = '\t') 


