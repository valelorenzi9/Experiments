###################
# Import packages #
###################

# Required packages
packages <- c("Seurat", "Signac", "argparse")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


##################################
# Define command line parameters #
##################################

parser <- ArgumentParser()

parser$add_argument("--seurat_object", type="character",
                    help = "Path to Seurat object with chromatin accessibility information.")

parser$add_argument("--labels", type = "character", 
		    help = "Metadata column that contains the cell type labels of interest")

parser$add_argument("--outdir", type = "character",
                    help = "Path to where you want to store the differentially accessible regions.")

args <- parser$parse_args()

seurat_object <- args$seurat_object

labels <- args$labels

outdir <- args$outdir

print("Selected arguments: ")
print(paste0("Seurat Object: ", seurat_object))
print(paste0("Cell type labels column: ", labels))
print(paste0("Outdir: ", outdir))

###########################################################
# Find differentially accessible peaks between cell types #
###########################################################

# Load Seurat object
SeuratObject <- readRDS(file = seurat_object)
print("Seurat object loaded")

# Set default assay to peaks (in case it is not already) 
DefaultAssay(SeuratObject) <- 'peaks' 

# Set idents to the cell type labels (in case they were not already)
Idents(SeuratObject) <- SeuratObject@meta.data$labels 

# Use Seurat's approach with logistic regression for differential accessibility and add the total number of fragments as a latent variable to mitigate the effect of differential sequencing depth
da_peaks <- FindAllMarkers(object = SeuratObject, test.use = 'LR', only.pos = F, random.seed = 1, min.pct = 0.05, latent.vars = 'peak_region_fragments') 
print("Computed differentially accessible peaks between clusters") 

# Save differentially accessible peaks in a .csv file 
write.csv(da_peaks, file = paste0(outdir, "da_peaks_LR.csv"))
print("Saved differentially accessible peaks")

