###################
# Import packages #
###################

# Required packages
packages <- c("Seurat", "SeuratDisk", "argparse")

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
                    help = "Path to Seurat object you wish to convert to Anndata")

parser$add_argument("--object_name",
                    type="character",
                    help = "Name of the saved object")

args <- parser$parse_args()

seurat_object <- args$seurat_object

object_name <- args$object_name

####################################
# Convert Seurat object to Anndata #
####################################

# Load Seurat object
SeuratObject <- readRDS(file = seurat_object)

# Save intermediate format (h5Seurat)
SaveH5Seurat(SeuratObject, filename = paste0(object_name, ".h5Seurat"))

# Convert and save to h5ad file 
Convert(paste0(object_name, ".h5Seurat"), dest = "h5ad")

print("Seurat object was successfully converted to Anndata")
