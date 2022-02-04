#!/usr/bin/env python3

import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys
import scipy
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
from matplotlib.backends.backend_pdf import PdfPages

# Function to make diagnostics plots to evaluate the quality of the sample
def qc_plots_sample(adata, sample, outdir):
    plt.tight_layout()

    # Violin plot of n_genes, n_counts and percent_mito
    fig1 = sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'], multi_panel=True, jitter = 0.4, show = False)
    plt.savefig(outdir + sample + '_qc_violin.pdf', dpi = 300)

    # Scatterplots of n_counts vs n_genes and n_counts vs percent_mito
    fig2 = sc.pl.scatter(adata, x='n_counts', y='percent_mito', show = False)
    plt.savefig(outdir + sample + '_qc_scatter1.pdf', dpi = 300)
    fig3 = sc.pl.scatter(adata, x='n_counts', y='n_genes', show = False)
    plt.savefig(outdir + sample + '_qc_scatter2.pdf', dpi = 300)
    
    
    # Histograms of n_genes and percent_mito
    histogram = plt.figure()
    axis1 = histogram.add_subplot(211)
    axis1.hist(adata.obs['n_genes'], bins = 100)
    axis1.axvline(1000, linestyle = '--', color = 'red')
    
    axis2 = histogram.add_subplot(212)
    axis2.hist(adata.obs['percent_mito'], bins = 100, cumulative=True)
    axis2.axvline(0.1, linestyle = '--', color = 'red')
    axis2.axvline(0.2, linestyle = '--', color = 'darkred')
    axis2.axhline(adata.n_obs*0.99, linestyle = '-', color = 'green')

    plt.tight_layout()
    histogram.savefig(outdir + sample + '_qc_histograms.pdf', dpi = 300)


# Function to run gene centric analysis to identify genes that behave like known cell cycle genes 
def per_gene_analysis(adata):
    bdata = adata.copy()
    # Normalize total counts per cell
    sc.pp.normalize_per_cell(bdata, counts_per_cell_after=1e4)
    # Logarithmize the data matrix
    sc.pp.log1p(bdata)

    # Extract highly variable genes
    sc.pp.highly_variable_genes(bdata)
    highly_variable_genes = bdata.var["highly_variable"]
    bdata = bdata[:, highly_variable_genes]

    # Traspose matrix for a GENE-centered analysis
    bdata = bdata.copy().T

    # Scale data to unit variance and zero mean
    sc.pp.scale(bdata, max_value=10)

    # Scatter plot in PCA coordinates
    sc.tl.pca(bdata)
    bdata.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat
    # Plot the variance ratio
    sc.pl.pca_variance_ratio(bdata, log=True)

    return bdata
