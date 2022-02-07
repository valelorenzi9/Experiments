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

import seaborn as sns

# Function to plot QC metrics on a per-sample basis 
def qc_plots_sample(adata, sample, outdir):
    plt.tight_layout()
    
    # Set multi-panel figure
    f, axs = plt.subplots(2, 3, figsize=(14, 9))

    # Violin plot of n_genes, n_counts and percent_mito
    v1 = sc.pl.violin(adata, ['n_genes'], show = False, ax=axs[0][0])
    v1.set_title('n_genes violin plot')
    v2 = sc.pl.violin(adata, ['n_counts'], show = False, ax=axs[0][1])
    v2.set_title('n_counts violin plot')
    v3 = sc.pl.violin(adata, ['percent_mito'], show = False, ax=axs[0][2])
    v3.set_title('percent_mito violin plot')

    # Scatterplots of n_counts vs n_genes and n_counts vs percent_mito
    sc.pl.scatter(adata, x='n_counts', y='percent_mito', show = False, title = 'n_counts vs percent_mito scatterplot', ax=axs[1][0])
    sc.pl.scatter(adata, x='n_counts', y='n_genes', show = False, title = 'n_counts vs n_genes scatterplot', ax=axs[1][1])
    
    # Histograms of n_genes and percent_mito
    sns.histplot(adata.obs['n_genes'], bins = 100, kde=True, ax=axs[1][2])
    plt.axvline(1500, linestyle = '--', color = 'red')
    plt.axvline(8000, linestyle = '--', color = 'green')
    plt.title('n_genes histogram')
    
    f.tight_layout()
    f.savefig(outdir + sample + '_qc_plots.pdf')

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

# Benjamini-Hochberg and Bonferroni FDR helper functions.
def bh(pvalues):
    """
    Computes the Benjamini-Hochberg FDR correction.

    Input:
        * pvals - vector of p-values to correct
    """
    pvalues = np.array(pvalues)
    n = int(pvalues.shape[0])
    new_pvalues = np.empty(n)
    values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
    values.sort()
    values.reverse()
    new_values = []
    for i, vals in enumerate(values):
        rank = n - i
        pvalue, index = vals
        new_values.append((n/rank) * pvalue)
    for i in range(0, int(n)-1):
        if new_values[i] < new_values[i+1]:
            new_values[i+1] = new_values[i]
    for i, vals in enumerate(values):
        pvalue, index = vals
        new_pvalues[index] = new_values[i]
    return new_pvalues

def bonf(pvalues):
    """
    Computes the Bonferroni FDR correction.

    Input:
        * pvals - vector of p-values to correct
    """
    new_pvalues = np.array(pvalues) * len(pvalues)
    new_pvalues[new_pvalues>1] = 1
    return new_pvalues


