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

    # Compute a neighborhood graph of observations
    sc.pp.neighbors(bdata, n_pcs=20)
    # Embed the neighborhood graph using UMAP
    sc.tl.umap(bdata)
    # Cluster GENES into subgroups using leiden: resolution < 1 to find less clusters
    sc.tl.leiden(bdata, resolution=0.5)

    # Locate ccs cluster
    bdata.obs['known_cyclers'] = [i in ['CDK1','MKI67','CCNB2','PCNA'] for i in bdata.obs_names]
    bdata.obs['known_cyclers'] = bdata.obs['known_cyclers'].astype(int)
    sc.pl.umap(bdata, color=['known_cyclers', 'leiden'], color_map='OrRd')
    print(bdata.obs.loc[[i in ['CDK1','MKI67','CCNB2','PCNA'] for i in bdata.obs_names],'leiden'])
    
    if 'CDK1' in bdata.obs_names:
        ccgs_cl = bdata.obs.loc['CDK1',['leiden']][0]
        print("Cell cycle genes cluster is "+ccgs_cl)
    
         # Add unstructured dict-like annotation for ccgs
        adata.uns['ccgs'] = list(bdata.obs[bdata.obs['leiden']==ccgs_cl].index)
    
        # Remove cc genes
        print('Total number of genes before ccg filter: {:d}'.format(adata.n_vars))
        adata = adata[:,[i not in adata.uns['ccgs'] for i in adata.var_names]]
        print('Total number of genes after ccg filter: {:d}'.format(adata.n_vars))
    else: 
        print("WARNING: CDK1 not present in bdata.obs_names, so not removing cell cycle genes")

    return adata

# Function to normalize and log-transform 
def normalize_log_transform(adata):
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    return adata

# Function to filter HVGs and compute PCA on them 
def hvgs_pca_umap(adata):
    bdata = adata.copy()
    sc.pp.highly_variable_genes(bdata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    for col in ['highly_variable','means', 'dispersions', 'dispersions_norm']:
        adata.var[col] = bdata.var[col]
    bdata = bdata[:, bdata.var['highly_variable']]
    sc.pp.scale(bdata, max_value=10)
    sc.tl.pca(bdata, svd_solver='arpack', n_comps=50)
    #fill NaNs with False so that subsetting to HVGs is possible
    adata.var['highly_variable'].fillna(value=False, inplace=True)
    adata.obsm['X_pca'] = bdata.obsm['X_pca'].copy()
    adata.uns['pca'] = bdata.uns['pca'].copy()
    adata.varm['PCs'] = np.zeros(shape=(adata.n_vars, 50))
    adata.varm['PCs'][adata.var['highly_variable']] = bdata.varm['PCs']
    sc.pl.pca_variance_ratio(adata, log=True)

    # Decide number of PCs to compute neighbourhood graph 
    n_pcs = int(input('Desired number of PCs to consider for neighbourhood graph:'))
    sc.pp.neighbors(adata, n_pcs = n_pcs)
    sc.tl.umap(adata)
    return adata

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

# Function to load doublet scores computed with scrublet 
def load_scrublet(adata, samples, scrublet_dir): 
    scorenames = ['scrublet_score','scrublet_cluster_score','zscore','bh_pval','bonf_pval']

    scrdf = []
    for sample in samples:
        scrdf.append(pd.read_csv(scrublet_dir + sample + '.csv', header=0, index_col=0))
    scrdf = pd.concat(scrdf)
    scrdf.index = [i.replace('-1', '') for i in scrdf.index]
    for score in scorenames:
        adata.obs[score] = scrdf[score]
    adata.obs['is_doublet'] = adata.obs['bonf_pval'] < 0.01

    # doublets %
    print("Percentage of doublets in dataset: {}".format(adata.obs['is_doublet'].sum() / adata.shape[0]))

    # Fix is_doublet data type to enable plotting correctly 
    adata.obs['is_doublet'] = adata.obs['is_doublet'].astype(int)
