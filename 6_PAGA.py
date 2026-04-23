# Import required libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
import sys

# Read annotated epithelial cell data
adata = sc.read_h5ad("Epi_CA_celltype.h5ad")
print(adata)

# Subset data: keep NC, CIN, NECC classes and select specific cell types
adata_subset = adata[adata.obs['class'].isin(['NC', 'CIN', 'NECC'])].copy()
adata_subset = adata_subset[adata_subset.obs['celltype.class'].isin([
    'CA_NE', 'Epi_Squamous', 'Epi_Columnar', 'Epi_Glandular',
    'CA_Columnar', 'Epi_Cililated', 'CA_Glandular', 'CA_Ciliated', 'CA_Squamous'
])].copy()
adata = adata_subset

# Configure plotting settings
sc.settings.verbosity = 3
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(7, 7), facecolor='white')

# Normalization and dimensionality reduction
sc.pp.normalize_total(adata, target_sum=1e4)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)

# Diffusion maps and graph-based layout
sc.tl.diffmap(adata)
sc.pp.neighbors(adata, n_neighbors=30, use_rep='X_diffmap')
sc.tl.draw_graph(adata)

# Define new categorical groups based on patient characteristics
adata.obs['new_group'] = np.select(
    [
        adata.obs['Patient'].isin(['CIN1', 'CIN2', 'CIN4']),   # CIN with chr1q21 loss
        adata.obs['Patient'].eq('CIN3'),                       # CIN with chr1q21 gain
        adata.obs['class'].eq('NC'),                           # Normal cervix
        adata.obs['class'].eq('NECC')                          # Neuroendocrine cervical carcinoma
    ],
    [
        'CIN_chr1q21_loss',
        'CIN_chr1q21_gain',
        'NC',
        'NECC'
    ],
    default='Other'
)

# Convert to categorical and check value counts
adata.obs['new_group'] = adata.obs['new_group'].astype('category')
print(adata.obs['new_group'].value_counts())

# PAGA (Partition-based graph abstraction) analysis
sc.tl.paga(adata_subset, groups='new_group')

# Plot PAGA graph
sc.pl.paga(adata_subset,
           threshold=0.25,
           node_size_scale=3,
           edge_width_scale=0.15,
           solid_edges="connectivities",
           fontsize=5,
           labels=None,
           show=True)
