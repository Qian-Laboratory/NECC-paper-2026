# Import necessary libraries
import gc
import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt
import pandas as pd

# Read the AnnData object containing epithelial tumor cells
adata = sc.read_h5ad("Epi_CA_celltype.h5ad")

# Set variable names and normalize
adata.var_names = adata.var.features
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# PCA and neighborhood graph
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)

# UMAP visualization
sc.tl.umap(adata)
sc.pl.umap(adata, color="Sample")

# Add genomic positions from GTF file for CNV analysis
cnv.io.genomic_position_from_gtf(gtf_file="/share/reference/10X/refdata-gex-GRCh38-2020-A/genes/genes.gtf",
                                 adata=adata,
                                 gtf_gene_id='gene_name',
                                 inplace=True)

# Preview genomic position data
adata.var.loc[:, ["features", "chromosome", "start", "end"]].head()

# Extract and optionally save gene positions
data = adata.var.loc[:, ["features", "chromosome", "start", "end"]]
# pd.DataFrame(data).to_csv('/path/to/gene_pos.txt', index=True, mode='a')

# Run inferCNV to infer copy number variations
cnv.tl.infercnv(adata=adata,
                reference_key='class.new',
                reference_cat=['Normal.Cervix'],
                n_jobs=40,
                window_size=100)

# Additional CNV analysis steps
cnv.tl.pca(adata)
cnv.pp.neighbors(adata)
cnv.tl.leiden(adata)
cnv.tl.umap(adata)
cnv.tl.cnv_score(adata)

# Reorder patient categories for visualization
new_categories = pd.Index(["NC1", "NC2", "NC3", "NC4", "NC5", "NC6",
                           "CIN1", "CIN2", "CIN3", "CIN4",
                           "SCC1", "SCC2", "SCC3", "SCC4", "SCC5", "SCC6", "SCC7", "SCC8", "SCC9",
                           "ADC1", "ADC2", "ADC3", "ADC4",
                           "NECC2", "NECC3", "NECC4", "NECC5", "NECC6", "NECC7"])

adata.obs['Patient'] = adata.obs['Patient'].cat.reorder_categories(new_categories, ordered=True)

# Create figure object and set size
fig = plt.figure()
fig.set_size_inches(8, 6)
print(fig.get_size_inches())

# Generate chromosome heatmap grouped by Patient
cnv.pl.chromosome_heatmap(adata, groupby="Patient",
                          figsize=(16, 10),
                          show=False)

# Generate summary chromosome heatmap
cnv.pl.chromosome_heatmap_summary(adata,
                                  figsize=(8, 5),
                                  groupby="Patient",
                                  show=False)
