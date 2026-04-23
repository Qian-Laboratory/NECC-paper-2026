# ============================================
# Single-cell RNA-seq Data Integration and Analysis
# ============================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(future)

# Set parallelization options
options(future.globals.maxSize = 20 * 1000 * 1024^2)

# Split by sample
merge <- readRDS("integrated_final.rds")
DefaultAssay(merge) <- 'RNA'
merge_split <- SplitObject(merge, split.by = "sample")

# Process each sample: regress out technical variables
features <- SelectIntegrationFeatures(object.list = merge_split,
                                      nfeatures = 2200)
merge_split <- lapply(X = merge_split, FUN = function(x){
    x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
    x <- ScaleData(x, vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
                   features = features)
    x <- RunPCA(x, npcs = 30, features = features)
    return(x)
})

# Find integration anchors (RPCA method)
anchors <- FindIntegrationAnchors(merge_split,
                                  k.anchor = 5,
                                  dims = 1:30,
                                  k.filter = 200,
                                  k.score = 30,
                                  anchor.features = features,
                                  reduction = 'cca')

# Integrate data
merge_integrated <- IntegrateData(anchorset = anchors,
                                dims = 1:30,
                                k.weight = 100)

# Dimensionality reduction and clustering
DefaultAssay(merge_integrated) <- 'integrated'
var.features  <- VariableFeatures(merge_integrated)
merge_integrated <- ScaleData(merge_integrated)
merge_integrated <- RunPCA(merge_integrated, features = var.features, npcs = 30)
merge_integrated <- RunUMAP(merge_integrated, dims = 1:30, umap.method = 'uwot')
merge_integrated <- FindNeighbors(merge_integrated, reduction = "pca", dims = 1:30, force.recalc = TRUE)
merge_integrated <- FindClusters(merge_integrated, algorithm = 1,
                               resolution = c(0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.6, 2, 2.5, 3, 3.5, 4, 4.5, 5))

# Visualization
DimPlot(merge_integrated, reduction = 'umap', group.by = "sample")
DimPlot(merge_integrated, reduction = 'umap', group.by = "origin")
DimPlot(merge_integrated, reduction = 'umap', label = TRUE, group.by = "Phase")
DimPlot(merge_integrated, reduction = 'umap', label = TRUE, group.by = "scrublet_callrate")
DimPlot(merge_integrated, reduction = 'umap', label = TRUE, group.by = "DF.classifications")

# Feature plots
FeaturePlot(merge_integrated, features = c('nFeature_RNA', 'percent.mt', 'G2M.Score', 'S.Score'))
FeaturePlot(merge_integrated, features = c('IFN.Score1', 'Stress.Score1', 'Hypoxia.Score1'))

# Cell type annotation using marker genes
marker_genes <- c("MYH11", "RGS5", "ACTA2", "LUM",           # Mesenchymal
                  "CD24", "EPCAM",                          # Epithelial
                  "KRT6A", "FOXJ1", "CEACAM6", "MSLN", "SCGB3A1", "SBSN",  # Epithelial subsets
                  "SYP",                                    # Neuroendocrine
                  "VWF",                                    # Endothelial
                  "CSF3R",                                  # Granulocyte
                  "LYZ", "CD68",                            # Myeloid
                  "CD3D", "TRDC", "CD8A", "CD4", "FOXP3",   # T cells
                  "NCAM1", "NCR1",                          # NK cells
                  "CD79A", "JCHAIN", "KIT", "MS4A2")        # B cells / Mast cells

DefaultAssay(merge_integrated) <- "RNA"
Idents(merge_integrated) <- "integrated_snn_res.1.2"

DotPlot(merge_integrated, features = marker_genes) +
    scale_color_gradient(low = "#00FFFF", high = "red") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save integrated object
saveRDS(merge_integrated, "merge_celltype_integrated.rds")
