# Load required libraries
library(spacexr)
library(Matrix)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(purrr)
library(stringr)

# 1. RCTD for deconv
# Load reference single-cell data
merge_integrated <- readRDS("merge_celltype_integrated.rds")


# Subset paired samples and downsample
i == "NECC1"
sample <- subset(merge_integrated, Patient == i)

# Create Reference object for RCTD
counts <- sample$RNA@counts %>% as.matrix()
meta.df <- sample@meta.data %>%
  rownames_to_column("barcode") %>%
  select(barcode, nCount_RNA, celltype)

cell_types <- meta.df$celltype
names(cell_types) <- meta.df$barcode
nCount_RNA <- meta.df$nCount_RNA
names(nCount_RNA) <- meta.df$barcode

reference <- Reference(counts, cell_types, nCount_RNA)

# Load 10X Visium spatial data
vNECC1 <- Load10X_Spatial(
  "vNECC1/outs",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "vNECC1",
  filter.matrix = TRUE,
  to.upper = FALSE
)

# Extract spot coordinates
spot <- colnames(vNECC1)
coord <- read.csv("/vNECC1/outs/spatial/tissue_positions.csv")
coord <- coord %>% filter(barcode %in% spot) %>%
  select(barcode, pxl_row_in_fullres, pxl_col_in_fullres)
colnames(coord) <- c("barcodes", "xcoord", "ycoord")
coord <- coord %>% column_to_rownames("barcodes")

# Extract counts and nUMI
counts_sp <- vNECC1$Spatial@counts %>% as.data.frame()
nUMI <- colSums(counts_sp)

# Create SpatialRNA object
puck <- SpatialRNA(coord, counts_sp, nUMI)


# Create RCTD object
myRCTD <- create.RCTD(puck, reference, max_cores = 1, CELL_MIN_INSTANCE = 5)

# Run RCTD in full mode
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

# Save results
saveRDS(myRCTD, "vNECC1_RCTD.rds")

# Extract weights and normalize
results <- myRCTD@results
norm_weights <- sweep(results$weights, 1, rowSums(results$weights), "/")
norm_weights <- as.data.frame(as.matrix(norm_weights))

# Add results to Seurat object metadata
vNECC1@meta.data <- vNECC1@meta.data[, c(1:3)]
meta <- vNECC1@meta.data %>% rownames_to_column("spot")
df <- norm_weights %>% rownames_to_column("spot")
vNECC1@meta.data <- meta %>% left_join(df, by = "spot") %>% column_to_rownames("spot")

# Calculate tumor proportion
vNECC1@meta.data <- vNECC1@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(Tumor = CA_NE + CA_Columnar + CA_Glandular + CA_Squamous + CA_Ciliated) %>%
  column_to_rownames("cell")

# Save final object
saveRDS(vNECC1, "vNECC1_after_decov.rds")

# Plot spatial scatterpie
library(RColorBrewer)

# Define color palette
my_col <- c("#E8602D", "#F3E55C", "#359566", "#a6dbbb", "#967bce", "#ecd9f1", "#075A84")

prop <- vNECC1@meta.data %>%
  select(T, NK, B, Plasma, Monocyte, Macrophage, DC, Neutrophil, Mast,
         Fibroblast, SMC, Pericyte, Endo, Tumor) %>%
  as.matrix()

plotSpatialScatterpie(
  x = vNECC1,
  y = prop,
  cell_types = colnames(prop),
  img = FALSE,
  scatterpie_alpha = 1,
  pie_scale = 0.6,
  degrees = -90,
  axis = "h"
) + scale_fill_manual(values = my_col)


# 2. Spati: Spatial analysis of cell type distribution in NECC1 sample


# Load required libraries
library(Giotto)
library(spati)
library(tidyverse)

# Load pre-processed data
vNECC1 <- readRDS("vNECC1_after_decov.rds")

# Convert normalized weights to integer counts
num.df <- round(10 * norm_weights)
n.sum <- rowSums(num.df) %>% data.frame()
colnames(n.sum) <- 'n'

num.df <- num.df %>% rownames_to_column("spot")
n.sum <- n.sum %>% rownames_to_column("spot")
num.df <- num.df %>% left_join(n.sum, by = "spot")

# Expand to individual cell-level data
cell.combine <- data.frame(cell = "")
for(i in rownames(num.df)) {
    n <- num.df[i,]$n
    cell.barcode <- paste0(num.df[i,]$spot, "_", 1:n)
    cell.df <- data.frame(cell = cell.barcode)
    cell.combine <- rbind(cell.df, cell.combine)
}

cell.combine$spot <- substr(cell.combine$cell, 1, 18)
cell.all <- cell.combine %>% left_join(num.df, by = "spot") %>% na.omit()

# Assign cell types to each cell
spot.cell.combine <- NULL
for(i in unique(cell.all$spot)) {
    spot.one <- cell.all %>% filter(spot == i)
    
    celltype.v <- c(
        rep("T", unique(spot.one$T)),
        rep("B", unique(spot.one$B)),
        rep("Plasma", unique(spot.one$Plasma)),
        rep("NK", unique(spot.one$NK)),
        rep("Mast", unique(spot.one$Mast)),
        rep("Tumor", unique(spot.one$Tumor)),
        rep("DC", unique(spot.one$DC)),
        rep("Macrophage", unique(spot.one$Macrophage)),
        rep("Monocyte", unique(spot.one$Monocyte)),
        rep("Neutrophil", unique(spot.one$Neutrophil)),
        rep("SMC", unique(spot.one$SMC)),
        rep("Fibroblast", unique(spot.one$Fibroblast)),
        rep("Pericyte", unique(spot.one$Pericyte)),
        rep("Endo", unique(spot.one$Endo))
    )
    
    cell.df <- data.frame(celltype = celltype.v)
    spot.cell <- cbind(spot.one, cell.df)
    spot.cell.combine <- rbind(spot.cell, spot.cell.combine)
}

# Load spatial coordinates
coord <- read.csv("vNECC1CA/outs/spatial/tissue_positions.csv")
df <- spot.cell.combine %>% select(cell, spot, n)
coord <- df %>% left_join(coord, by = c("spot" = "barcode"))

# Generate spatial points within each spot
coord.pxl <- NULL
for(i in unique(coord$spot)) {
    coord.one <- coord %>% filter(spot == i)
    pxl_row <- unique(coord.one$pxl_row_in_fullres)
    pxl_col <- unique(coord.one$pxl_col_in_fullres)
    range.pxl_row <- seq(from = pxl_row, to = pxl_row + 1, by = 0.05)
    range.pxl_col <- seq(from = pxl_col, to = pxl_col + 1, by = 0.05)
    n <- unique(coord.one$n)
    
    sample.row <- sample(range.pxl_row, size = n)
    sample.col <- sample(range.pxl_col, size = n)
    sample_pxl.df <- data.frame(sample_pxl_row = sample.row, sample_pxl_col = sample.col)
    coord.sample.pxl <- cbind(coord.one, sample_pxl.df)
    coord.pxl <- rbind(coord.sample.pxl, coord.pxl)
}

# Extract locations for spati analysis
locs <- coord.pxl %>% select(sample_pxl_row, sample_pxl_col)

# Create Delaunay triangulation for spatial analysis
S <- spati(locs, center = TRUE)

# Define cell type colors and plot
L <- list(
    B = list(pt = grepl("B", spot.cell.combine$celltype), col = '#917099', rdot = 10),
    Plasma = list(pt = grepl("Plasma", spot.cell.combine$celltype), col = '#8861AC', rdot = 10),
    T = list(pt = grepl("T", spot.cell.combine$celltype), col = '#DEB969', rdot = 10),
    NK = list(pt = grepl("NK", spot.cell.combine$celltype), col = "#E7E099", rdot = 10),
    Mast = list(pt = grepl("Mast", spot.cell.combine$celltype), col = "#F79C5D", rdot = 10),
    Tumor = list(pt = grepl("Tumor", spot.cell.combine$celltype), col = "yellow", rdot = 10),
    Macrophage = list(pt = grepl("Macrophage", spot.cell.combine$celltype), col = "#E39970", rdot = 10),
    Monocyte = list(pt = grepl("Monocyte", spot.cell.combine$celltype), col = "#BFA5CF", rdot = 10),
    Neutrophil = list(pt = grepl("Neutrophil", spot.cell.combine$celltype), col = "#FDA746", rdot = 10),
    DC = list(pt = grepl("DC", spot.cell.combine$celltype), col = "#FE8205", rdot = 10),
    Fibroblast = list(pt = grepl("Fibroblast", spot.cell.combine$celltype), col = "#EB494A", rdot = 10),
    SMC = list(pt = grepl("SMC", spot.cell.combine$celltype), col = "#F99392", rdot = 10),
    Pericyte = list(pt = grepl("Pericyte", spot.cell.combine$celltype), col = "#40A635", rdot = 10),
    Endo = list(pt = grepl("Endo", spot.cell.combine$celltype), col = "#E83C2D", rdot = 10)
)

# Generate plots
spati_plot(S,L=L)

# make the color more faint
spati_plot(S,L=L,alpha=0.25)

# zoom in
spati_plot(S,L=L,radius=200)

# move the center
spati_plot(S,L=L,radius=500,center=c(1000,200))

# show with a zoom panel
spati_plot2(S,L=L,zradius=500,zcenter=c(1000,200))

spati_plot2(S,L=L,zradius=355,zcenter=c(-720,585))
spati_edges(S,col="black")
spati_plot(S,L=L[c(names(L))])
#spati_plot(S,L=L[c("Tumor","Stroma")])
# define a new adjacency matrix of edges
A355 <- filter_long_edges( S$A, 355) # elimate edges with greater than 50 µm
spati_edges(S,E=A355$E)


Y <- nh_sum(
A355,   # adjacency matrix
S$P,   # cell coordinates
V=cbind(                    # multiple neighborhood summaries
all=rep(1,nrow(S$P)),
Tumor=L$Tumor$pt,
Macrophage=L$Macrophage$pt,
Monocyte=L$Monocyte$pt,
Pericyte=L$Pericyte$pt,
Plasma=L$Plasma$pt,
Neutrophil=L$Neutrophil$pt,
DC=L$DC$pt,
T=L$T$pt,
Endo=L$Endo$pt,
Fibroblast=L$Fibroblast$pt,
NK=L$NK$pt,
SMC=L$SMC$pt,
B=L$B$pt,
Mast=L$Mast$pt
),
W=1:10   # neighborhood of 1,2,..,10 layers of cells (for example)
)

str(Y)

tumor_area <- ae_open( A355,ae_close( A355, L$Tumor$pt))

spati_plot2(S,L=L)
spati_plot2(S,L=c(
    list(tumor_area=list(pt=tumor_area,dot=15,col="green")),
    L))


tl <- sapply(-5:5,function(i) list(pt=(tumor_layers==i),
                                    col=hsv(0.45+i/15)),simplify=F)
names(tl) <- paste0("layer ",-5:5)
spati_plot2(S,L=tl)

ltab <- table( tumor_layers, spot.cell.combine$celltype == "Tumor")
Tumor.layer.data <- data.frame(ratio = ltab[,2]/rowSums(ltab),celltype = "Tumor") %>%
    rownames_to_column("layer")

ltab <- table( tumor_layers, spot.cell.combine$celltype == "Macrophage")
Macrophage.layer.data <- data.frame(ratio = ltab[,2]/rowSums(ltab),celltype = "Macrophage") %>%
    rownames_to_column("layer")
 
ltab <- table( tumor_layers, spot.cell.combine$celltype == "SMC")
SMC.layer.data <- data.frame(ratio = ltab[,2]/rowSums(ltab),celltype = "SMC") %>%
    rownames_to_column("layer")

ltab <- table( tumor_layers, spot.cell.combine$celltype == "Fibroblast")
Fibroblast.layer.data <- data.frame(ratio = ltab[,2]/rowSums(ltab),celltype = "Fibroblast") %>%
    rownames_to_column("layer")

ltab <- table( tumor_layers, spot.cell.combine$celltype == "Endo")
Endo.layer.data <- data.frame(ratio = ltab[,2]/rowSums(ltab),celltype = "Endo") %>%
    rownames_to_column("layer")

ltab <- table( tumor_layers, spot.cell.combine$celltype == "Pericyte")
Pericyte.layer.data <- data.frame(ratio = ltab[,2]/rowSums(ltab),celltype = "Pericyte") %>%
    rownames_to_column("layer")

ltab <- table( tumor_layers, spot.cell.combine$celltype == "T")
T.layer.data <- data.frame(ratio = ltab[,2]/rowSums(ltab),celltype = "T") %>%
    rownames_to_column("layer")

ltab <- table( tumor_layers, spot.cell.combine$celltype == "DC")
DC.layer.data <- data.frame(ratio = ltab[,2]/rowSums(ltab),celltype = "DC") %>%
    rownames_to_column("layer")

ltab <- table( tumor_layers, spot.cell.combine$celltype == "Mast")
Mast.layer.data <- data.frame(ratio = ltab[,2]/rowSums(ltab),celltype = "Mast") %>%
    rownames_to_column("layer")

ltab <- table( tumor_layers, spot.cell.combine$celltype == "B")
B.layer.data <- data.frame(ratio = ltab[,2]/rowSums(ltab),celltype = "B") %>%
    rownames_to_column("layer")


ltab <- table( tumor_layers, spot.cell.combine$celltype == "Monocyte")
Monocyte.layer.data <- data.frame(ratio = ltab[,2]/rowSums(ltab),celltype = "Monocyte") %>%
    rownames_to_column("layer")

ltab <- table( tumor_layers, spot.cell.combine$celltype == "NK")
NK.layer.data <- data.frame(ratio = ltab[,2]/rowSums(ltab),celltype = "NK") %>%
    rownames_to_column("layer")


ltab <- table( tumor_layers, spot.cell.combine$celltype == "Plasma")
Plasma.layer.data <- data.frame(ratio = ltab[,2]/rowSums(ltab),celltype = "Plasma") %>%
    rownames_to_column("layer")


layer.data <- Tumor.layer.data %>% rbind(
                                      Neutrophil.layer.data,
                                      T.layer.data,
                                      NK.layer.data,
                                      B.layer.data,
                                      Plasma.layer.data,
                                      Mast.layer.data,
                                      SMC.layer.data,
                                      Endo.layer.data,
                                      Fibroblast.layer.data,
                                      Pericyte.layer.data,
                                      DC.layer.data,
                                      Monocyte.layer.data,
                                      Macrophage.layer.data
                                            )
#layer.data <- stroma.layer.data %>% rbind(Tumor.layer.data)
layer.data$layer <- as.numeric(layer.data$layer)
saveRDS(layer.data,"NECC1CA_layer.data.rds")



# 3.Neighborhood Analysis

# Load required packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)


# Read Seurat object
vNECC1 <- readRDS("vNECC1_after_decov.rds")

# Read spatial coordinates
library(readr)
los <- read.csv("/vNECC1/outs/spatial/tissue_positions.csv")
los <- los %>%
  dplyr::select(barcode, array_col, array_row) %>%
  dplyr::rename(spot = barcode, xcoord = array_col, ycoord = array_row)

los$xcoord <- as.numeric(los$xcoord)
los$ycoord <- as.numeric(los$ycoord)

# Scale coordinates
los <- los %>%
  column_to_rownames("spot") %>%
  scale(scale = c(1, -1)) %>%
  as.data.frame() %>%
  rownames_to_column("spot") %>%
  filter(spot %in% colnames(vNECC1))

# Get cell type names (columns 4 to 67 in meta.data)
celltype.name <- colnames(vNECC1@meta.data)[4:67]

# Get original coordinates
coord <- vNECC1@images$vNECC1@coordinates %>% rownames_to_column("spot")


# Calculate Neighborhood Cell Type Composition for Each Spot
neighbor.combined <- NULL
meta <- vNECC1@meta.data %>% rownames_to_column("spot")

for (i in coord$spot) {
  temp <- coord %>% filter(spot == i)
  row0 <- temp$row
  col0 <- temp$col
  
  # Define neighborhood range (row ±1, col ±2)
  row.neighbor <- seq(row0 - 1, row0 + 1)
  col.neighbor <- seq(col0 - 2, col0 + 2)
  
  temp2 <- coord %>% filter(row %in% row.neighbor & col %in% col.neighbor)
  spot.neighbor <- temp2$spot
  
  meta.neighbor <- meta %>% filter(spot %in% spot.neighbor)
  
  # Calculate neighbor sums
  neighbor_sums <- colSums(meta.neighbor[, celltype.name, drop = FALSE], na.rm = TRUE)
  
  spot.self <- meta %>% filter(spot == i)
  
  # Neighbor sum minus self = true neighbors
  for (ct in celltype.name) {
    spot.self[[paste0(ct, ".Neighbor")]] <- neighbor_sums[ct] - spot.self[[ct]]
  }
  
  neighbor.combined <- rbind(spot.self, neighbor.combined)
}

# Convert to matrix for correlation calculation
mtx <- neighbor.combined %>% column_to_rownames("spot")


#Calculate Correlations Between Cell Types
neighbor.cor <- cor(mtx)


# 4.Cottrazm
library(Cottrazm)
library(Seurat)
library(ggplot2)
library(ggpubr)

print('STPreProcess')
InDir = "vsNECC1CA/outs/"
Sample = "NECC1"
OutDir = paste("~/1_Project/NECC/data_final/NEW/sp/Cottrazm","/",Sample,"/",sep = "")

TumorST <-
  STPreProcess(
    InDir = InDir,
    OutDir = OutDir,
    Sample = Sample
    )


STModiCluster_PyFree <- function(TumorST,
                                 res = 1.5,
                                 Sample = "Sample",
                                 OutDir = NULL) {
  if (is.null(OutDir)) {
    OutDir <- file.path(getwd(), Sample)
    dir.create(OutDir, recursive = TRUE, showWarnings = FALSE)
  }

  # Use Spatial counts directly as "Morph" assay
  MorphMatirxSeurat <- CreateSeuratObject(counts = TumorST@assays$Spatial@counts)
  TumorST@assays$Morph <- MorphMatirxSeurat@assays$RNA

  # Normalize, find variable features, scale, PCA, neighbors, UMAP, clustering
  TumorST <- NormalizeData(TumorST, assay = "Morph")
  TumorST <- FindVariableFeatures(TumorST, assay = "Morph", mean.function = ExpMean, dispersion.function = LogVMR)
  TumorST <- ScaleData(TumorST, assay = "Morph")
  TumorST <- RunPCA(TumorST, assay = "Morph", npcs = 50, verbose = FALSE)
  TumorST <- FindNeighbors(TumorST, reduction = "pca", dims = 1:50, assay = "Morph")
  TumorST <- RunUMAP(TumorST, reduction = "pca", dims = 1:50, assay = "Morph")
  TumorST <- FindClusters(TumorST, resolution = res, algorithm = 1, graph.name = "Morph_snn")
  TumorST@meta.data$seurat_clusters <- TumorST@meta.data[, paste0("Morph_snn_res.", res)]

  # Define cluster colors
  cluster_colors <- c(
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
    "#33A02C", "#B2DF8A", "#55B1B1", "#8DD3C7", "#A6761D",
    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"
  )

  # SpatialDimPlot
  pdf(file.path(OutDir, paste0(Sample, "_Spatial_SeuratCluster.pdf")), width = 7, height = 7)
  p <- SpatialDimPlot(TumorST, group.by = "seurat_clusters", cols = cluster_colors, pt.size.factor = 1, alpha = 0.8) +
    scale_fill_manual(values = cluster_colors) +
    labs(title = paste("Resolution =", res))
  print(p)
  dev.off()

  # UMAP plot
  pdf(file.path(OutDir, paste0(Sample, "_UMAP_SeuratCluster.pdf")), width = 7, height = 7)
  p <- DimPlot(TumorST, group.by = "seurat_clusters", cols = cluster_colors) +
    scale_fill_manual(values = cluster_colors) +
    labs(title = paste("Resolution =", res))
  print(p)
  dev.off()

  # Add NormalScore
  NormalFeatures <- c("PTPRC","CD2","CD3D","CD3E","CD3G","CD5","CD7","CD79A","MS4A1","CD19")
  present_features <- NormalFeatures[NormalFeatures %in% rownames(TumorST@assays$Morph@data)]
  TumorST@meta.data$NormalScore <- colMeans(TumorST@assays$Morph@data[present_features, , drop = FALSE])

  # NormalScore violin plot
  pdf(file.path(OutDir, paste0(Sample, "_NormalScore.pdf")), width = 6, height = 4)
  p <- VlnPlot(TumorST, features = "NormalScore", group.by = "seurat_clusters", pt.size = 0, cols = cluster_colors) +
    geom_boxplot(width = 0.1) +
    ggpubr::stat_compare_means() +
    ggplot2::theme(legend.position = "none")
  print(p)
  dev.off()

  # Identify NormalCluster (highest NormalScore)
  cluster_means <- tapply(TumorST@meta.data$NormalScore, TumorST@meta.data$seurat_clusters, mean)
  NormalCluster <- names(which.max(cluster_means))
  message("NormalCluster = ", NormalCluster)

  # Save cell annotation for InferCNV
  cellAnnotation <- data.frame(CellID = rownames(TumorST@meta.data),
                               DefineTypes = TumorST@meta.data$seurat_clusters)
  dir.create(file.path(OutDir, "InferCNV"), showWarnings = FALSE)
  write.table(cellAnnotation, file.path(OutDir, "InferCNV/CellAnnotation.txt"),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  return(TumorST)
}


print('STModiCluster')
res = 1.5
TumorST <- STModiCluster_PyFree(TumorST,
                                 res = res,
                                 Sample = Sample,
                                 OutDir = OutDir )


print('STCNV' )
STInferCNV <-
  STCNV(
    TumorST = TumorST ,
    OutDir  = OutDir,
    assay = "Spatial"
    )


TumorST <-
  STCNVScore2(
    TumorST = TumorST,
    assay = "Spatial",
    Sample = Sample,
    OutDir = OutDir
    )

### 3.1.4 Define of tumor boundary
TumorSTn <-
  BoundaryDefine(
    TumorST = TumorST,
    MalLabel = c('all_observations.all_observations.1.1.1.1','all_observations.all_observations.1.1.1.2','all_observations.all_observations.1.1.2.1','all_observations.all_observations.1.1.2.2',
                'all_observations.all_observations.1.2.1.1','all_observations.all_observations.1.2.1.2'),
    OutDir = OutDir,
    Sample = Sample
  )


TumorST <-
  BoundaryPlot(
    TumorSTn = TumorSTn,
    TumorST = TumorST,
    OutDir = OutDir,
    Sample = Sample
  )


options(repr.plot.width =6, repr.plot.height =6)
  TumorST_DefineColors <- data.frame(
    Location = c("Mal", "Bdy", "nMal"),
    colors = c("#CB181D", "#1f78b4", "#fdb462")
  )
  TumorST@meta.data$Location <- factor(TumorST@meta.data$Location, levels = TumorST_DefineColors$Location)
  TumorST_DefineColors$colors <- as.character(TumorST_DefineColors$colors)

  #pdf(paste(OutDir, Sample, "_BoundaryDefine.pdf", sep = ""), width = 7, height = 7)
  p <- SpatialDimPlot(TumorST, group.by = "Location", cols = TumorST_DefineColors$colors,pt.size.factor = 1.7) +
    scale_fill_manual(values = TumorST_DefineColors$colors )
  print(p)
dev.off()

SpatialFeaturePlot(TumorST,features = 'cnv_score',pt.size.factor = 1.8)

saveRDS(TumorST,"/Cottrazm/NECC1/TumorST_NECC1.rds")


# 5. Integration of Spatial Transcriptomics Data
# Load Required Libraries
library(Seurat)
library(Matrix)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyverse)


# Load Data Objects (RDS files)
# Update paths to relative paths for GitHub portability
vNECC1 <- readRDS("data/TumorST_NECC1_HM.rds")
vNECC2 <- readRDS("data/TumorST_NECC2_HM.rds")
vNECC3 <- readRDS("data/TumorST_NECC3_HM.rds")
vNECC4 <- readRDS("data/TumorST_NECC4_HM.rds")

# Preprocessing: Add Slice Metadata
vNECC1[["slice"]] <- "vNECC1"
vNECC2[["slice"]] <- "vNECC2"
vNECC3[["slice"]] <- "vNECC3"
vNECC4[["slice"]] <- "vNECC4"

# Filter Empty Genes (Remove genes with zero expression across all spots)
filter_empty_genes <- function(seu_obj) {
  counts <- GetAssayData(seu_obj, assay = "Spatial", slot = "counts")
  keep_genes <- rowSums(as.matrix(counts)) != 0
  return(seu_obj[keep_genes, ])
}

vNECC1 <- filter_empty_genes(vNECC1)
vNECC2 <- filter_empty_genes(vNECC2)
vNECC3 <- filter_empty_genes(vNECC3)
vNECC4 <- filter_empty_genes(vNECC4)

# Assign Sample IDs
vNECC1[["sample"]] <- "NECC1"
vNECC2[["sample"]] <- "NECC2"
vNECC3[["sample"]] <- "NECC3"
vNECC4[["sample"]] <- "NECC4"

vNECC <- merge(vNECC1,c(vNECC2,vNECC3,vNECC4))

DefaultAssay(vNECC) <- 'Spatial'
NECC.split <- SplitObject(vNECC,split.by = "sample")
NECC.split <- lapply(X = NECC.split, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize",
                     scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2500)
})
features <- SelectIntegrationFeatures(object.list = NECC.split,
                                      nfeatures = 2500)
NECC.split <- lapply(X = NECC.split, FUN = function(x){
    x <- ScaleData(x, vars.to.regress = c("nCount_Spatial","nFeature_Spatial",
                                         "percent.mt"),
                 features = features.sub)
     x <- RunPCA(x, npcs = 30,features = features)
})

anchors <- FindIntegrationAnchors(NECC.split,
                                  k.anchor = 5,
                                  dims = 1:30,
                                  k.filter =200,
                                  k.score = 30,
                                  anchor.features = features,
                                  reduction = 'cca')

NECC_integrated <- IntegrateData(anchorset = anchors,
                                dims = 1:30,
                                k.weight = 100)

DefaultAssay(NECC_integrated) <- 'integrated'
var.features  <- VariableFeatures(NECC_integrated)
NECC_integrated <- ScaleData(NECC_integrated)
                           #   vars.to.regress = c("G2M.Score","S.Score"))
NECC_integrated <- RunPCA(NECC_integrated,
                          features = var.features,
                         dims = 1:30)
NECC_integrated <- FindNeighbors(NECC_integrated, reduction = "pca", dims = 1:30)
NECC_integrated <- FindClusters(NECC_integrated, resolution = c(0.1, 0.2, 0.3,0.4, 0.5))
NECC_integrated <- RunUMAP(NECC_integrated, reduction = "pca", dims = 1:30)

options(repr.plot.width =6, repr.plot.height =5)
DimPlot(object = NECC_integrated , reduction = 'umap', group.by = "sample")
DimPlot(object = NECC_integrated , reduction = 'umap', group.by = "Location")
DimPlot(object = NECC_integrated , reduction = 'umap',label = T, group.by = "Phase")
FeaturePlot(NECC_integrated , features = 'nFeature_Spatial')
FeaturePlot(NECC_integrated , features = 'percent.mt')
FeaturePlot(NECC_integrated , features = 'S.Score')
FeaturePlot(NECC_integrated , features = 'G2M.Score')
