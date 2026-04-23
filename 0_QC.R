# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(DoubletFinder)

# Read sample metadata
sample.df <- read.csv("samplesheet_gse.csv")

# Function to run Scrublet for doublet detection
run_scrublet <- function(seu, npcs = 20, doublet_rate = 0.05){
  cells <- rownames(seu@meta.data)
  mat <- t(as.matrix(seu[["RNA"]]@counts[,cells]))
  scrub <- reticulate::import("scrublet")$Scrublet(mat, expected_doublet_rate = doublet_rate)
  res <- scrub$scrub_doublets(n_prin_comps = as.integer(npcs))
  seu@meta.data$scrublet_score <- res[[1]]
  seu@meta.data$scrublet_call <- res[[2]]
  return(seu)
}

# Load and process each sample
sample.combined <- NULL
for(batch_id in unique(sample.df$Batch)){
  cat("Loading data:", batch_id, "\n")
  
  # Read 10X data
  sample.count <- Read10X_h5(paste0("sc", batch_id, "gex/raw_feature_bc_matrix_cb_filtered.h5"))
  colnames(sample.count) <- paste0(batch_id, "_", colnames(sample.count))
  
  # Create Seurat object
  sample <- CreateSeuratObject(sample.count, project = batch_id, min.cells = 3, min.features = 10)
  sample <- NormalizeData(sample)
  
  # Run doublet detection
  dbl_rate <- sample.df$DblRate[sample.df$Batch == batch_id][1]
  sample <- run_scrublet(sample, npcs = 20, doublet_rate = dbl_rate)
  
  # Run DoubletFinder
  sample <- NormalizeData(sample)
  sample <- FindVariableFeatures(sample, nfeatures = 2000)
  sample <- ScaleData(sample, vars.to.regress = "nCount_RNA")
  sample <- RunPCA(sample)
  
  nExp <- round(dbl_rate * ncol(sample))
  sample <- doubletFinder_v3(sample, PCs = 1:20, nExp = nExp, pN = 0.25, pK = 0.09)
  
  # Merge samples
  if(is.null(sample.combined)){
    sample.combined <- sample
  } else {
    sample.combined <- merge(sample.combined, sample)
  }
}

# Quality control filtering
cell.use <- sample.combined@meta.data %>%
  rownames_to_column("Cell") %>%
  filter(nFeature_RNA > 250, nCount_RNA > 500, percent.mt < 20) %>%
  pull(Cell)

merge <- subset(sample.combined, cell = cell.use)

# Add metadata from sample sheet
df <- merge@meta.data %>%
  rownames_to_column("cell") %>%
  left_join(sample.df, by = c("orig.ident" = "Sample")) %>%
  column_to_rownames("cell")
merge@meta.data <- df

# Add module scores (IFN, Stress, Hypoxia)
interferon.genes <- c("STAT1", "IRF1", "MX1", "OAS1", "ISG15", "IFIT1", "IFIT2", "IFIT3")
stress.genes <- c("JUN", "FOS", "ATF3", "EGR1", "HSPA1A", "HSPA1B", "HSP90AA1")
hypoxia.genes <- c("VEGFA", "SLC2A1", "LDHA", "PDK1", "HK2", "PGK1")

merge <- AddModuleScore(merge, features = list(interferon.genes), name = "IFN.Score")
merge <- AddModuleScore(merge, features = list(stress.genes), name = "Stress.Score")
merge <- AddModuleScore(merge, features = list(hypoxia.genes), name = "Hypoxia.Score")

# Save final object
saveRDS(merge, "integrated_final.rds")
