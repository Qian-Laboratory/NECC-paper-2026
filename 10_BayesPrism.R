# Load required libraries
library(BayesPrism)
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)

# Read integrated Seurat object
merge.clean <- readRDS("~/1_Project/NECC/data_final/all/merge_integrated_celltype.rds")

# Subset to cancer samples (ADC, SCC, NECC)
merge.CA <- subset(merge.clean, class %in% c("ADC", "SCC", "NECC"))

# Downsample to 20,000 cells for computational efficiency
cell <- colnames(merge.CA)
cell.sample <- sample(cell, size = 20000)
merge.CA.sample <- subset(merge.CA, cell = cell.sample)

# Set cell state and cell type labels
merge.CA.sample$cell.state.labels <- merge.CA.sample$celltype.class
merge.CA.sample$cell.type.labels <- merge.CA.sample$celltype.class

# Extract single-cell count matrix (cells as rows, genes as columns)
sc.dat <- t(merge.CA.sample$RNA@counts) %>% as.matrix()

# Load TCGA bulk RNA-seq data
expdata <- readRDS("/TCGA_CESC/TCGA_CESC_expdata.rds")

# Load gene annotation
anno.df <- read.csv("/TCGA_CESC/anno_counts.csv", header = FALSE) %>% select(-V3)
colnames(anno.df) <- c("geneID", "gene")

# Extract TPM matrix
library(SummarizedExperiment)
tpm_matrix <- assays(expdata)$"tpm_unstrand"

# Process TPM data: annotate and remove duplicates
tpm.df <- as.data.frame(tpm_matrix) %>% rownames_to_column("geneID")
tpm.anno.df <- anno.df %>% left_join(tpm.df, by = "geneID") %>% select(-geneID)
tpm.unique.df <- tpm.anno.df[!duplicated(tpm.anno.df$gene), ]
tpm.unique.df <- tpm.unique.df %>% rownames_to_column("n") %>% select(-n)
tpm.mtx <- tpm.unique.df %>% column_to_rownames("gene")

# Bulk data matrix (samples as rows, genes as columns)
bk.dat <- t(tpm.mtx)

# Extract cell type and state labels
cell.type.labels <- merge.CA.sample@meta.data$cell.type.labels
cell.state.labels <- merge.CA.sample@meta.data$cell.state.labels

# Plot correlation between cell states
cell.state.corr.plot <- plot.cor.phi(input = sc.dat,
                                      input.labels = cell.state.labels,
                                      title = "cell state correlation",
                                      cexRow = 0.4, cexCol = 0.4,
                                      margins = c(2, 2))

# Plot correlation between cell types
cell.type.corr.plot <- plot.cor.phi(input = sc.dat,
                                    input.labels = cell.type.labels,
                                    title = "cell type correlation",
                                    cexRow = 0.5, cexCol = 0.5)

# Set up parallel processing
plan("multicore", workers = 5)
plan()
options(future.globals.maxSize = 20 * 1000 * 1024^2)

# Detect outliers in single-cell data
sc.stat <- plot.scRNA.outlier(
  input = sc.dat,
  cell.type.labels = cell.type.labels,
  species = "hs",
  return.raw = TRUE
)

# Detect outliers in bulk data relative to single-cell reference
bk.stat <- plot.bulk.outlier(
  bulk.input = bk.dat,
  sc.input = sc.dat,
  cell.type.labels = cell.type.labels,
  species = "hs",
  return.raw = TRUE
)

# Clean up genes: remove ribosomal, mitochondrial, and other unwanted genes
sc.dat.filtered <- cleanup.genes(input = sc.dat,
                                  input.type = "count.matrix",
                                  species = "hs",
                                  gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"),
                                  exp.cells = 5)

# Plot single-cell vs bulk gene expression comparison
plot.bulk.vs.sc(sc.input = sc.dat.filtered,
                bulk.input = bk.dat)

# Increase memory limit for large-scale computation
options(future.globals.maxSize = 35 * 1000 * 1024^2)

# Perform differential expression analysis between cell types/states
diff.exp.stat <- get.exp.stat(sc.dat = sc.dat[, colSums(sc.dat > 0) > 3],
                              cell.type.labels = cell.type.labels,
                              cell.state.labels = cell.state.labels,
                              pseudo.count = 0.1,
                              cell.count.cutoff = 50,
                              n.cores = 3)

# Select marker genes based on differential expression statistics
sc.dat.filtered.pc.sig <- select.marker(sc.dat = sc.dat.filtered.pc,
                                        stat = diff.exp.stat,
                                        pval.max = 0.01,
                                        lfc.min = 0.1)

# Initialize BayesPrism object for deconvolution
myPrism <- new.prism(
  reference = sc.dat.filtered.pc,
  mixture = bk.dat,
  input.type = "count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key = NULL,
  outlier.cut = 0.01,
  outlier.fraction = 0.1,
)

# Run BayesPrism deconvolution
bp.res <- run.prism(prism = myPrism, n.cores = 5)

# Save results
saveRDS(bp.res, "/TCGA/bp.res.rds")

# Extract posterior mean of cell type fractions (theta)
theta <- get.fraction(bp = bp.res,
                      which.theta = "final",
                      state.or.type = "type")

# Extract coefficient of variation (CV) of cell type fractions
theta.cv <- bp.res@posterior.theta_f@theta.cv
