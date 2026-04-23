# Load required libraries
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(slingshot)

# Read integrated T cell data
T_integrated <- readRDS("T_integrated_celltype.rds")

# Subset CD8+ T cell subtypes
CD8T <- subset(T_integrated, celltype.sub %in% c("CD8_TEM_GZMK", "CD8_TEM_IL7R",
                                                  "CD8_TEMRA_CX3CR1", "CD8_TEX_CXCL13",
                                                  "CD8_TN", "CD8_TRM_GNLY", "CD8_TRM_IFNG",
                                                  "CD8_TRM_IL7R"))

# Extract UMAP embeddings and cell identities
rd <- Embeddings(CD8T, "umap")
Idents(CD8T) <- "celltype.new.sub"
cl <- Idents(CD8T)

# Define lineage starting from naive CD8 T cells
lineage <- getLineages(rd, cl, start.clus = c("CD8_TN"))

# Store slingshot object
CD8T@misc[['slingshot']] <- list("dr" = "pca", "data" = ss)

# Get pseudotime values
pst <- slingPseudotime(CD8T@misc[['slingshot']]$data)
CD8T@meta.data[, colnames(pst)] <- as.data.frame(pst)

# Convert pseudotime to percentage scale for each lineage
for (cv in colnames(pst)) {
  maxtime <- max(CD8T@meta.data[, cv], na.rm = TRUE)
  cv_pct <- paste0(cv, "_pct")
  CD8T@meta.data[, cv_pct] <- (CD8T@meta.data[, cv] / maxtime) * 100
}

# Calculate mean percentage across lineages
metadata.curve <- CD8T@meta.data %>%
  tibble::rownames_to_column(var = "Cell") %>%
  rowwise() %>%
  dplyr::mutate(curve_mean = mean(c(Lineage1_pct, Lineage2_pct, Lineage3_pct), na.rm = TRUE)) %>%
  tibble::column_to_rownames(var = "Cell")
CD8T@meta.data <- metadata.curve

# Function to plot pseudotime trajectory
plotPseudoTime <- function(object,
                          group.by = NULL,
                          reduction = 'umap',
                          dims = 1:2,
                          pt.size = pt.size,
                          no.curve = FALSE,
                          colors = NULL) {
  object[['ident']] <- Idents(object = object)
  group.by <- group.by %||% 'ident'
  dims <- paste0(Seurat::Key(object = object[[reduction]]), dims)

  # Extract curve coordinates
  curved <- bind_rows(lapply(names(object@misc$slingshot$data@metadata$curves), function(x) {
    c <- slingCurves(object@misc$slingshot$data)[[x]]
    d <- as.data.frame(c$s[c$ord, dims])
    d$curve <- x
    return(d)
  }))

  data <- FetchData(object = object, vars = c(dims, group.by))
  p <- ggplot(data, aes_string(x = dims[1], y = dims[2])) +
    geom_point(aes(color = !!sym(group.by)), size = pt.size) +
    theme(legend.position = "top")
  
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }

  if (is.character(data[[group.by]]) | is.factor(data[[group.by]])) {
    p <- p + guides(col = guide_legend(nrow = 2))
  } else {
    p <- p + scale_color_distiller(palette = "YlOrRd", na.value = 'grey90')
  }
  
  if (!no.curve) {
    p <- p + geom_path(aes_string(dims[1], dims[2], linetype = "curve"), curved, size = 1)
  }
  return(p)
}

# Plot pseudotime trajectory
options(repr.plot.width = 10.5, repr.plot.height = 4)
plotPseudoTime(CD8T, group.by = "celltype_sub",
               reduction = "umap",
               dims = 1:2,
               pt.size = 0.5,
               colors = c(rep('grey80', 8))) +
  theme_classic()

# Define gene signatures
activated.gene <- c("FAS", "FASLG", "CD44", "CD69", "CD38", "NKG7", "KLRB1", "KLRD1", "KLRF1", "KLRG1",
                    "KLRK1", "FCGR3A", "CX3CR1", "FGFBP2")

Cytotoxicity.gene <- c("GZMA", "GZMB", "GZMH", "GZMK", "GZMH", "GNLY", "PRF1", "IFNG", "TNF",
                       "SERPINB1", "SERPINB6", "SERPINB9", "CTSA", "CTSB", "CTSC", "CTSD",
                       "CTSW", "CST3", "CST7", "CSTB", "LAMP1", "LAMP3", "CAPN2")

exh.gene <- c("PDCD1", "LAYN", "HAVCR2", "LAG3", "CD244", "CTLA4",
              "TIGIT", "TOX", "BTLA", "ENTPD1")

# Add module scores for exhaustion
CD8T <- AddModuleScore(
  object = CD8T,
  features = list(exh.gene),
  ctrl = length(exh.gene),
  name = 'Exh.Score',
  assay = "RNA")

# Add module scores for activation
CD8T <- AddModuleScore(
  object = CD8T,
  features = list(activated.gene),
  ctrl = length(activated.gene),
  name = 'Activated.Score',
  assay = "RNA")

# Add module scores for cytotoxicity
CD8T <- AddModuleScore(
  object = CD8T,
  features = list(Cytotoxicity.gene),
  ctrl = length(Cytotoxicity.gene),
  name = 'Cytotoxicity.Score',
  assay = "RNA")
