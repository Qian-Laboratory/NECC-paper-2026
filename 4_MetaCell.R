# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(metacell)

# Read NE (neuroendocrine) data and markers
NE <- readRDS("~/1_Project/NECC/data_final/Epi/NE.cell.rds")
marker <- readRDS("~/1_Project/NECC/data_final/Epi/NE.new.marker.rds")


# Collect variable genes from clusters
gene.all <- NULL
for(i in unique(marker$cluster)) {
    temp <- marker %>% filter(cluster == i) %>% head(500)
    gene <- rownames(temp)
    gene.all <- c(gene, gene.all)
}
var.genes <- unique(gene.all)

# Define genes to exclude (stress, hypoxia, immunoglobulin, ribosomal, mitochondrial, etc.)
hypoxia.genes <- c("VEGFA", "SLC2A1", "PGAM1", "ENO1", "LDHA", "TPI1", "P4HA1", "MRPS17",
                   "CDKN3", "ADM", "NDRG1", "TUBB6", "ALDOA", "MIF", "ACOT7")

ig_genes <- c(grep("^IGJ", var.genes, v = TRUE),
              grep("^IGH", var.genes, v = TRUE),
              grep("^IGK", var.genes, v = TRUE),
              grep("^IGL", var.genes, v = TRUE))

hgGenes <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
igGenes <- c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", "JCHAIN",
             "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC")

TRV.genes <- c(grep("^TRDV", var.genes, v = TRUE),
               grep("^TRAV", var.genes, v = TRUE),
               grep("^TRBV", var.genes, v = TRUE),
               grep("^TRGV", var.genes, v = TRUE))

# Combine all bad genes to remove
bad_genes <- unique(c(grep("^MT-", var.genes, v = TRUE),
                      grep("^RPL", var.genes, v = TRUE),
                      grep("^HIST", var.genes, v = TRUE),
                      grep("^RPS", var.genes, v = TRUE),
                      ig_genes, hgGenes, TRV.genes, igGenes, hypoxia.genes,
                      cc.genes$g2m.genes, cc.genes$s.genes))

var.genes <- var.genes[!var.genes %in% bad_genes]
length(var.genes)

# Initialize Metacell directories
if (!dir.exists("scdb")) dir.create("scdb")
scdb_init("scdb", force_reinit = TRUE)

if (!dir.exists("figs")) dir.create("figs")
scfigs_init("figs")

# Create gene set from variable genes
var.genes <- structure(rep(1:length(var.genes)), names = var.genes)
var.genes <- gset_new_gset(sets = var.genes, desc = "seurat variable genes")
scdb_add_gset("NE", var.genes)

# Convert Seurat to SingleCellExperiment and import as matrix
sce <- as.SingleCellExperiment(NE)
mat <- scm_import_sce_to_mat(sce)
scdb_add_mat("NE", mat)

# Build balanced kNN graph
mcell_add_cgraph_from_mat_bknn(mat_id = "NE",
                               gset_id = "NE",
                               graph_id = "NE_k100",
                               K = 100,
                               dsamp = FALSE)

# Co-clustering
mcell_coclust_from_graph_resamp(coc_id = "NE_n1000",
                                graph_id = "NE_k100",
                                min_mc_size = 20,
                                p_resamp = 0.75,
                                n_resamp = 1000)

# Generate initial metacells
mcell_mc_from_coclust_balanced(coc_id = "NE_n1000",
                               mat_id = "NE",
                               mc_id = "NE",
                               K = 30,
                               min_mc_size = 20,
                               alpha = 2)

# Filter and refine metacells
mcell_plot_outlier_heatmap(mc_id = "NE", mat_id = "NE", T_lfc = 3)
mcell_mc_split_filt(new_mc_id = "NE",
                    mc_id = "NE",
                    mat_id = "NE",
                    T_lfc = 3,
                    plot_mats = TRUE)

# Set metacell colors
mc_f <- scdb_mc("NE")
my.col <- c("darkgray", "burlywood1", "chocolate4", "orange", "red", "purple", "blue", "darkgoldenrod3", "cyan")
mc_f@colors <- colorRampPalette(my.col)(ncol(mc_f@mc_fp))
scdb_add_mc("NE", mc_f)

# Generate markers
mcell_gset_from_mc_markers(gset_id = "NE_markers", mc_id = "NE")

# Plot marker heatmaps
options(repr.plot.width = 20, repr.plot.height = 50)

# Single-cell level heatmap
mcell_mc_plot_marks(mc_id = "NE",
                    gset_id = "NE_markers",
                    mat_id = "NE",
                    plot_cells = TRUE)

# Metacell level heatmap with reordering
mcell_mc_plot_marks(mc_id = "NE",
                    gset_id = "NE_markers",
                    mat_id = "NE",
                    reorder_marks = TRUE,
                    plot_cells = TRUE)

# Extract gene expression folds
gene_folds <- mc@mc_fp

# Plot heatmap for selected genes
options(repr.plot.width = 18, repr.plot.height = 18)
mtx <- mc@mc_fp[NE.cells_heat_marks_gene_names, ]
pheatmap(mtx, cluster_rows = FALSE, cluster_cols = FALSE, scale = "column")

# Configure and plot 2D metacell projection
tgconfig::set_param("mcell_mc2d_height", 1000, "metacell")
tgconfig::set_param("mcell_mc2d_width", 1000, "metacell")
mcell_mc2d_plot(mc2d_id = "NE")

# Export metacell table
mcell_mc_export_tab(mc_id = "NE", gstat_id = "NE", mat_id = "NE", T_fold = 2)

# Hierarchical clustering of metacells
mc_hc <- mcell_mc_hclust_confu(mc_id = "NE", graph_id = "NE_k100")
mc_sup <- mcell_mc_hierarchy(mc_id = "NE", mc_hc = mc_hc, T_gap = 0.04)

# Plot metacell hierarchy
mcell_mc_plot_hierarchy(mc_id = "NE",
                        graph_id = "NE_k100",
                        mc_order = mc_hc$order,
                        sup_mc = mc_sup,
                        width = 2800,
                        heigh = 2000,
                        min_nmc = 2)

# Load saved metacell object
load("~/1_Project/NECC/script/new/metacell/scdb/mc.NE.Rda")
mc.NE <- object

# Calculate log2 fold change
lfp <- log2(mc.NE@mc_fp)
lfp <- log2(mc_f@mc_fp)

# Barplot for specific gene (CHGA) across metacells
barplot(lfp[, 1:10],
        col = mc_f@colors,
        las = 2,
        main = "CHGA",
        cex.main = 3,
        cex.axis = 1,
        ylab = "log2FC",
        xlab = "metacells")

# Function to create scatter plot comparing two genes
plt <- function(gene1, gene2, lfp, colors) {
    plot(lfp[gene1, ], lfp[gene2, ],
         pch = 21, cex = 3, bg = colors,
         xlab = gene1, ylab = gene2)
    text(lfp[gene1, ], lfp[gene2, ], colnames(lfp))
}

# Plot gene expression comparisons
options(repr.plot.width = 4, repr.plot.height = 4.5)
plt(gene1 = 'ASCL1', gene2 = 'ZEB1', lfp = lfp, colors = mc.NE@colors)
plt(gene1 = 'ASCL1', gene2 = 'NEUROD1', lfp = lfp, colors = mc.NE@colors)
