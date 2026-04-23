# Load required libraries
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(copykat)

# Read merged and integrated Seurat object
merge.clean <- readRDS("~/1_Project/NECC/data_final/all/merge_integrated_celltype.rds")

# Subset reference immune cells
ref <- subset(merge.clean, celltype.new %in% c("CD4T", "CD8T", "Macrophage", "Monocyte", "DC", "B", "Mast"))

# Subset epithelial tumor cells
Epi <- subset(merge.clean, celltype.new %in% c("Epi.Tumor"))

# Merge reference and epithelial cells
merge <- merge(ref, Epi)

# Further subset to include only CIN, SCC, ADC, NECC samples
merge.copykat <- subset(merge, class.new %in% c("CIN", "SCC", "ADC", "NECC"))

# Run copykat analysis for each sample
for (i in unique(merge.copykat$sample)) {
    print(i)
    
    # Subset data for current sample
    sample.copykat <- subset(merge.copykat, sample == i)
    
    # Define normal cell types as reference
    normal_types <- c("CD4T", "CD8T", "Macrophage", "Monocyte", "DC", "B", "Mast")
    
    # Identify normal and cancer cells
    cells.nocancer <- WhichCells(sample.copykat, expression = celltype.new %in% normal_types, invert = FALSE)
    cells.cancer <- WhichCells(sample.copykat, expression = celltype.new %in% c("Epi.Tumor"), invert = FALSE)
    
    # Subset to include both cancer and normal cells
    sample.copykat <- subset(sample.copykat, cells = c(cells.cancer, cells.nocancer))
    
    print(sample.copykat)
    
    # Extract raw count matrix
    exp.rawdata <- sample.copykat@assays$RNA@counts
    
    print("start copykat")
    
    # Run copykat to infer copy number alterations
    sample.copykat.test <- copykat(rawmat = exp.rawdata,
                                   id.type = "S",
                                   ngene.chr = 5,
                                   win.size = 25,
                                   KS.cut = 0.1,
                                   sam.name = i,
                                   distance = "euclidean",
                                   norm.cell.names = cells.nocancer,
                                   output.seg = "FLASE",
                                   LOW.DR = 0.05,
                                   UP.DR = 0.2,
                                   cell.line = "no",
                                   plot.genes = "FALSE",
                                   genome = "hg20",
                                   n.cores = 10)
    
    # Format prediction results
    sample.pred.test <- tibble(sample.copykat.test$prediction) %>%
        dplyr::mutate(batch = i)
    
    # Merge with cell type metadata
    df <- sample.copykat@meta.data %>%
        rownames_to_column("cell.names") %>%
        select(cell.names, celltype.new)
    
    sample.pred.test <- df %>% left_join(sample.pred.test, by = "cell.names")
    
    print("finish copykat")
    
    # Save prediction results
    write.csv(sample.pred.test, paste0("~/NECC/data_new/copykat/copykat", "_", i, "_pred.test.csv"))
    print("saved")
}

# Function to process patient CNA data for chromosome 1
process_patient_cna <- function(data, patient_id) {
    # Filter for chromosome 1 data
    data_chr1 <- filter(data, chromosome_name == "1")
    
    # Extract expression matrix (remove metadata columns)
    expr_mat <- data_chr1[, c(-1:-5)][, -2]
    expr_mat <- expr_mat %>% column_to_rownames("hgnc_symbol")
    
    # Calculate mean CNA per gene across all cells
    cna_mean <- rowMeans(expr_mat) %>% data.frame()
    colnames(cna_mean) <- "CNA.mean"
    
    # Extract basic information (chromosome, gene, band, etc.)
    info <- data_chr1[, 2:7]
    
    # Merge info with mean CNA values
    result <- cbind(info, cna_mean)
    rownames(result) <- NULL
    result$Patient <- patient_id
    
    # Summarize by chromosomal band and get top 5 with lowest CNA mean
    result <- result %>%
        group_by(band) %>%
        summarize(CNA.mean = mean(CNA.mean)) %>%
        arrange(CNA.mean) %>%
        head(5)
    
    result <- result %>% mutate(Patient = patient_id)
    
    return(result)
}

# Define sample names
sample.name <- c("NECC1", "NECC2", "NECC3", "NECC4", "NECC5", "NECC6", "NECC7",
                 "CIN1", "CIN2", "CIN3", "CIN4")

# Initialize list to store combined results
chr1.df.list <- list()

# Loop through each sample to process CNA data
for (sample in sample.name) {
    # Construct variable name dynamically
    var_name <- paste0(sample, "_copykat_CNA_raw_results_gene_by_cell")
    
    # Retrieve data object by name
    data <- get(var_name)
    
    # Process the data
    df_temp <- process_patient_cna(data, sample)
    
    # Create individual variable (e.g., df.NECC7, df.CIN1, etc.)
    assign(paste0("df.", sample), df_temp)
    
    # Store in list for combined dataframe
    chr1.df.list[[sample]] <- df_temp
}

# Combine all data frames into one
chr1.df <- do.call(rbind, chr1.df.list)

options(repr.plot.width =4, repr.plot.height =7)
my_col <- c(
  "#9FA3A8", "#E0D4CA", "#5F3D69", "#C5DEBA", "#58A4C3", "#E4C755", "#F7F398",
  "#E95C59", "#E59CC4", "#AB3282", "#23452F", "#BD956A", "#8C549C", "#585658",
  "#AA9A59", "#E63863", "#E39A35", "#C1E6F3", "#6778AE", "#91D0BE", "#B53E2B",
  "#712820", "#DCC1DD", "#CCE0F5", "#CCC9E6", "#625D9E", "#68A180", "#3A6963",
  "#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3", "#476D87",
  "#968175"
)

chr1.df  %>% ggplot(aes(band,
                     CNA.mean,fill = Patient))+
    geom_col(color = 'black',width = 0.8,size = 0.5)+
    coord_flip()+
    theme_bw()+
theme(panel.border = element_rect(fill=NA,color="black", size=0.8, linetype="solid"),
      panel.grid =element_blank(),
          axis.text.x = element_text(color = 'black',size = 11),
          axis.text.y = element_text(color = 'black',size = 11))+
scale_fill_manual(values = my_col)+
ylab("CNA")+xlab("chr1 band")+
#scale_y_continuous(expand = c(0,0),limits = c(0,0.17))  +
    ggtitle("CNA Chr1(hg19)")
dev.off()
