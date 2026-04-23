library(SCEVAN)

library(EnsDb.Hsapiens.v86)
load("~/1_Project/NECC/scevan/sysdata.rda")
source("~/1_Project/NECC/scevan/SCEVAN.R")


merge.clean <- readRDS("~/1_Project/NECC/data_final/all/merge_integrated_celltype.rds")
obj <- subset(merge.clean , Patient %in% c("NECC3"))

DefaultAssay(obj) <- 'RNA'
count  <- obj$RNA@counts

count_mtx <- count
sample= "NECC3"
par_cores = 10
norm_cell = NULL
SUBCLONES = TRUE
beta_vega = 0.5
ClonalCN = TRUE
plotTree = TRUE
AdditionalGeneSets = NULL
SCEVANsignatures = TRUE
organism = "human"

dir.create(file.path("/SCEVAN/output/NECC3"), showWarnings = FALSE)

start_time <- Sys.time()
norm_cell = NULL
normalNotKnown <- FALSE

res_proc <- preprocessingMtx(count_mtx,sample, par_cores=par_cores, findConfident = normalNotKnown,
                             AdditionalGeneSets = AdditionalGeneSets, SCEVANsignatures = SCEVANsignatures, organism = organism)
res_class <- classifyTumorCells(res_proc$count_mtx_norm, res_proc$count_mtx_annot, sample, par_cores=par_cores,
                                ground_truth = NULL,  norm_cell_names = norm_cell, SEGMENTATION_CLASS = TRUE, SMOOTH = TRUE, beta_vega = beta_vega)
print(paste("found", length(res_class$tum_cells), "tumor cells"))

classDf <- data.frame(class = rep("filtered", length(colnames(res_proc$count_mtx))), row.names = colnames(res_proc$count_mtx))
  classDf[colnames(res_class$CNAmat)[-(1:3)], "class"] <- "normal"
  classDf[res_class$tum_cells, "class"] <- "tumor"
  classDf[res_class$confidentNormal, "confidentNormal"] <- "yes"
  
  end_time<- Sys.time()

  print(paste("time classify tumor cells: ", end_time -start_time))

write.table(classDf,"/SCEVAN/output/NECC3/NECC3_class.Df.csv")

if(ClonalCN) getClonalCNProfile(res_class, res_proc, sample, par_cores, organism = organism)
 mtx_vega <- segmTumorMatrix(res_proc, res_class, sample, par_cores, beta_vega)
 

if (SUBCLONES) {
    res_subclones <- subcloneAnalysisPipeline(res_proc$count_mtx, res_class, res_proc,mtx_vega, sample, par_cores, classDf, beta_vega, plotTree, organism)
    #res_subclones <- subcloneAnalysisPipeline(count_mtx, res_class, res_proc,mtx_vega, sample, par_cores, classDf, 3, plotTree)
    FOUND_SUBCLONES <- res_subclones$FOUND_SUBCLONES
    classDf <- res_subclones$classDf
  }else{
    FOUND_SUBCLONES <- FALSE
  }
  
  #if(!FOUND_SUBCLONES) plotCNAlineOnlyTumor(sample) getClonalCNProfile(sample,)
  
  if(!FOUND_SUBCLONES) plotCNclonal(sample,ClonalCN, organism)


#save CNA matrix
#CNAmtx <- res_class$CNAmat[,-c(1,2,3)]
#save(CNAmtx, file = paste0("./output/",sample,"_CNAmtx.RData"))

#save annotated matrix
count_mtx_annot <- res_proc$count_mtx_annot
save(count_mtx_annot, file = paste0("./output/",sample,"_count_mtx_annot.RData"))


#remove intermediate files
mtx_vega_files <- list.files(path = "./output/", pattern = "_mtx_vega")
sapply(mtx_vega_files, function(x) file.remove(paste0("./output/",x)))

res_subclones <- subclonesTumorCells(res_class$tum_cells,
                                   res_class$CNAmat,
                                   sample,
                                   par_cores,
                                   beta_vega,
                                   res_proc,
                                   NULL,
                                   mtx_vega,
                                   organism = organism)
