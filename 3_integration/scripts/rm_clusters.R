library(Seurat)
library(ggplot2)
library(dplyr)
library(glue)
library(S4Vectors)
library(GenomicRanges)
library(patchwork)
library(biomaRt)
library(harmony)
library(reticulate)
Sys.setenv(RETICULATE_PYTHON='/home/vthihong/miniconda3/envs/reticulate/bin/python')
reticulate::use_python('/home/vthihong/miniconda3/envs/reticulate/bin/python')
reticulate::use_condaenv('/home/vthihong/miniconda3/envs/reticulate')
leidenalg <- import("leidenalg")
library(rtracklayer)
library(plyranges)

options(future.globals.maxSize = 8000 * 1024^2)

harmonized_data <- readRDS("/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/min.cells1_min.features1_rmIndivDblts_rmClustAbv55pct_harmonized_data.Rds")
metrics <- readRDS("/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/min.cells1_min.features1_rmIndivDblts_metrics.Rds")
harmonized_data@meta.data <- inner_join(harmonized_data@meta.data, metrics,
                                        by = c("barcodes" = "barcode", "samples" = "samples")) # add metadata
rownames(harmonized_data@meta.data) <- harmonized_data@meta.data$rownames

Idents(harmonized_data) <- as.factor(harmonized_data@meta.data$`RNA_snn_res.1.8`)
Idents(harmonized_data) <- factor(Idents(harmonized_data), levels=1:max(as.numeric(Idents(harmonized_data))))
harmonized_data$seurat_clusters <- factor(Idents(harmonized_data), levels=1:max(as.numeric(Idents(harmonized_data))))

harmonized_data@meta.data$to_remove <- ifelse(harmonized_data@meta.data$`RNA_snn_res.1.8` %in% c(20, 36, 37, 46, 48, 49, 50), 1, 0) # remove clusters that have significantly different nFeatures and nCounts

harmonized_data <- subset(x = harmonized_data, subset = to_remove == 0)

harmonized_data <- RunHarmony(harmonized_data, c("sex", "bmi", "age", "chemistry", "study", "treatments", "source"), assay.use='RNA', plot_convergence = TRUE)

harmonized_data <- harmonized_data %>%
  RunUMAP(reduction = 'harmony', dims = 1:30) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:30) %>%
  FindClusters(algorithm = 4, resolution = 1, method = 'igraph') %>% #Algorithm 4 is Leiden clustering
  FindClusters(algorithm = 4, resolution = 0.8, method = 'igraph') %>%
  FindClusters(algorithm = 4, resolution = 0.6, method = 'igraph') %>% 
  FindClusters(algorithm = 4, resolution = 1.8, method = 'igraph')

saveRDS(harmonized_data, "/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/filtered_round1__min.cells1_min.features1_rmIndivDblts_rmClustAbv55pct_harmonized_data.Rds")

