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
Sys.setenv(RETICULATE_PYTHON='/root/miniconda3/envs/reticulate/bin/python')
reticulate::use_python('/root/miniconda3/envs/reticulate/bin/python')
reticulate::use_condaenv('/root/miniconda3/envs/reticulate')
leidenalg <- import("leidenalg")

harmonized_data <- readRDS("/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/202504_freeze/min.cells1_min.features1_rmIndivDblts_rmClustAbv65pct_harmonized_data.Rds")
print(harmonized_data)
metrics <- readRDS("/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/202504_freeze/min.cells1_min.features1_rmIndivDblts_rmClustAbv65pct_metrics.Rds")
harmonized_data@meta.data <- inner_join(harmonized_data@meta.data, metrics,
                                        by = c("barcodes" = "barcode", "samples" = "samples")) # add metadata
rownames(harmonized_data@meta.data) <- harmonized_data@meta.data$rownames
dim(harmonized_data@meta.data)

Idents(harmonized_data) <- as.factor(harmonized_data@meta.data$`RNA_snn_res.1.8`)
Idents(harmonized_data) <- factor(Idents(harmonized_data), levels=1:max(as.numeric(Idents(harmonized_data))))
harmonized_data$seurat_clusters <- factor(Idents(harmonized_data), levels=1:max(as.numeric(Idents(harmonized_data))))

harmonized_data@meta.data$to_remove <- ifelse(harmonized_data@meta.data$`RNA_snn_res.1.8` %in% c(29, 45, 48, 49, 52, 53), 1, 0) # remove clusters that have significantly different nFeatures and nCounts

harmonized_data <- subset(x = harmonized_data, subset = to_remove == 0)

harmonized_data <- RunHarmony(harmonized_data, c("sex", "bmi", "age", "chemistry", "study", "treatments", "source"), assay.use='RNA', plot_convergence = TRUE)

harmonized_data <- harmonized_data %>%
  RunUMAP(reduction = 'harmony', dims = 1:30) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:30) %>%
  FindClusters(algorithm = 4, resolution = 1, method = 'igraph') %>% #Algorithm 4 is Leiden clustering
  FindClusters(algorithm = 4, resolution = 1.8, method = 'igraph') %>%
  FindClusters(algorithm = 4, resolution = 3, method = 'igraph')

saveRDS(harmonized_data, "/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/202504_freeze/filtered_round1__min.cells1_min.features1_rmIndivDblts_rmClustAbv65pct_harmonized_data.Rds")

