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

# doublet finder results round 2
dblf <- read.table("/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/doubletfinder_round2/nonDup_proteinCoding/indivDblts_clusterRateAbv55pct.txt", header = F)
colnames(dblf) <- c("barcodes")

# HPAP
rna <- readRDS("/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/doubletfinder_round1/nonDup_proteinCoding/HPAP.rna.rds")
print("info prior to filtering based on samples (HPAP):")
print(rna)

rna@meta.data$to_remove <- ifelse(rna@meta.data$doubletfinder == "doublet", 1, 0)
rna@meta.data$to_remove <- ifelse(rna@meta.data$rownames %in% dblf$barcodes, 1, rna@meta.data$to_remove)

rna <- subset(x = rna, subset = to_remove == 0)

print("info prior to filtering based on min.cells and min.features (HPAP):")
print(rna)
print("number of samples:")
print(length(unique(rna@meta.data$samples)))

counts <- GetAssayData(rna, slot="counts", assay="RNA")

hpap <- CreateSeuratObject(counts=counts) 
print("info after filtering based on min.cells and min.features (HPAP):")
print(hpap)

# IIDP
rna <- readRDS("/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/doubletfinder_round1/nonDup_proteinCoding/IIDP.rna.rds")
print("info prior to filtering based on samples (IIDP):")
print(rna)

rna@meta.data$to_remove <- ifelse(rna@meta.data$doubletfinder == "doublet", 1, 0)
rna@meta.data$to_remove <- ifelse(rna@meta.data$rownames %in% dblf$barcodes, 1, rna@meta.data$to_remove)

rna <- subset(x = rna, subset = to_remove == 0)

print("info prior to filtering based on min.cells and min.features (IIDP):")
print(rna)
print("number of samples:")
print(length(unique(rna@meta.data$samples)))

counts <- GetAssayData(rna, slot="counts", assay="RNA")

iidp <- CreateSeuratObject(counts=counts) 
print("info after filtering based on min.cells and min.features (IIDP):")
print(iidp)

# Prodo
rna <- readRDS("/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/doubletfinder_round1/nonDup_proteinCoding/Prodo.rna.rds")
print("info prior to filtering based on samples (Prodo):")
print(rna)

rna@meta.data$to_remove <- ifelse(rna@meta.data$doubletfinder == "doublet", 1, 0)
rna@meta.data$to_remove <- ifelse(rna@meta.data$rownames %in% dblf$barcodes, 1, rna@meta.data$to_remove)

rna <- subset(x = rna, subset = to_remove == 0)

print("info prior to filtering based on min.cells and min.features (Prodo):")
print(rna)
print("number of samples:")
print(length(unique(rna@meta.data$samples)))

counts <- GetAssayData(rna, slot="counts", assay="RNA")

prodo <- CreateSeuratObject(counts=counts)
print("info after filtering based on min.cells and min.features (Prodo):")
print(prodo)

# merge data
print("START MERGING...")
merged_data <- merge(hpap, y=c(iidp, prodo), project='SeuratProj') #this may run into Error: vector::reserve because there are too many objects - to test

#normalize and scale
merged_data <- NormalizeData(merged_data, normalization.method = 'LogNormalize', scale.factor = 10000)
merged_data <- FindVariableFeatures(merged_data, selection.method = 'vst', nfeatures = 2000)
merged_data <- ScaleData(merged_data, verbose = FALSE) %>% 
  RunPCA(pc.genes = merged_data@var.genes, npcs = 50, verbose = FALSE)
#ElbowPlot(merged_data, ndims = 50)
merged_data <- RunUMAP(merged_data, reduction = 'pca', dims = 1:30)

saveRDS(merged_data, "/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/min.cells1_min.features1_rmIndivDblts_rmClustAbv55pct_merged_data.Rds")

# integration
print("START INTEGRATING")
meta <- read.table("/nfs/turbo/umms-scjp-pank/4_integration/data/202503_freeze/20250314_meta_runs.txt", header = T)

merged_data@meta.data$rownames <- rownames(merged_data@meta.data)
merged_data@meta.data$samples <- sub("(.*)-[^-]*$", "\\1", merged_data@meta.data$rownames)
merged_data@meta.data$barcodes <- sub(".*-", "", merged_data@meta.data$rownames)
merged_data@meta.data <- inner_join(merged_data@meta.data, meta, by = c("samples" = "srr")) # add metadata
rownames(merged_data@meta.data) <- merged_data@meta.data$rownames

#correct for sex, bmi, age, studies, treatments, donor IDs, and chemistry
merged_data <- RunHarmony(merged_data, c("sex", "bmi", "age", "chemistry", "study", "treatments", "source"), assay.use='RNA', plot_convergence = TRUE)

merged_data <- merged_data %>%
  RunUMAP(reduction = 'harmony', dims = 1:30) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:30) %>%
  FindClusters(algorithm = 4, resolution = 1, method = 'igraph') %>% #Algorithm 4 is Leiden clustering
  FindClusters(algorithm = 4, resolution = 0.8, method = 'igraph') %>%
  FindClusters(algorithm = 4, resolution = 0.6, method = 'igraph') %>% 
  FindClusters(algorithm = 4, resolution = 1.8, method = 'igraph')


saveRDS(merged_data, "/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/min.cells1_min.features1_rmIndivDblts_rmClustAbv55pct_harmonized_data.Rds")

