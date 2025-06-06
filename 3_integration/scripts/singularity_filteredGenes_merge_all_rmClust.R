library(Seurat)
library(ggplot2)
suppressPackageStartupMessages(library(dplyr))
library(glue)
suppressPackageStartupMessages(library(S4Vectors))
suppressPackageStartupMessages(library(GenomicRanges))
library(patchwork)
library(biomaRt)
library(harmony)
library(reticulate)
library(optparse)

Sys.setenv(RETICULATE_PYTHON='/root/miniconda3/envs/reticulate/bin/python')
reticulate::use_python('/root/miniconda3/envs/reticulate/bin/python')
reticulate::use_condaenv('/root/miniconda3/envs/reticulate')
leidenalg <- import("leidenalg")

# doublet finder results round 2
dblf <- read.table("/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/doubletfinder_round2/nonDup_proteinCoding/indivDblts_clusterRateAbv65pct.txt", header = F) # this file is a combination of doublets markes at both rounds and in cluster with doublet rates > 65%
colnames(dblf) <- c("barcodes")

# HPAP
rna <- readRDS("/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/doubletfinder_round1/nonDup_proteinCoding/singularity_HPAP.rna.rds")
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
rna <- readRDS("/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/doubletfinder_round1/nonDup_proteinCoding/singularity_IIDP.rna.rds")
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
rna <- readRDS("/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/doubletfinder_round1/nonDup_proteinCoding/singularity_Prodo.rna.rds")
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

prodo <- CreateSeuratObject(counts=counts) #, min.cells = ncol(rna)*0.01/100, min.features = nrow(rna)*3/100)
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

saveRDS(merged_data, "/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/202504_freeze/min.cells1_min.features1_rmIndivDblts_rmClustAbv55pct_merged_data.Rds")

# integration
print("START INTEGRATING")
meta <- read.table("/nfs/turbo/umms-scjp-pank/4_integration/data/202503_freeze/20250410_meta_runs.txt", header = T)

rna@meta.data$rownames <- rownames(rna@meta.data)
rna@meta.data$samples <- sub("(.*)-[^-]*$", "\\1", rna@meta.data$rownames)
rna@meta.data$barcodes <- sub(".*-", "", rna@meta.data$rownames)
rna@meta.data <- inner_join(rna@meta.data, meta, by = c("samples" = "sample_id")) # add metadata
rna@meta.data <- inner_join(rna@meta.data, dblf, by = c("rownames" = "barcodes")) # add doublet status
rownames(rna@meta.data) <- rna@meta.data$rownames

hto <- read.table("/nfs/turbo/umms-scjp-pank/2_IIDP/results/gencode_v39_private/GSE251912/results/singlets/HTO_barcode_maps.txt", header = T, sep = "\t")
tmp <- left_join(rna@meta.data, hto[, c("SRR", "barcode", "Treatment")], by = c("samples" = "SRR", "barcodes" = "barcode")) # add HTO demultiplexed results
treatments <- ifelse(is.na(tmp$Treatment), tmp$treatments, tmp$Treatment)
rna@meta.data$treatments <- treatments

#correct for sex, bmi, age, studies, treatments, donor IDs, and chemistry
merged_data <- RunHarmony(merged_data, c("sex", "bmi", "age", "chemistry", "study", "treatments", "source"), assay.use='RNA', plot_convergence = TRUE)

merged_data <- merged_data %>%
  RunUMAP(reduction = 'harmony', dims = 1:30) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:30) %>%
  FindClusters(algorithm = 4, resolution = 1, method = 'igraph') %>% #Algorithm 4 is Leiden clustering
  FindClusters(algorithm = 4, resolution = 1.8, method = 'igraph')


saveRDS(merged_data, "/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/202504_freeze/min.cells1_min.features1_rmIndivDblts_rmClustAbv65pct_harmonized_data.Rds")

