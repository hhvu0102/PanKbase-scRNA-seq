#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(glue)
library(optparse)
library(Matrix)
library(Seurat)


option_list <- list(
  make_option(c("--matrix"), action = 'store', type = 'character', help = '[Required] Matrix for HTO'),
  make_option(c("--features"), action = 'store', type = 'character', help = '[Required] Features for HTO'),
  make_option(c("--barcodes"), action = 'store', type = 'character', help = '[Required] Barcodes for HTO'),
  make_option(c("--rna"), action = 'store', type = 'character', help = '[Required] Seurat object of GEX file'),
  make_option(c("--prefix"), action = 'store', type = 'character', default = 'seurat.', help = '[Required] Prefix of output files'),
  make_option(c("--outdir"), action = 'store', type = 'character', default = 'seurat.', help = '[Required] Output dir')
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T, description = '')
opts <- parse_args(option_parser)


load_HTO <- function(barcode.path, features.path, matrix.path) {
    mat.htos <- readMM(file = matrix.path)
    feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
    barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
    colnames(mat.htos) = barcode.names$V1
    rownames(mat.htos) = feature.names$V1

    return(mat.htos)
}

hto_dir = "/nfs/turbo/umms-scjp-pank/2_IIDP/results/gencode_v39_private/GSE251912/data/GEO/"
gex_dir = "/nfs/turbo/umms-scjp-pank/2_IIDP/results/gencode_v39_private/GSE251912/results/temp_rds/"


hto_bc <- opts$barcodes
hto_fea <- opts$features
hto_matrix <- opts$matrix
RNA <- opts$rna
PREFIX <- opts$prefix
OUTDIR <- opts$outdir

HTO_BC <- paste0(hto_dir, hto_bc)
HTO_FEATURE <- paste0(hto_dir, hto_fea)
HTO_MATRIX <- paste0(hto_dir, hto_matrix)

mat.htos <- load_HTO(HTO_BC, HTO_FEATURE, HTO_MATRIX)

rna <- readRDS(paste0(gex_dir, RNA))

joint.bcs <- intersect(colnames(mat.htos), colnames(rna))
mat.htos = mat.htos[rownames(mat.htos) != "unmapped", ]
mat.htos = mat.htos[rowSums(mat.htos) > 0, ]

# Subset RNA and HTO counts by joint cell barcodes
rna <- rna[, joint.bcs]
mat.htos <- as.matrix(mat.htos[, joint.bcs])

# Setup Seurat object
islet.hashtag <- CreateSeuratObject(counts = rna@assays$RNA@counts)
islet.hashtag <- NormalizeData(islet.hashtag)
islet.hashtag <- FindVariableFeatures(islet.hashtag, selection.method = "mean.var.plot")
islet.hashtag <- ScaleData(islet.hashtag, features = VariableFeatures(islet.hashtag))

# Add HTO data as a new assay independent from RNA
islet.hashtag[["HTO"]] <- CreateAssayObject(counts = mat.htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
islet.hashtag <- NormalizeData(islet.hashtag, assay = "HTO", normalization.method = "CLR")
islet.hashtag

islet.hashtag <- HTODemux(islet.hashtag, assay = "HTO", positive.quantile = 0.99)
table(islet.hashtag$HTO_classification.global)
# Group cells based on the max HTO signal
Idents(islet.hashtag) <- "HTO_maxID"

for (i in 1:nrow(mat.htos)) {
    p1 <- RidgePlot(islet.hashtag, assay = "HTO", features = rownames(islet.hashtag[["HTO"]])[i], ncol = 1)
    png(paste0(OUTDIR, PREFIX, rownames(islet.hashtag[["HTO"]])[i], '_ridge-plot.png'), height=8, width=18, units='in', res=300)
    print(p1)
    dev.off()
}
Idents(islet.hashtag) <- "HTO_classification.global"

png(paste0(OUTDIR, PREFIX, 'violin-plot.png'), height=3, width=9, units='in', res=300)
VlnPlot(islet.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()

saveRDS(islet.hashtag, paste0(OUTDIR, PREFIX, 'hashtag.Rds'))

# First, we will remove negative cells from the object
islet.hashtag.subset <- subset(islet.hashtag, idents = "Negative", invert = TRUE)

# Calculate a UMAP embedding of the HTO data
DefaultAssay(islet.hashtag.subset) <- "HTO"
islet.hashtag.subset <- ScaleData(islet.hashtag.subset, features = rownames(islet.hashtag.subset),
    verbose = FALSE)
islet.hashtag.subset <- RunPCA(islet.hashtag.subset, features = rownames(islet.hashtag.subset), approx = FALSE)
islet.hashtag.subset <- RunUMAP(islet.hashtag.subset, dims = 1:nrow(mat.htos))

png(paste0(OUTDIR, PREFIX, 'dimplot-withDoublets.png'), height=8, width=7, units='in', res=300)
DimPlot(islet.hashtag.subset)
dev.off()

# Extract the singlets
islet.singlet <- subset(islet.hashtag, idents = "Singlet")

# Select the top 1000 most variable features
islet.singlet <- FindVariableFeatures(islet.singlet, selection.method = "mean.var.plot")
islet.singlet <- ScaleData(islet.singlet, features = VariableFeatures(islet.singlet))
islet.singlet <- RunPCA(islet.singlet, features = VariableFeatures(islet.singlet))
islet.singlet <- FindNeighbors(islet.singlet, reduction = "pca", dims = 1:20)
islet.singlet <- FindClusters(islet.singlet, resolution = 0.6, verbose = FALSE)
islet.singlet <- RunUMAP(islet.singlet, reduction = "pca", dims = 1:20)

# Projecting singlet identities on UMAP visualization
png(paste0(OUTDIR, PREFIX, 'dimplot-HTOclass.png'), height=5, width=8, units='in', res=300)
DimPlot(islet.singlet, group.by = "HTO_classification")
dev.off()

saveRDS(islet.singlet, paste0(OUTDIR, PREFIX, 'singlet.Rds'))
