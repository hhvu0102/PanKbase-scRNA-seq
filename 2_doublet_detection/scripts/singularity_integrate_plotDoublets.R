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


option_list <- list(
  make_option(c("--matrix"), action = 'store', type = 'character', help = '[Required] Matrix'),
  make_option(c("--features"), action = 'store', type = 'character', help = '[Required] Features'),
  make_option(c("--barcodes"), action = 'store', type = 'character', help = '[Required] Barcodes'),
  make_option(c("--resolution"), action = 'store', type = 'numeric', default = 0.1, help = '[Optional] Resolution to use in clustering (default: 0.1)'),
  make_option(c("--pcs"), action = 'store', type = 'numeric', default = 20, help = '[Optional] Number of top PCs to use in clustering (default: 20)'),
  make_option(c("--prefix"), action = 'store', type = 'character', default = 'seurat.', help = '[Optional] Prefix of output files'),
  make_option(c("--dfdir"), action = 'store', type = 'character', help = '[Required] Directory of doubletfinder results'),
  make_option(c("--threshold"), action = 'store', type = 'numeric', default = 20, help = '[Optional] Threshold for doublet rates to remove clusters (default: 20)')
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T, description = '')
opts <- parse_args(option_parser)

RNA_MTX <- opts$matrix # "/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/mm_merge/HPAP_merged.matrix.mtx"
RNA_FEATURES <- opts$features # "/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/mm_merge/HPAP_merged.features.tsv"
RNA_BARCODES <- opts$barcodes # "/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/mm_merge/HPAP_merged.barcodes.tsv"
dir <- opts$dfdir # /nfs/turbo/umms-scjp-pank/4_integration/data/202503_freeze/doubletfinder_round1/
thres <- opts$threshold
PREFIX <- opts$prefix

load_mm <- function(matrix_file, features_file, barcodes_file) {
    tmp <- as(Matrix::readMM(matrix_file), 'dgCMatrix')
    features <- read.table(features_file, as.is=T, sep='\t', head=F)
    features <- features$V2
    barcodes <- read.table(barcodes_file, as.is=T, head=F)[,1]
    dimnames(tmp) <- list(features, barcodes)
    return(tmp)
}

files <- list.files(paste0(dir, "/doubletfinder/"))
dblf <- data.frame("barcodes" = NA, "doubletfinder" = NA, "donor" = NA)
for (i in files) {
    donor <- gsub(".doubletfinder_assignments.txt", "", i)
    df <- read.table(paste0(dir, "/doubletfinder/", i), header = F)
    colnames(df) <- c("barcodes", "doubletfinder")
    df$donor <- donor
    dblf <- rbind(dblf, df)
}
dblf$barcodes <- paste0(dblf$donor, "-", dblf$barcodes)
dblf <- dblf[!is.na(dblf$doubletfinder),]

meta <- read.table("/nfs/turbo/umms-scjp-pank/4_integration/data/202503_freeze/20250314_meta_runs.txt", header = T)

mm <- load_mm(RNA_MTX, RNA_FEATURES, RNA_BARCODES)

rna <- CreateSeuratObject(counts = mm, min.cells=1, min.features=1, assay = "RNA", project='RNA')
rna <- NormalizeData(rna, verbose=F)
rna <- FindVariableFeatures(rna, selection.method='vst', nfeatures=2000, verbose=F)
rna <- ScaleData(rna, verbose=F)
rna <- RunPCA(rna, npcs=200, verbose=F)

rna@meta.data$rownames <- rownames(rna@meta.data)
rna@meta.data$samples <- sub("(.*)-[^-]*$", "\\1", rna@meta.data$rownames)
rna@meta.data$barcodes <- sub(".*-", "", rna@meta.data$rownames)
rna@meta.data <- inner_join(rna@meta.data, meta, by = c("samples" = "srr")) # add metadata
rna@meta.data <- inner_join(rna@meta.data, dblf, by = c("rownames" = "barcodes")) # add doublet status
rownames(rna@meta.data) <- rna@meta.data$rownames

if (PREFIX == "IIDP") {
	hto <- read.table("/nfs/turbo/umms-scjp-pank/2_IIDP/results/gencode_v39_private/GSE251912/results/singlets/HTO_barcode_maps.txt", header = T, sep = "\t")
	tmp <- left_join(rna@meta.data, hto[, c("SRR", "barcode", "Treatment")], by = c("samples" = "SRR", "barcodes" = "barcode")) # add HTO demultiplexed results
	treatments <- ifelse(is.na(tmp$Treatment), tmp$treatments, tmp$Treatment)
	rna@meta.data$treatments <- treatments
}
#correct for sex, bmi, age, studies, treatments and chemistry
rna <- RunHarmony(rna, c("sex", "bmi", "age", "chemistry", "study", "treatments"), assay.use='RNA', plot_convergence = TRUE)

rna <- rna %>%
  RunUMAP(reduction = 'harmony', dims = 1:30) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:30) %>%
  FindClusters(algorithm = 4, resolution = 1.8, method = 'igraph') #Algorithm 4 is Leiden clustering
saveRDS(rna, paste0("/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/doubletfinder_round1/nonDup_proteinCoding/singularity_", PREFIX, ".rna.rds"))
