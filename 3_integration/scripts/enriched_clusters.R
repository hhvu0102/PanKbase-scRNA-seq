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
Sys.setenv(RETICULATE_PYTHON='/home/vthihong/miniconda3/envs/reticulate/bin/python')
reticulate::use_python('/home/vthihong/miniconda3/envs/reticulate/bin/python')
reticulate::use_condaenv('/home/vthihong/miniconda3/envs/reticulate')
leidenalg <- import("leidenalg")

option_list <- list(
  make_option(c("--sample"), action = 'store', type = 'character', help = '[Required] Sample SRR'),
  make_option(c("--outdir"), action = 'store', type = 'character', help = '[Required] Outdir')
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T, description = '')
opts <- parse_args(option_parser)

sample <- opts$sample
outdir <- opts$outdir

test_enrichment <- function(cluster, metadata = metadata, percentile = 0.9, feature = "nCount_RNA") {
    metadata <- metadata
    df <- data.frame(samples = unique(metadata$samples), q_percentile = NA, p = NA, x = NA, n = NA, p_val = NA)

    set.seed(1234)
    for (s in unique(metadata$samples)) { # per sample analyis
        #get 90th percentile of nUMIs across all cells in sample
        q_percentile <- quantile(metadata[metadata$samples == s, feature], percentile)

        # expected rate to observe cells of sample with 'feature' > percentile
        p <- 1 - percentile

        #in other clusters, what rate of cells in sample have 'feature' > percentile
        ## n is the number of cells in sample s and is in cluster
        n <- length(metadata[metadata$samples == s & metadata$seurat_clusters %in% cluster, feature])
        if (n > 50) {
            #if n > 50, downsample to n=50
            n <- 50
            #randomly sample 50 cells in sample s and is in cluster
            cells <- sample(metadata[metadata$samples == s & metadata$seurat_clusters %in% cluster, "rownames"], 50)
            #how many cells in the 50 cells sampled above has 'feature' > q0.9
            x <- length(which(metadata[metadata$samples == s &
                                       metadata$seurat_clusters %in% cluster &
                                       metadata$rownames %in% cells, feature] > q_percentile))
        } else {
            x <- length(which(metadata[metadata$samples == s &
                                       metadata$seurat_clusters %in% cluster, feature] > q_percentile))
        }

        df[df$samples == s, "q_percentile"] <- q_percentile
        df[df$samples == s, "p"] <- p
        df[df$samples == s, "x"] <- x
        df[df$samples == s, "n"] <- n

        #Are there more cells with very high coverage in a cluster for a sample than expected
        if (n > 0) {
            test <- binom.test(x, n, p, "greater")
            df[df$samples == s, "p_val"] <- test$p.value
        } else {
            df[df$samples == s, "p_val"] <- NA
        }
    }

    df$cluster <- cluster
    df$adj_pval <- p.adjust(df$p_val, method = "BH", n = length(df$p_val))
    df <- inner_join(df, distinct(metadata[, c("samples", "treatments")]))

    return(df)
}

metadata <- read.table("/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/doublet_like/202504_freeze/metadata.txt", header = T)

sample_info <- read.table("/nfs/turbo/umms-scjp-pank/4_integration/data/202503_freeze/20250410_meta_runs.txt", header = T)
treatment <- sample_info[sample_info$sample_id == sample, "treatments"]

rna <- readRDS(paste0("/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/doubletfinder_round2/nonDup_proteinCoding/pass-qc-nuclei-counts-with-doublets/rds/", sample, ".rna.rds"))

rna@meta.data$samples <- sample
rna@meta.data$rownames <- rownames(rna@meta.data)
rna@meta.data$treatments <- treatment

# at this point, only doublets marked at DoubletFinder round 1 had been removed
markers <- c("ENSG00000254647.7 (INS)",
             "ENSG00000115263.15 (GCG)", "ENSG00000157005.4 (SST)", "ENSG00000108849.8 (PPY)",
             "ENSG00000001626.17 (CFTR)", "ENSG00000204983.15 (PRSS1)", "ENSG00000113721.14 (PDGFRB)",
             "ENSG00000130300.9 (PLVAP)", "ENSG00000173372.17 (C1QA)")

p1 <- DimPlot(rna, label = T) + ggtitle(sample)
p2 <- DimPlot(rna, group.by = "doublet")
p <- FeaturePlot(rna, features = markers) + p1 + p2
png(paste0(outdir, "/", sample, "__", treatment, "__featureplots_postDFround1.png"), width = 15, height = 12, units = "in", res = 300)
print(p)
dev.off()

# remove doublets marked by DF round 2 and cells in clusters with rate > 65%
rna@meta.data$to_remove <- ifelse(rna@meta.data$barcodes %in%
                                  metadata[metadata$samples == sample & metadata$treatments == treatment, "barcodes"],
                                 0, 1)
rna <- subset(x = rna, subset = to_remove == 0)
rna <- RunUMAP(rna, dims=1:20, verbose = F)
rna <- FindNeighbors(rna, dims = 1:20, k.param = 20, verbose = F)
rna <- FindClusters(rna, resolution = 3, n.start = 100, verbose = F)
Idents(rna) <- as.factor(rna@meta.data$`RNA_snn_res.3`)
Idents(rna) <- factor(Idents(rna), levels=0:max(as.numeric(Idents(rna))))
rna$seurat_clusters <- factor(Idents(rna), levels=0:max(as.numeric(Idents(rna))))

## plot based on remained cells
p1 <- DimPlot(rna, label = T) + ggtitle(sample)
p <- FeaturePlot(rna, features = markers) + p1
png(paste0(outdir, "/", sample, "__", treatment, "__featureplots_postRmDoubletClust.png"), width = 15, height = 12, units = "in", res = 300)
print(p)
dev.off()

png(paste0(outdir, "/", sample, "__", treatment, "__dotplot_postRmDoubletClust.png"), width = 15, height = 6, units = "in", res = 300)
print(DotPlot(rna, features = markers) + coord_flip())
dev.off()

png(paste0(outdir, "/", sample, "__", treatment, "__dimplot_postRmDoubletClust.png"), width = 5, height = 4, units = "in", res = 300)
print(DimPlot(rna, label = T) + ggtitle(sample))
dev.off()

# enrichment using binominal test to identify clusters with high ratio (high means rate > 0.1) of cells whose nUMIs > q0.9 of the sample
i = 0
df <- as.data.frame(test_enrichment(i, metadata = rna@meta.data, percentile = 0.9))
for (i in 1:max(as.numeric(rna$seurat_clusters))) {
    a <- test_enrichment(i, metadata = rna@meta.data, percentile = 0.90)
    df <- rbind(df, a)
}
df$adj_pval <-  p.adjust(df$p_val, method = "bonferroni", n = length(df$p_val))
write.table(df, paste0(outdir, "/", sample, "__", treatment, "__nCount_binomTest.txt"), quote = F, row.names = F, sep = "\t")

i = 0
dfFeature <- as.data.frame(test_enrichment(i, metadata = rna@meta.data, percentile = 0.9, feature = "nFeature_RNA"))
for (i in 1:max(as.numeric(rna$seurat_clusters))) {
    a <- test_enrichment(i, metadata = rna@meta.data, percentile = 0.9, feature = "nFeature_RNA")
    dfFeature <- rbind(dfFeature, a)
}
dfFeature$adj_pval <-  p.adjust(dfFeature$p_val, method = "bonferroni", n = length(dfFeature$p_val))
write.table(dfFeature, paste0(outdir, "/", sample, "__", treatment, "__nFeature_binomTest.txt"), quote = F, row.names = F, sep = "\t")

png(paste0(outdir, "/", sample, "__", treatment, "__vlnplotCount_postRmDoubletClust.png"), width = 10, height = 3, units = "in", res = 300)
print(VlnPlot(rna, "nCount_RNA") + geom_hline(yintercept = df$q_percentile[1]))
dev.off()
png(paste0(outdir, "/", sample, "__", treatment, "__vlnplotFeature_postRmDoubletClust.png"), width = 10, height = 3, units = "in", res = 300)
print(VlnPlot(rna, "nFeature_RNA") + geom_hline(yintercept = dfFeature$q_percentile[1]))
dev.off()

data <- DotPlot(rna, features = markers)
a <- data$data[data$data$pct.exp > 50,]
b <- data.frame(table(a$id, a$features.plot)) %>% group_by(Var1) %>% summarise(n = sum(Freq))
c <- unique(b[b$n >= 2,]$Var1)
print(c)

c1 <- length(c)
l <- length(df[df$adj_pval < 0.05 & !is.na(df$p_val), "cluster"])

if (l > 0 & c1 > 0) { 
        clust <- df[df$adj_pval < 0.05 & !is.na(df$p_val), "cluster"]
        clust <- intersect(clust, c)
} else {
	clust <- "na"
}

print(l)
print(df[df$adj_pval < 0.05 & !is.na(df$p_val), "cluster"])
print(c1)
print(clust)

if (length(clust) > 0 & length(grep("na", clust)) == 0) { 
	cells <- rownames(rna@meta.data[rna@meta.data$seurat_clusters %in% clust,])
	png(paste0(outdir, "/", sample, "__", treatment, "__featureplots_enrichedClusters.png"), width = 15, height = 12, units = "in", res = 300)
	p1 <- DimPlot(rna, label = T, repel = F, cells = cells)
	p2 <- FeaturePlot(rna, features = markers, cells = cells)
	print(p2+p1)
	dev.off()


	cells <- rownames(rna@meta.data[!(rna@meta.data$seurat_clusters %in% clust),])
	png(paste0(outdir, "/", sample, "__", treatment, "__featureplots_nonenrichedClusters.png"), width = 15, height = 12, units = "in", res = 300)
	p1 <- DimPlot(rna, label = T, repel = F, cells = cells)
	p2 <- FeaturePlot(rna, features = markers, cells = cells)
	p2+p1
	dev.off()

} 


p <- DotPlot(rna, features = markers)
p1 <- p$data[p$data$id %in% clust,]
p1 <- p1[p1$avg.exp > 0,]
write.table(p1, paste0(outdir, "/", sample, "__", treatment, "__expGenes_enrichedClusters.txt"), quote = F, row.names = F, sep = "\t")

if (length(clust) > 0 &length(grep("na", clust)) == 0) {
	cells <- data.frame(x = rownames(rna@meta.data[rna@meta.data$seurat_clusters %in% clust,]))
	cells$x <- paste0(sample, "_", cells$x)
	write.table(cells, paste0(outdir, "/", sample, "__enrichedCells.txt"), quote = F, row.names = F, sep = "\t", col.names = F)
} else {
	cells <- "x"
	write.table(cells, paste0(outdir, "/", sample, "__enrichedCells.txt"), quote = F, row.names = F, sep = "\t", col.names = F)
}

if (length(clust) > 0 & length(grep("na", clust)) == 0) {
	rna@meta.data$enriched_cells <- ifelse(rownames(rna@meta.data) %in%
                             rownames(rna@meta.data[!(rna@meta.data$seurat_clusters %in% clust),]),
                             "False", "True")
	png(paste0(outdir, "/", sample, "__", treatment, "__featureplots_markingCells.png"), width = 15, height = 12, units = "in", res = 300)
	p1 <- DimPlot(rna, label = T, repel = F)
	p3 <- DimPlot(rna, label = F, repel = F, group.by = "enriched_cells")
	p2 <- FeaturePlot(rna, features = markers)
	print(p2+p1+p3)
	dev.off()
}
