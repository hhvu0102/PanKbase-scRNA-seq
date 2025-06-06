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


options(future.globals.maxSize = 8000 * 1024^2)

harmonized_data <- readRDS("/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/202504_freeze/filtered_round1__min.cells1_min.features1_rmIndivDblts_rmClustAbv65pct_harmonized_data.Rds")
print(harmonized_data)
enriched_cells <- read.table("/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/doublet_like/202504_freeze/enriched_cells/allEnrichedClust.txt", header = F)
enriched_cells$V1 <- gsub("_", "-", enriched_cells$V1)

Idents(harmonized_data) <- as.factor(harmonized_data@meta.data$`RNA_snn_res.4`)
Idents(harmonized_data) <- factor(Idents(harmonized_data), levels=1:max(as.numeric(Idents(harmonized_data))))
harmonized_data$seurat_clusters <- factor(Idents(harmonized_data), levels=1:max(as.numeric(Idents(harmonized_data))))

harmonized_data$enriched_cells <- ifelse(rownames(harmonized_data@meta.data) %in% enriched_cells$V1, "True", "False")

markers <- c('INS', #Beta
             'GCG', #Alpha
             'SST', #Delta
             'PPY', #Gamma
             'GHRL', #Epsilon
             'CFTR', #Ductal
             'PRSS1',#Acinar
             'PDGFRB', #Stellate
             'PLVAP' #Endothelial
             )
data <- DotPlot(harmonized_data, features = markers, cluster.idents = TRUE, col.min = 0, scale = T)

a <- data$data[data$data$pct.exp > 50,]
b <- data.frame(table(a$id, a$features.plot)) %>% group_by(Var1) %>% summarise(n = sum(Freq))
cluster_data <- unique(b[b$n >= 2,]$Var1)
print(cluster_data)

metadata <- harmonized_data@meta.data

### for Fishers exact test

df <- data.frame(cluster = unique(harmonized_data$seurat_clusters),
                 ncells = NA,
                 a = NA, b = NA, c = NA, d = NA, p_val = NA)

for (s in unique(harmonized_data$seurat_clusters)) {
    df[df$cluster == s, "ncells"] <- length(metadata[metadata$seurat_clusters == s, "barcodes"])

    a <- length(metadata[metadata$seurat_clusters == s & metadata$enriched_cells == "True", "barcodes"])
    n <- length(metadata[metadata$seurat_clusters == s, "barcodes"])
    b <- n - a
    c <- length(metadata[metadata$enriched_cells == "True", "barcodes"]) - a
    d <- length(metadata[metadata$seurat_clusters != s, "barcodes"]) - c

    df[df$cluster == s, "a"] <- a
    df[df$cluster == s, "b"] <- b
    df[df$cluster == s, "c"] <- c
    df[df$cluster == s, "d"] <- d

    t <- matrix(c(a, c, b, d), ncol = 2, dimnames = list(nCount = c("n_enriched_cells", "n_nonenriched_cells"),
                                                        nCells = c("nCells_inCluster", "nCells_outCluster")))

    #Are there more cells with very high coverage in a cluster for a sample than expected
    if (n > 0) {
        test <- fisher.test(t(t), alternative = "greater")
        df[df$cluster == s, "p_val"] <- test$p.value
    } else {
        df[df$cluster == s, "p_val"] <- NA
    }
}
df$adj_pval <- p.adjust(df$p_val, method = "BH", n = length(df$p_val))
df$fold <- (df$a/(df$a+df$b)) / (df$c/(df$c+df$d))
df$rate <- df$a/df$b
res <- df[df$adj_pval < 0.05,]
print(res$cluster)

c1 <- c(res$cluster)
print(intersect(c1, cluster_data))

# these 2 steps are to remove doublet-like cells and clusters
harmonized_data$to_remove <- ifelse(harmonized_data$enriched_cells == "True", 1, 0)
harmonized_data$to_remove <- ifelse(harmonized_data$seurat_clusters %in% intersect(c1, cluster_data), 1, harmonized_data$to_remove)

table(harmonized_data$to_remove)

harmonized_data <- subset(x = harmonized_data, subset = to_remove == 0)

hto <- read.table("/nfs/turbo/umms-scjp-pank/2_IIDP/results/gencode_v39_private/GSE251912/results/singlets/HTO_barcode_maps.txt", header = T, sep = "\t")
tmp <- left_join(harmonized_data@meta.data, hto[, c("SRR", "barcode", "Treatment")], by = c("samples" = "SRR", "barcodes" = "barcode")) # add HTO demultiplexed results
treatments <- ifelse(is.na(tmp$Treatment), tmp$treatments, tmp$Treatment)
harmonized_data@meta.data$treatments <- treatments

harmonized_data <- RunHarmony(harmonized_data, c("sex", "bmi", "age", "chemistry", "study", "treatments", "source"), assay.use='RNA', plot_convergence = TRUE)

harmonized_data <- harmonized_data %>%
  RunUMAP(reduction = 'harmony', dims = 1:30) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:30) %>%
  FindClusters(algorithm = 4, resolution = 1, method = 'igraph') %>% #Algorithm 4 is Leiden clustering
  FindClusters(algorithm = 4, resolution = 1.8, method = 'igraph') %>%
  FindClusters(algorithm = 4, resolution = 3, method = 'igraph')

saveRDS(harmonized_data, "/nfs/turbo/umms-scjp-pank/4_integration/results/202503_freeze/202504_freeze/20250424_filtered_round2__min.cells1_min.features1_rmIndivDblts_rmClustAbv65pct_harmonized_data_enrichedCellClusters.Rds")

