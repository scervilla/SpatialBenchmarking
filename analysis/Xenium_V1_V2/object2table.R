library(Seurat)
library(dplyr)
library(ggplot2)
library(readxl)

# Paths
path_panels <- "data/raw/gene_panel/"
path_seurat <- "data/objects/seurat_aligned/"
path_annotation <- "data/output/SingleR/"

# Load gene panels
genes_cosmx <- read_xlsx(file.path(path_panels, "LBL-11178-05-Human-Universal-Cell-Characterization-Panel-Gene-Target-List.xlsx"), sheet = 2, skip = 1) %>% head(1000) %>% pull(`Display Name`)
genes_v1 <- read.csv(file.path(path_panels, "Xenium_hMulti_v1_metadata.csv")) %>% pull(Gene)
genes_v2 <- read.csv(file.path(path_panels, "XeniumPrimeHuman5Kpan_tissue_pathways_metadata.csv")) %>% pull(gene_name)


int_v1_cosmx <- intersect(genes_v1, genes_cosmx)
int_v2_cosmx <- intersect(genes_v2, genes_cosmx)
int_v1_v2 <- intersect(genes_v1, genes_v2)
int_cosmx <- intersect(int_v1_v2, genes_cosmx)

# Algnied cells
for (sample in samples) {
  for (version in versions) {
    print(sample)
    
    object <- readRDS(paste0( paste0(path_seurat, sample, "_", version, "_aligned.rds")))
    
    if (version == "V1") {
      intersect_cosmx <- compute_counts(object@assays$Xenium$counts, int_v1_cosmx)
    } else {
      intersect_cosmx <- compute_counts(object@assays$Xenium$counts, int_v2_cosmx)
    }
    
    intersect_v1v2 <- compute_counts(object@assays$Xenium$counts, int_v1_v2)
    intersect_all <- compute_counts(object@assays$Xenium$counts, int_cosmx)
    total <- compute_counts(object@assays$Xenium$counts, rownames(object))
    
    object$nCount_intersect_cosmx <- intersect_cosmx$count
    object$nFeature_intersect_cosmx <- intersect_cosmx$feature
    object$nCount_intersect_v1v2 <- intersect_v1v2$count
    object$nFeature_intersect_v1v2 <- intersect_v1v2$feature
    object$nCount_intersect_all <- intersect_all$count
    object$nFeature_intersect_all <- intersect_all$feature
    object$nCount_total <- total$count
    object$nFeature_total <- total$feature
    
    object@meta.data[, c("x", "y")] <- object@images[[1]]$centroids@coords
    
    ct_labels <- readRDS(paste0(path_annotation, "predictions_", sample, "_", version, "_aligned_final.rds"))$pruned.labels
    
    df_xenium <- object@meta.data %>%
      select(nCount_intersect_cosmx, nFeature_intersect_cosmx,
             nCount_intersect_v1v2, nFeature_intersect_v1v2,
             nCount_intersect_all, nFeature_intersect_all,
             nCount_total, nFeature_total,
             nCount_ControlProbe, nFeature_ControlProbe,
             cell_area, nucleus_area, x_centroid, y_centroid, x, y) %>%
      mutate(
        ct = ct_labels,
        platform = version,
        sample = sample
      )
    
    
    write.table(df_xenium,  paste0("data/output/df/df_", sample, "_v1_v2.txt"))
  }
}

# Complete list of cells
for (sample in samples) {
  for (version in versions) {
    print(sample)
    
    object <- readRDS(paste0( paste0(path_seurat, sample, "_", version, "_aligned.rds")))
    
    if (version == "V1") {
      intersect_cosmx <- compute_counts(object@assays$Xenium$counts, int_v1_cosmx)
    } else {
      intersect_cosmx <- compute_counts(object@assays$Xenium$counts, int_v2_cosmx)
    }
    
    intersect_v1v2 <- compute_counts(object@assays$Xenium$counts, int_v1_v2)
    intersect_all <- compute_counts(object@assays$Xenium$counts, int_cosmx)
    total <- compute_counts(object@assays$Xenium$counts, rownames(object))
    
    object$nCount_intersect_cosmx <- intersect_cosmx$count
    object$nFeature_intersect_cosmx <- intersect_cosmx$feature
    object$nCount_intersect_v1v2 <- intersect_v1v2$count
    object$nFeature_intersect_v1v2 <- intersect_v1v2$feature
    object$nCount_intersect_all <- intersect_all$count
    object$nFeature_intersect_all <- intersect_all$feature
    object$nCount_total <- total$count
    object$nFeature_total <- total$feature
    
    object@meta.data[, c("x", "y")] <- object@images[[1]]$centroids@coords
    
    df_xenium <- object@meta.data %>%
      select(nCount_intersect_cosmx, nFeature_intersect_cosmx,
             nCount_intersect_v1v2, nFeature_intersect_v1v2,
             nCount_intersect_all, nFeature_intersect_all,
             nCount_total, nFeature_total,
             cell_area, nucleus_area, x_centroid, y_centroid, x, y) %>%
      mutate(
        platform = version,
        sample = sample
      )
    
    write.table(df_xenium,  paste0("data/output/df/df_", sample, "_v1_v2_whole.txt"))
  }
}