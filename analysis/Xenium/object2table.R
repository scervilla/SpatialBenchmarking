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

samples <- c("Bladder", "Breast", "Lung", "Colorectal", "Lymphoma", "Ovarian")


compute_counts <- function(counts, genes) {
  list(
    count = colSums(counts[genes, , drop = FALSE]),
    feature = colSums(counts[genes, , drop = FALSE] > 0)
  )
}

for (sample in samples) {
  print(sample)
  
  cosmx <- readRDS(file.path(path_seurat, paste0(sample, "_alignment_cosmx.rds")))
  xenium <- readRDS(file.path(path_seurat, paste0(sample, "_alignment_xenium.rds")))
  genes_intersect <- intersect(rownames(cosmx), rownames(xenium))
  
  # COSMX counts
  cosmx_counts <- cosmx@assays$RNA$counts
  cosmx$nCount_intersect <- compute_counts(cosmx_counts, genes_intersect)$count
  cosmx$nFeature_intersect <- compute_counts(cosmx_counts, genes_intersect)$feature
  cosmx$nCount_intersect_all <- compute_counts(cosmx_counts, int_cosmx)$count
  cosmx$nFeature_intersect_all <- compute_counts(cosmx_counts, int_cosmx)$feature
  cosmx$nCount_specific <- compute_counts(cosmx_counts, setdiff(rownames(cosmx), genes_intersect))$count
  cosmx$nFeature_specific <- compute_counts(cosmx_counts, setdiff(rownames(cosmx), genes_intersect))$feature
  cosmx$nCount_total <- colSums(cosmx_counts)
  cosmx$nFeature_total <- colSums(cosmx_counts > 0)
  
  # XENIUM counts
  xenium_counts <- xenium@assays$Xenium$counts
  xenium$nCount_intersect <- compute_counts(xenium_counts, genes_intersect)$count
  xenium$nFeature_intersect <- compute_counts(xenium_counts, genes_intersect)$feature
  xenium$nCount_intersect_all <- compute_counts(xenium_counts, int_cosmx)$count
  xenium$nFeature_intersect_all <- compute_counts(xenium_counts, int_cosmx)$feature
  xenium$nCount_specific <- compute_counts(xenium_counts, setdiff(rownames(xenium), genes_intersect))$count
  xenium$nFeature_specific <- compute_counts(xenium_counts, setdiff(rownames(xenium), genes_intersect))$feature
  xenium$nCount_total <- colSums(xenium_counts)
  xenium$nFeature_total <- colSums(xenium_counts > 0)
  
  # Metadata
  df_cosmx <- cosmx@meta.data %>%
    select(fov, Area.um2, seurat_clusters) %>%
    mutate(
      nCount_intersect = cosmx$nCount_intersect,
      nFeature_intersect = cosmx$nFeature_intersect,
      nCount_specific = cosmx$nCount_specific,
      nFeature_specific = cosmx$nFeature_specific,
      nCount_total = cosmx$nCount_total,
      nFeature_total = cosmx$nFeature_total,
      nCount_intersect_all = cosmx$nCount_intersect_all,
      nFeature_intersect_all = cosmx$nFeature_intersect_all,
      nucleus_area = NA,
      x_centroid = cosmx@images$fov@boundaries$centroids@coords[,1],
      y_centroid = cosmx@images$fov@boundaries$centroids@coords[,2],
      platform = "cosmx",
      ct = readRDS( paste0(path_annotation,"predictions_", sample, "_cosmx_hpca_final.rds"))$pruned.labels
    )
  
  df_xenium <- xenium@meta.data %>%
    select(fov, cell_area, seurat_clusters) %>%
    rename(Area.um2 = cell_area) %>%
    mutate(
      nCount_intersect = xenium$nCount_intersect,
      nFeature_intersect = xenium$nFeature_intersect,
      nCount_specific = xenium$nCount_specific,
      nFeature_specific = xenium$nFeature_specific,
      nCount_total = xenium$nCount_total,
      nFeature_total = xenium$nFeature_total,
      nCount_intersect_all = xenium$nCount_intersect_all,
      nFeature_intersect_all = xenium$nFeature_intersect_all,
      nucleus_area = xenium$nucleus_area,
      x_centroid = xenium@images$fov@boundaries$centroids@coords[,1],
      y_centroid = xenium@images$fov@boundaries$centroids@coords[,2],
      platform = "xenium",
      ct = readRDS(paste0(path_annotation, "predictions_", sample, "_xenium_hpca_final.rds"))$pruned.labels
    )
  
  
  
  df <- bind_rows(
    df_cosmx %>% mutate(fov = as.character(fov)),
    df_xenium %>% mutate(fov = as.character(fov))
  ) %>% mutate(sample = sample)
  
  write.table(df, paste0("data/output/df/df_", sample, "_xenium_cosmx.txt"))
}
