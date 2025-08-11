
library(dplyr)
library(Seurat)
path_objects <- "data/objects/"
path_panels <- "data/raw/gene_panel/"
genes_v1 <- read.csv(paste0(path_panels, "Xenium_hMulti_v1_metadata.csv")) %>% pull(Gene)
genes_v2 <- read.csv(paste0(path_panels, "XeniumPrimeHuman5Kpan_tissue_pathways_metadata.csv")) %>% pull(gene_name)

intersect_genes <- intersect(genes_v1, genes_v2)

df <- data.frame()
for (sample in c("Breast", "Lung", "Colorectal", "Lymphoma", "Ovarian")) {
  print(sample)
  for (platform in c("V1", "V2")) {
    print(sample)
    object <- readRDS(paste0(path_objects, "/seurat_aligned/", sample, "_", platform, "_aligned.rds"))
    
    gex <- sum(object@assays$Xenium$counts)
    negativex <- sum(object@assays$ControlProbe$counts)
    
    gex_int <- sum(object@assays$Xenium$counts[intersect_genes, ])
    
    df <- rbind(df, data.frame(
      neg_rate1 = (negativex / (negativex + gex)) * 100,
      neg_rate2 = (negativex / (negativex + gex_int)) * 100,
      neg_rate3 = (negativex / (negativex + gex)) *
        (nrow(object@assays$Xenium$counts) / nrow(object$ControlProbe$counts)) * 100,
      platform = platform,
      sample = sample
    ))
  }
}

write.table(df, "data/output/neg_rate/neg_rate_v1v2.txt")