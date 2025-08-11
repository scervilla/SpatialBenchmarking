library(dplyr)
library(Seurat)

path_objects <- "data/objects/"
path_panels <- "data/raw/gene_panel/"
genes_xenium <- read.csv(paste0(path_panels, "Xenium_hMulti_v1_metadata.csv")) %>% pull(Gene)
genes_cosmx <- readxl::read_xlsx(paste0(path_panels, "LBL-11178-05-Human-Universal-Cell-Characterization-Panel-Gene-Target-List.xlsx"),
                                 sheet = 2, skip = 1) %>% head(n=1000)%>% pull(`Display Name`)
intersect_genes <- intersect(genes_xenium, genes_cosmx)

tiledbURI <- "data/raw/CosMx/4fa4fec4-3ee0-4763-bf32-ed32a891af55_TileDB/"
tiledb_scdataset <- tiledbsc::SOMACollection$new(uri = tiledbURI,
                                                 verbose = T)

df <- data.frame()

for (sample in c("Bladder", "Breast", "Lung", "Colorectal", "Lymphoma", "Ovarian")) {
  print(sample)
  
  # COSMX
  cosmx <- readRDS(paste0(path_objects, "seurat_aligned/", sample, "_alignment_cosmx.rds"))
  gex <- sum(apply(cosmx@assays$RNA$counts, 2, sum))
  gex_int <- sum(apply(cosmx@assays$RNA$counts[intersect_genes, ], 2, sum))
  
  exp_neg <- tiledb_scdataset$somas$negprobes$X$members$counts$to_matrix(batch_mode = TRUE)[, colnames(cosmx)]
  negativex <- sum(apply(exp_neg, 2, sum))
  
  df <- rbind(df, data.frame(
    neg_rate1 = (negativex / (negativex + gex)) * 100,
    neg_rate2 = (negativex / (negativex + gex_int)) * 100,
    neg_rate3 = (negativex / (negativex + gex)) * (nrow(cosmx@assays$RNA@layers$counts) / nrow(exp_neg)) * 100,
    platform = "cosmx",
    sample = sample
  ))
  
  # XENIUM
  xenium <- readRDS(paste0(path_objects, "seurat_aligned/", sample, "_alignment_xenium.rds"))
  gex <- sum(apply(xenium@assays$Xenium@layers$counts, 2, sum))
  gex_int <- sum(apply(xenium@assays$Xenium$counts[intersect_genes, ], 2, sum))
  negativex <- sum(apply(xenium@assays$ControlProbe$counts, 2, sum))
  
  df <- rbind(df, data.frame(
    neg_rate1 = (negativex / (negativex + gex)) * 100,
    neg_rate2 = (negativex / (negativex + gex_int)) * 100,
    neg_rate3 = (negativex / (negativex + gex)) * (nrow(xenium@assays$Xenium@layers$counts) / nrow(xenium@assays$ControlProbe$counts)) * 100,
    platform = "xenium",
    sample = sample
  ))
}

write.table(df, "data/output/neg_rate/neg_rate_cosmx_xenium.txt")
