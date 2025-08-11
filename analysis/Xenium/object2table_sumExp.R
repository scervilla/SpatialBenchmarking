library(Seurat)

path_seurat <- "data/objects/seurat_aligned/"
for (sample in c("Bladder","Colorectal", "Breast", "Lung", "Lymphoma", "Ovarian")) {
  object <- readRDS(paste0(path_seurat, sample, "_alignment_xenium.rds"))  
  df <- data.frame(n=rowSums(object@assays$Xenium$counts))
  write.table(df, paste0("data/output/df/", sample, "_xenium_aligned_sumExp.txt"))
}