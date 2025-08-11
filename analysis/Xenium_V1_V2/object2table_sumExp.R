library(Seurat)


path_seurat <- "data/objects/seurat_aligned/"
for (sample in c("Colorectal", "Breast", "Lung", "Lymphoma", "Ovarian")) {
  for (version in c("V1", "V2")) {
    object <- readRDS(paste0(path_seurat, sample, "_", version, "_aligned.rds"))  
    df <- as.data.frame(rowSums(object@assays$Xenium$counts))
    write.table(df, paste0("data/output/df/", sample, "_", version, "_aligned_sumExp.txt"))
  }
}