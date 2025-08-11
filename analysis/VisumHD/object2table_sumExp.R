library(Seurat)

path_seurat <- "data/objects/seurat_aligned/"
for (sample in c("Lymphoma", "Lung", "Colorectal", "Bladder")) {
  object <- readRDS(paste0(path_seurat, sample, "_HD.rds"))  
  df <- as.data.frame(rowSums(object@assays$Spatial.008um$counts))
  write.table(df, paste0("data/output/df/", sample, "_HD_sumExp.txt"))
}