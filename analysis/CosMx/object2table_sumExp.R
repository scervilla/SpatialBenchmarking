library(Seurat)

path_seurat <- "data/objects/seurat_aligned/"
for (sample in c("Bladder","Colorectal", "Breast", "Lung", "Lymphoma", "Ovarian")) {
  object <- readRDS(paste0(path_seurat, sample, "_alignment_cosmx.rds"))  
  df <- data.frame(n=rowSums(object@assays$RNA$counts))
  write.table(df, paste0("data/output/df/", sample, "_cosmx_aligned_sumExp.txt"))
}



path_seurat <- "Desktop/IJC/projects/benchmark-scsp/seurat_aligned/"
write.table(df, paste0("Desktop/IJC/projects/benchmark-scsp/df/", sample, "_cosmx_aligned_sumExp.txt"))

write.table(df, paste0("Desktop/IJC/projects/benchmark-scsp/df/", sample, "_xenium_aligned_sumExp.txt"))

path_seurat <- "Desktop/IJC/pr"
cytassist <- readRDS("Desktop/IJC/projects/CytAssist/Colorectal_cyt/outs/RDS/Colorectal_cyt.rds")
visium <- readRDS("Desktop/IJC/projects/benchmark-scsp/Visium/Colorectal.rds")


write.table(df_exp_cyt, paste0("Desktop/IJC/projects/benchmark-scsp/df/", "df_", sample, "_cyt_sumExp.txt"))

write.table(df_exp_vis, paste0("Desktop/IJC/projects/benchmark-scsp/df/", "df_", sample, "_vis_sumExp.txt"))
