library(Seurat)
library(dplyr)


path_seurat <- "data/objects/seurat/"
path_annotation <- "data/output/SingleR/"
path_outputs <- "data/output/df/"

for (sample in c( "Lymphoma", "Breast", "Bladder", "Colorectal", "Lung")) {
  print(sample)
  visium <- readRDS(paste0(path_seurat, sample, "_vis.rds"))
  cytassist <- readRDS(paste0(path_seurat, sample, "_cyt.rds"))
  
  
  DefaultAssay(cytassist) <- "RNA"
  DefaultAssay(visium) <- "RNA"
  

  visium$n_genes <- apply(visium@assays$RNA@counts, 2, function(x) {
    sum(x > 0)
  })
  
  visium$n_umi <- apply(visium@assays$RNA@counts, 2, function(x) {
    sum(x)
  })
  
  cytassist$n_genes <- apply(cytassist@assays$RNA@layers$counts, 2, function(x) {
    sum(x > 0)
  })
  
  cytassist$n_umi <- apply(cytassist@assays$RNA@layers$counts, 2, function(x) {
    sum(x)
  })
  
  
  cytassist$platform <- "Cytassist"
  visium$platform <- "Visium"
  
  df <- rbind(visium@meta.data[,c("n_genes", "n_umi", "platform")],
              cytassist@meta.data[,c("n_genes", "n_umi", "platform")])
  df$sample <- sample
  write.table(paste0(df, paste0(path_outputs, "df_", sample, "_cyt_vis.txt")))

  
  intersect_genes <- intersect(rownames(visium),
                               rownames(cytassist))

 # df_exp <- data.frame(row.names = intersect_genes,
 #                      visium=apply(visium@assays$RNA$counts[intersect_genes, ], 1, sum),
#                       cytassist=apply(cytassist@assays$RNA$counts[intersect_genes, ], 1, sum))
  
 # write.table(df_exp, paste0(path_outputs, "df_", sample, "_cyt_vis_sumExp.txt"))
  
  df_exp_cyt <- data.frame(n=rowSums(cytassist@assays$RNA$counts))
  write.table(df_exp_cyt, paste0(path_outputs, "df_", sample, "_cyt_sumExp.txt"))
  
  df_exp_vis <- data.frame(n=rowSums(visium@assays$RNA$counts))
  write.table(df_exp_cyt, paste0(path_outputs, "df_", sample, "_vis_sumExp.txt"))



  df_genes <-  data.frame(
    category = c("Both", "Only in Visium", "Only in Cytassist"),
    count = c(
      length(intersect_genes),
      sum(!rownames(visium) %in% intersect_genes),
      sum(!rownames(cytassist) %in% intersect_genes)
    ), sample = sample)
  write.table(df_genes, paste0(path_outputs, "genes_", sample, "_cyt_vis.txt"))
}





