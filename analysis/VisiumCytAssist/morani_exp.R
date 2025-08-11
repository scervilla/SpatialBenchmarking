library(Seurat)
library(semla)
library(dplyr)
library(ggplot2)

samples <- c("Breast", "Bladder", "Colorectal", "Lung")
object_path <- "data/objects/"

for (sample in samples) {
  message(sample)
  
  visium_path <- paste0(object_path, "/seurat/", sample, "_vis.rds")
  cytassist_path <- paste0(object_path, "/seurat/", sample, "_vis.rds")
  
  visium <- readRDS(visium_path)
  cytassist <- readRDS(cytassist_path)
  
  genes_keep <- rownames(cytassist@assays$RNA)[
    Matrix::rowSums(cytassist@assays$RNA@layers$counts) >= 5
  ]
  
  DefaultAssay(cytassist) <- "RNA"
  DefaultAssay(visium) <- "RNA"
  
  cytassist <- subset(cytassist, features = genes_keep)
  
  visium <- visium %>% UpdateSeuratForSemla() %>% SCTransform(assay = "RNA")
  cytassist <- cytassist %>% UpdateSeuratForSemla() %>% SCTransform(assay = "RNA")
  
  spatgenes <- CorSpatialFeatures(visium, features = rownames(visium@assays$SCT), assay_use = "SCT") %>% as.data.frame()
  rownames(spatgenes) <- spatgenes$gene
  
  spatgenes_cyt <- CorSpatialFeatures(cytassist, features = rownames(cytassist), assay_use = "SCT") %>% as.data.frame()
  rownames(spatgenes_cyt) <- spatgenes_cyt$gene
  
  intersect_genes <- intersect(rownames(spatgenes_cyt), rownames(spatgenes))
  
  df_autocor <- data.frame(
    row.names = intersect_genes,
    visium = spatgenes[intersect_genes, "cor"],
    cytassist = spatgenes_cyt[intersect_genes, "cor"]
  )
  
  write.table(df_autocor,paste0("data/output/morani/morani_", sample, "_cyt_vis.txt"))
} 