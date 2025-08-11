library(Seurat)
library(dplyr)

path_objects <- "data/objects/seurat/"
path_output <- "data/output/"

sample <- "Ovarian"
object <- readRDS(paste0(path_objects, sample, "_HD.rds"))
ct <- readRDS(paste0(path_output, "SingleR/predictions_", sample, "_VisiumHD_final.rds"))

object$ct <- ct$pruned.labels
Idents(object) <- "ct"
object <- subset(object, idents="Macrophage")


object <- SCTransform(object, assay = "Spatial.008um")
object <- RunPCA(object)
object <- RunUMAP(object, dims = 1:20, reduction = "pca", reduction.name = "umap")
DimPlot(object)


object <- FindNeighbors(object, reduction = "pca", dims = 1:20)
object <- FindClusters(object, resolution = 0.3)
write.table(object@meta.data[,c("ct","seurat_clusters")], 
            paste0(path_output, "subtype/", sample, "_subtype_HD.txt"))

markers <- FindAllMarkers(object, only.pos = T, min.pct = 0.1, logfc.threshold = 1)
write.table(markers, paste0(path_output, "subtype/", sample, "_subtype_HD_markers.txt"))

