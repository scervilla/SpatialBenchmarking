library(Seurat)
library(dplyr)

path_objects <- "data/objects/seurat_aligned/"
path_output <- "data/output/"

sample <- "Ovarian"
object <- readRDS(paste0(path_objects, sample, "_xenium_in_hd.rds"))
ct <- readRDS(paste0(path_output, "SingleR/predictions_", sample, "_xenium_xen_hd_hpca.rds"))

object$ct <- ct$pruned.labels
Idents(object) <- "ct"
object <- subset(object, idents="Macrophage")


xenium <- SCTransform(xenium, assay = "Xenium")
object <- RunPCA(object)
object <- RunUMAP(object, dims = 1:20, reduction = "pca", reduction.name = "umap")
DimPlot(object)


object <- FindNeighbors(object, reduction = "pca", dims = 1:20)
object <- FindClusters(object, resolution = 0.5)
write.table(object@meta.data[,c("ct","seurat_clusters")], 
            paste0(path_output, "subtype/", sample, "_subtype_xenium.txt"))

markers <- FindAllMarkers(object, only.pos = T, min.pct = 0.1, logfc.threshold = 1)
write.table(markers, paste0(path_output, "subtype/", sample, "_subtype_xenium_markers.txt"))

