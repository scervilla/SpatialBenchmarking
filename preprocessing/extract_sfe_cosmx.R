library(Seurat)


args = commandArgs(trailingOnly=TRUE)
sample <- args[1]
platform <- "cosmx"

file_path <- "data/object/sfe/"
path_objects <- "data/objects/seurat_aligned/"

object <- readRDS(paste0(path_objects,  sample, "_alignment_", platform, ".rds"))

saveRDS(object@assays$SCT$data, paste0(file_path, sample, "_alignment_", platform, "_sct_exp.rds"))
write.table(object@images[[1]]@boundaries$centroids@coords,paste0(file_path,   sample, "_alignment_", platform,  "_sct_coord.rds"))
