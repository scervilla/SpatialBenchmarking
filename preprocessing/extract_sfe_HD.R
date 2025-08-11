library(Seurat)


args = commandArgs(trailingOnly=TRUE)
sample <- args[1]

file_path <- "data/object/sfe/"
path_objects <- "data/objects/seurat/"

object <- readRDS(paste0(path_objects, sample, "_HD_sct.rds"))

saveRDS(object@assays$SCT$data, paste0(file_path,  sample,  "_HD_sct_exp.rds"))
write.table(object@images[[1]]@boundaries$centroids@coords,paste0(file_path,  sample, "_HD_sct_coord.rds"))
