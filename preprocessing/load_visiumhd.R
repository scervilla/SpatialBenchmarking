library(Seurat)
library(patchwork)
library(dplyr)

thresholds <- c(Ovarian = 10,
                Lymphoma = 30,
                Lung = 3,
                Colorectal = 15,
                Bladder = 5)

main_path <- "data/raw/VisiumHD/"
for (sample in  c( "Lung", "Colorectal", "Bladder", "Ovarian")) {
  print(sample)
  localdir <- paste0(main_path, sample, "/outs/")
  object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8))
  #mitochondrial percentage
  object <- PercentageFeatureSet(object, "^MT-", col.name = "percent_mito")
  cells_keep <- colnames(object)[object$nCount_Spatial.008um > thresholds[sample]]
  object <- subset(object, cells=cells_keep)
  
  write.table(object@meta.data, paste0(main_path, "/df/df_", sample, "_HD.txt"))
  saveRDS(object, paste0(main_path, "/RDS/", sample, "_HD.rds"))
  
}
