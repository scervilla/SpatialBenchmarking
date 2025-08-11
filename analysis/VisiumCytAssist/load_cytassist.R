library(semla)
library(Seurat)
library(dplyr)

path_folder <- "data/raw/CytAssist/"
for (sample in c("Lung", "Colorectal", "Lymphoma","Ovarian")) {
  print(sample)
  mat <- Seurat::Read10X_h5(paste0(path_folder, sample,"_cyt/outs/filtered_feature_bc_matrix.h5"))
  #Add Gene Expression
  seurat <- CreateSeuratObject(counts = mat$`Gene Expression`)
  img <- Seurat::Read10X_Image(image.dir = paste0(path_folder, sample,"_cyt/outs/spatial/"), image.name = "tissue_lowres_image.png")
  Seurat::DefaultAssay(object = img) <- 'RNA'
  img <- img[colnames(x = seurat)]
  seurat[[sample]] <- img
  seurat <- semla::UpdateSeuratForSemla(seurat)
  
  #Add proteins
  seurat[["Protein"]] <- CreateAssayObject(counts = mat$`Antibody Capture`)
  if (file.exists(paste0(path_folder, "filter_spots/", sample, ".csv"))) {
    meta <- read.csv(paste0(path_folder, "filter_spots/", sample, ".csv"))
    seurat <- subset(seurat, cells = setdiff(colnames(seurat), meta$Barcode)  )
  }
  
  if (!file.exists(paste0(path_folder, sample,"_cyt/outs/Figures"))) dir.create(paste0(path_folder, sample,"_cyt/outs/Figures"))
  SpatialFeaturePlot(seurat, features = c("nCount_RNA", "nFeature_RNA", "nCount_Protein", "nFeature_Protein"), 
                     pt.size.factor = 2, ncol = 2)
  ggsave(paste0(path_folder, sample,"_cyt/outs/Figures/QC_", sample, ".png"), width = 10, height = 10)
  
  seurat <- SCTransform(seurat, assay = "RNA")
  if (!file.exists(paste0(path_folder, sample,"_cyt/outs/RDS"))) dir.create(paste0(path_folder, sample,"_cyt/outs/RDS"))
  saveRDS(seurat, paste0(path_folder, sample,"_cyt/outs/RDS/", sample, "_cyt.rds"))
}
