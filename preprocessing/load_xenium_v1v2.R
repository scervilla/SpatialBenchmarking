library(progressr)
library(Seurat)
library(arrow)
library(stringr)

## Skip transcript location for the Seurat object
ReadXenium <- function (data.dir, outs = c("matrix"), type = "centroids", 
                        mols.qv.threshold = 20) {
  type <- match.arg(arg = type, choices = c("centroids", "segmentations"), 
                    several.ok = TRUE)
  outs <- match.arg(arg = outs, choices = c("matrix"), 
                    several.ok = TRUE)
  outs <- c(outs, type)
  has_dt <- requireNamespace("data.table", quietly = TRUE) && 
    requireNamespace("R.utils", quietly = TRUE)
  data <- sapply(outs, function(otype) {
    switch(EXPR = otype, matrix = {
      pmtx <- progressor()
      pmtx(message = "Reading counts matrix", class = "sticky", 
           amount = 0)
      matrix <- suppressWarnings(Read10X(data.dir = file.path(data.dir, 
                                                              "cell_feature_matrix/")))
      pmtx(type = "finish")
      matrix
    }, centroids = {
      pcents <- progressor()
      pcents(message = "Loading cell centroids", class = "sticky", 
             amount = 0)
      if (has_dt) {
        cell_info <- as.data.frame(data.table::fread(file.path(data.dir, 
                                                               "cells.csv.gz")))
      } else {
        cell_info <- read.csv(file.path(data.dir, "cells.csv.gz"))
      }
      cell_centroid_df <- data.frame(x = cell_info$x_centroid, 
                                     y = cell_info$y_centroid, cell = cell_info$cell_id, 
                                     stringsAsFactors = FALSE)
      pcents(type = "finish")
      cell_centroid_df
    }, segmentations = {
      psegs <- progressor()
      psegs(message = "Loading cell segmentations", class = "sticky", 
            amount = 0)
      if (has_dt) {
        cell_boundaries_df <- as.data.frame(data.table::fread(file.path(data.dir, 
                                                                        "cell_boundaries.csv.gz")))
      } else {
        cell_boundaries_df <- read.csv(file.path(data.dir, 
                                                 "cell_boundaries.csv.gz"), stringsAsFactors = FALSE)
      }
      names(cell_boundaries_df) <- c("cell", "x", "y")
      psegs(type = "finish")
      cell_boundaries_df
    }, stop("Unknown Xenium input type: ", otype))
  }, USE.NAMES = TRUE)
  return(data)
}

LoadXenium <- function (data.dir, fov = "fov", assay = "Xenium") 
{
  
  data <- ReadXenium(data.dir = data.dir, type = c("centroids", 
                                                   "segmentations"), )
  segmentations.data <- list(centroids = CreateCentroids(data$centroids), 
                             segmentation = CreateSegmentation(data$segmentations))
  coords <- CreateFOV(coords = segmentations.data, type = c("segmentation", 
                                                            "centroids"), assay = assay)
  xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], 
                                   assay = assay)
  if ("Blank Codeword" %in% names(data$matrix)) 
    xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Blank Codeword"]])
  else xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
  xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
  xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
  xenium.obj[[fov]] <- coords
  return(xenium.obj)
}



args <- commandArgs(trailingOnly = TRUE)
sample <- args [1]

sample_name <- stringr::str_split(sample, pattern = "_", simplify = T)[1,1]
version <- stringr::str_split(sample, pattern = "_", simplify = T)[1,2]

path_raw <- "data/raw/xenium_v1v2/"

path_run <- ifelse(version == "V1", 
                   "20250305__150014__COMPARATIVA_V1_vs_V2_MULTITISSUE_PANEL/", 
                   "20250224__143707__COMPARATIVA_V1_vs_V2_5K_PANEL/")

if (sample_name == "Colorectal") {sample_name <- "Colon"}

path_xenium <-  list.dirs(paste0(path_raw, path_run), full.name=T, recursive=F)
path_xenium <- path_xenium[grepl(stringr::str_to_upper(sample_name), path_xenium)]

#load xenium
xenium.obj <- LoadXenium(path_xenium, fov = "fov")
cell_info <- read_parquet(paste0(path_xenium, "/cells.parquet"))
xenium.obj <- AddMetaData(xenium.obj, cell_info)

# remove cells with 0 counts
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)

#QC
n_panel <- 377
#normalize the number of transcripts and genes by the size of the panel used
xenium.obj$nCount_Xenium_norm <- xenium.obj$nCount_Xenium/n_panel
xenium.obj$nFeature_Xenium_norm <- xenium.obj$nFeature_Xenium/n_panel

#compute the proportion of the nucleus within a cell 
xenium.obj$prop_nuc <- xenium.obj$nucleus_area / xenium.obj$cell_area

saveRDS(xenium.obj, paste0("~/projects/benchmark-scsp/seurat_aligned/", sample, ".rds"))
