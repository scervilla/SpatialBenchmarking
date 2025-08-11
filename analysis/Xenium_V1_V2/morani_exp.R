library(BiocParallel)
library(SpatialExperiment)
library(SpatialFeatureExperiment)
library(Voyager)

args = commandArgs(trailingOnly=TRUE)
sample <- args[1]

file_path <- "data/object/sfe/"
output_path <- "data/output/"

object <- readRDS(paste0(file_path,  sample, "_aligned_sct_exp.rds"))
coord <- read.table(paste0(file_path, sample, "_aligned_sct_coord.rds"))

sce = SingleCellExperiment(assays =list(SCT = object))
colData(sce) <- cbind(colData(sce), coord)
spe <- toSpatialExperiment(sce, spatialCoordsNames = c("x", "y")) 
sfe <- toSpatialFeatureExperiment(spe)

colGraph(sfe, "knn5") <- findSpatialNeighbors(sfe, method = "knearneigh", 
                                              dist_type = "idw", k = 5, 
                                              style = "W")  

sfe <- runMoransI(sfe, colGraphName = "knn5", BPPARAM = MulticoreParam(12), exprs_values = "SCT")
write.table(as.data.frame(rowData(sfe)), paste0(output_path, "morani/", sample,"_morani_aligned.txt"))
