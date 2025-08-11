library(Seurat)
library(ggplot2)
library(stringr)
library(data.table)
args = commandArgs(trailingOnly=TRUE)

sample <- args[1]
path_folder <- "data/raw/CytAssist/"
raw_data = list("Bladder" = "output-XETG00051__0003605__S150192__20230922__150140/",
                "Lung" = "output-XETG00051__0003605__S150193__20230922__150140/",
                "Breast" = "output-XETG00051__0003605__S150194__20230922__150140/",
                "Colorectal" = "output-XETG00051__0003606__S150195__20230922__150140/",
                "Lymphoma"= "output-XETG00051__0003606__S150196__20230922__150140/",
                "Ovarian" = "output-XETG00051__0003606__S150197__20230922__150140/")


# Load the Xenium data
xenium.obj <- LoadXenium(paste0(path_folder, raw_data[sample]), fov = "fov")
cell_info <- fread(paste0(path_folder, raw_data[sample],"cells.csv.gz"))
xenium.obj <- AddMetaData(xenium.obj, cell_info)

#QC
n_panel <- 377
xenium.obj$nCount_Xenium_norm <- xenium.obj$nCount_Xenium/n_panel
xenium.obj$nFeature_Xenium_norm <- xenium.obj$nFeature_Xenium/n_panel

#compute the proportion of the nucleus within a cell 
xenium.obj$prop_nuc <- xenium.obj$nucleus_area / xenium.obj$cell_area


xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0 & nFeature_Xenium > 0)

xenium.obj <- SCTransform(xenium.obj, assay = "Xenium")
saveRDS(xenium.obj, paste0(path_folder, "Results/", sample, "/", sample, "_seurat.rds"))
