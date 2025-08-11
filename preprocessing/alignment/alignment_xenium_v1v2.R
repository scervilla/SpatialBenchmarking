library(VoltRon)
source("Desktop/IJC/scripts/utils/voltron_xenium_light.R")


sample <- "Ovarian"
sample_name <- sample

path_raw <- "data/raw/xenium_v1v2/"
path_object <- "data/objects/"

path_run_v1 <- "20250305__150014__COMPARATIVA_V1_vs_V2_MULTITISSUE_PANEL/"
path_run_v2 <- "20250224__143707__COMPARATIVA_V1_vs_V2_5K_PANEL/"

if (sample_name == "Colorectal") {sample_name <- "Colon"}

path_xenium_v1 <-  list.dirs(paste0(path_raw, path_run_v1), full.name=T, recursive=F)
path_xenium_v1 <- path_xenium_v1[grepl(stringr::str_to_upper(sample_name), path_xenium_v1)]


path_xenium_v2 <-  list.dirs(paste0(path_raw, path_run_v2), full.name=T, recursive=F)
path_xenium_v2 <- path_xenium_v2[grepl(stringr::str_to_upper(sample_name), path_xenium_v2)]

resolution = 7
Xen_V1 <- importXenium(path_xenium_v1, 
                       sample_name = "XeniumV1", import_molecules = F, resolution_level = resolution)
Xen_V2 <- importXenium(path_xenium_v2, 
                       sample_name = "XeniumV2", import_molecules = F, resolution_level = resolution)


xen_reg <- registerSpatialData(object_list = list(Xen_V1, Xen_V2))


###
scale_coord <- function(path, df) {
  image <-  image_read(paste0(path, "/morphology_lowres.tif"))
  imageinfo <- unlist(magick::image_info(image)[c("height")])
  range_coords <- c(0,imageinfo)
  df[,"y"] <- -( - range_coords[2] + df[,"y"] - range_coords[1])
  scale_factor <- 13.6
  df[,c("x","y")] <- df[,c("x","y")] * scale_factor
  return(df)
}

dfV1 <- xen_reg$registered_spat[[1]]@samples$XeniumV1@layer$Section1@assay$Xenium@image$image_1@coords %>% as.data.frame()
rownames(dfV1) <- stringr::str_replace(rownames(dfV1), "_Assay1", "")
dfV1 <- scale_coord(path_xenium_v1, dfV1)

dfV2 <- xen_reg$registered_spat[[2]]@samples$XeniumV2@layer$Section1@assay$Xenium@image$image_1_reg@coords %>% as.data.frame()
rownames(dfV2) <- stringr::str_replace(rownames(dfV2), "_Assay1", "")
dfV2 <- scale_coord(path_xenium_v1, dfV2)


seg <- xen_reg$registered_spat[[2]]@samples$XeniumV2@layer$Section1@assay$Xenium@image$image_1_reg@segments
seg <- do.call(rbind, seg)
seg <- scale_coord(path_xenium_v1, seg)
rownames(seg) <- NULL

path_results <- paste0(path_object, "XeniumV1V2_alignment/")
dfV1 <- cbind(rownames(dfV1), dfV1)
colnames(dfV1)[1] <- "cell_id"
write_parquet(dfV1, paste0(path_results, sample, "_V1_coords_aligned.parquet"))
dfV2 <- cbind(rownames(dfV2), dfV2)
colnames(dfV2)[1] <- "cell_id"
write_parquet(dfV2, paste0(path_results, sample, "_V2_coords_aligned.parquet"))
write_parquet(seg, paste0(path_results, sample, "_V2_seg_aligned.parquet"))

