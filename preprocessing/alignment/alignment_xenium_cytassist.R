library(dplyr)
library(Seurat)
library(VoltRon)

path_raw <- "data/raw/"
path_object <- "data/objects/"
path_output <- "data/output/"


xenium_folder <- list("Bladder" = "output-XETG00051__0003605__S150192__20230922__150140/",
                      "Lung" = "output-XETG00051__0003605__S150193__20230922__150140/",
                      "Breast" = "output-XETG00051__0003605__S150194__20230922__150140/",
                      "Colorectal" = "output-XETG00051__0003606__S150195__20230922__150140/",
                      "Lymphoma" = "output-XETG00051__0003606__S150196__20230922__150140/",
                      "Ovarian" = "output-XETG00051__0003606__S150197__20230922__150140/"
)


sample <- ""

## Align samples 
Xen_R1 <- importXenium(path_raw, "Xenium/", xenium_folder[[sample]], sample_name = "XeniumR1")
Vis <- importVisium(paste0(path_raw, "CytAssist/", sample, "_cyt/outs/"), sample_name = "VisiumR1")

keypoints <- readRDS(paste0(path_object, "/CytAssist_alignment/landmarks/", sample, "_landmarks.rds"))
xen_reg <- registerSpatialData(object_list = list(Xen_R1, Vis), keypoints = keypoints)

rm(Xen_R1, Vis)



df1 <- xen_reg$registered_spat[[2]]@samples$VisiumR1@layer$Section1@assay$Visium@coords_reg %>% as.data.frame()
df2 <- xen_reg$registered_spat[[1]]@samples$XeniumR1@layer$Section1@assay$Xenium@coords %>% as.data.frame()
rm(xen_reg)

#Rotate and assign FOV
xenium <- readRDS(paste0(path_object, "seurat/", sample, "_seurat.rds"))
rownames(df2) <- stringr::str_replace(rownames(df2), "_Assay1", "")

## Assign fov to Xenium
cytassist <- readRDS(paste0(path_object, "seurat/", sample, "_cyt.rds"))
rownames(df1) <- stringr::str_replace(rownames(df1), "_Assay1", "")
rm(xen_reg)
df1 <- df1 %>% subset(x > min(df2$x) & x < max(df2$x) &
                        y > min(df2$y) & y < max(df2$y))

image_file <- paste0(path_raw, "Xenium/", xenium_folder[[sample]], "/morphology_lowres.tif")
image <- magick::image_read(image_file)
imageinfo <- unlist(magick::image_info(image)[c("height")])
range_coords <- c(0, imageinfo)

resolution_level <- 7
scaleparam <- 0.2125 * (2^(resolution_level - 1))

df1[, 2] <- range_coords[2] - df1[, 2] + range_coords[1]
df1 <- df1*scaleparam



euclidean_distance <- function(x1, y1, x2, y2) {
  # Compute the squared differences
  dx <- x2 - x1
  dy <- y2 - y1
  
  # Compute the squared Euclidean distance
  distance_squared <- dx^2 + dy^2
  
  # Return the square root to get the Euclidean distance
  distance <- sqrt(distance_squared)
  
  return(distance)
}

df_xen <- xenium@meta.data
df_xen$within <- "no"
for (cell in 1:nrow(df1)) {
  center <- center <- df1[cell,]
  tmp <- df_xen %>% filter(x_centroid >= (center$x-55) & x_centroid <= (center$x+55) &
                             y_centroid >= (center$y-55) & y_centroid <= (center$y+55) ) %>% 
    select(x_centroid, y_centroid)
  
  
  if (nrow(tmp) > 0) {
    dists <- euclidean_distance(center$x, center$y, tmp$x_centroid, tmp$y_centroid)
    keep_cells <- rownames(tmp)[dists <= 55/2]
    df_xen[keep_cells, "within"] <- rownames(df1)[cell]
  }
}
xenium@meta.data <- df_xen


ct <-  readRDS(paste0("data/output/SingleR/predictions_", sample, "_", platform, "_whole_final.rds"))
df_xen$cell_type <- ct[rownames(df_xen), "pruned.labels"]
df_xen <- df_xen[!is.na(df_xen$cell_type)]
df_xen$cell_type <- factor(df_xen$cell_type)

cell_type <- c()
for (spot in unique(df_xen$within)) {
  cell_type <- rbind(cell_type, df_xen %>% filter(within == spot) %>% pull(cell_type) %>% table)
}
cell_type <- cell_type[-1,]
rownames(cell_type) <- unique(df_xen$within)[-1]
cell_type <- round(cell_type / rowSums(cell_type), 3)
colnames(cell_type) <- colnames(cell_type)

cytassist <- subset(cytassist, cells=rownames(cell_type))
cytassist <- AddMetaData(cytassist, metadata = cell_type)

# Sum expression
exp_spot <- AggregateExpression(xenium,  group.by = "within")[[1]] %>% as.matrix()
exp_spot <- exp_spot[,colnames(cytassist)]
rownames(exp_spot) <- paste0(rownames(exp_spot), "_pseudo")
cytassist[["pseudo"]] <- CreateAssayObject(exp_spot)

saveRDS(cytassist, paste0(path_object, "seurat_aligned/",sample, "_cyt_xenium.rds"))
