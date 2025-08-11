# File: preprocessing/run_alignment.R

#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(dplyr)
library(data.table)
library(tiledbsc)

filter_cosmx <- function(cellCoords, sample) {
  filter_params <- list(
    Colorectal = list(tissue = "CxAms124 S2", expr = quote(y_slide_mm > 10 & x_slide_mm > 8)),
    Lung       = list(tissue = "CxAms124 S1", expr = quote(y_slide_mm > 6  & x_slide_mm > 7.5)),
    Breast     = list(tissue = "CxAms124 S1", expr = quote(y_slide_mm < 6.5)),
    Ovarian    = list(tissue = "CxAms124 S2", expr = quote(y_slide_mm < 8.5)),
    Lymphoma   = list(tissue = "CxAms124 S2", expr = quote(y_slide_mm > 7.5 & x_slide_mm < 8))
  )
  p <- filter_params[[sample]]
  cellCoords %>% filter(!!p$expr, Run_Tissue_name == p$tissue)
}

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]

alignment_params <- list(
  Bladder    = list(angle = -70,  shift_x = 2150, shift_y = 725),
  Breast     = list(angle = -165, shift_x = 2150, shift_y = 725),
  Colorectal = list(angle = 73,   shift_x = 400,  shift_y = 3050),
  Lung       = list(angle = 18,   shift_x = 720,  shift_y = 975),
  Lymphoma   = list(angle = -180, shift_x = 515,  shift_y = 930),
  Ovarian    = list(angle = -182, shift_x = 2400, shift_y = 470)
)

params <- alignment_params[[sample]]
angle_rad <- -params$angle * (pi / 180)

# Load CosMx Data
tiledbURI <- "data/raw/CosMx/4fa4fec4-3ee0-4763-bf32-ed32a891af55_TileDB/"
tiledb_scdataset <- tiledbsc::SOMACollection$new(uri = tiledbURI, verbose = TRUE)
cellCoords <- tiledb_scdataset$somas$RNA$obs$to_dataframe()
cellCoords <- filter_cosmx(cellCoords, sample)

# Process CosMx coordinates
df_cos <- cellCoords %>%
  mutate(
    x_slide_mm = x_slide_mm - min(x_slide_mm),
    y_slide_mm = y_slide_mm - min(y_slide_mm),
    x_slide_um = x_slide_mm * 1000,
    y_slide_um = y_slide_mm * 1000
  )
df_cos$filter <- "yes"

# Special handling for Bladder
if (sample == "Bladder") {
  df_cos_1 <- df_cos %>% filter(y_slide_um > 5000 & x_slide_um < 3500)
  block_1 <- rownames(df_cos_1)
  angle_radians_2 <- 2 * (pi / 180)
  df_cos_1 <- df_cos_1 %>% mutate(
    x_slide_um = x_slide_um * cos(angle_radians_2) - y_slide_um * sin(angle_radians_2),
    y_slide_um = x_slide_um * sin(angle_radians_2) + y_slide_um * cos(angle_radians_2)
  )
  df_cos[block_1,] <- df_cos_1
}

# Load Xenium Data
xenium_full_path <- paste0("results/xenium/", sample, "/", sample, "_seurat.rds")
xenium <- readRDS(xenium_full_path)
df_xen <- as.data.frame(xenium@images$fov@boundaries$centroids@coords)
colnames(df_xen) <- c("x_slide_um", "y_slide_um")
rownames(df_xen) <- colnames(xenium)

# Rotate Xenium coordinates
df_xen <- df_xen %>%
  mutate(
    x_slide_um = x_slide_um - min(x_slide_um),
    y_slide_um = y_slide_um - min(y_slide_um),
    x_rotated = x_slide_um * cos(angle_rad) - y_slide_um * sin(angle_rad),
    y_rotated = x_slide_um * sin(angle_rad) + y_slide_um * cos(angle_rad),
    x_rotated = x_rotated - min(x_rotated),
    y_rotated = y_rotated - min(y_rotated)
  )

# Additional filtering for Colorectal
if (sample == "Colorectal") {
  df_model <- df_xen[df_xen$y_slide_um < 10, ] %>% mutate(
    x_rotated = x_slide_um * cos(angle_rad) - y_slide_um * sin(angle_rad),
    y_rotated = x_slide_um * sin(angle_rad) + y_slide_um * cos(angle_rad)
  )
  df_model$x_rotated <- df_model$x_rotated - min(df_model$x_rotated)
  df_model$y_rotated <- df_model$y_rotated - min(df_model$y_rotated)
  model <- lm(y_rotated ~ x_rotated, data = df_model)
  slope <- coef(model)[2]
  intercept <- coef(model)[1]
  
  filter_out_cos <- sapply(1:nrow(df_cos), function(i) {
    df_cos[i, "y_slide_um_shift"] > slope * df_cos[i, "x_slide_um_shift"] + 7400
  })
  df_cos$filter[filter_out_cos] <- "no"
}

# Shift CosMx for alignment
if (sample == "Bladder") {
  df_cos$x_slide_um_shift <- df_cos$x_slide_um + 400
  df_cos$y_slide_um_shift <- df_cos$y_slide_um + 750
  df_cos_1 <- df_cos[block_1, ] %>% mutate(
    x_slide_um_shift = x_slide_um + 680,
    y_slide_um_shift = y_slide_um + 825
  )
  df_cos[block_1, ] <- df_cos_1
} else {
  x <- params$shift_x
  y <- params$shift_y
  df_cos <- df_cos %>%
    mutate(
      x_slide_um_shift = x_slide_um + x,
      y_slide_um_shift = y_slide_um + y
    )
}

# Assign FOV from CosMx to Xenium based on overlap
fovs <- unique(df_cos$fov)
df_xen$fov <- "out"
for (f in fovs) {
  f_coords <- df_cos %>% filter(fov == f)
  x_range <- range(f_coords$x_slide_um_shift)
  y_range <- range(f_coords$y_slide_um_shift)
  
  df_xen$fov[df_xen$x_rotated >= x_range[1] & df_xen$x_rotated <= x_range[2] &
               df_xen$y_rotated >= y_range[1] & df_xen$y_rotated <= y_range[2]] <- f
}

# Special FOV exclusion for Bladder
if (sample == "Bladder") {
  exclude <- c("51","52", "59", "66", "60")
  df_cos <- df_cos %>% filter(!fov %in% exclude)
  df_xen <- df_xen %>% filter(fov != "out" & !fov %in% exclude)
}

# Save raw alignment data
output_dir <- "data/outputs/"
dir.create(paste0(output_dir, "df"), recursive = TRUE, showWarnings = FALSE)
write.table(df_xen, paste0(output_dir, "df/", sample, "_xenium_aligned.txt"), row.names = FALSE)
write.table(df_cos, paste0(output_dir, "df/", sample, "_cosmx_aligned.txt"), row.names = FALSE)

# Filter coordinates
df_xen <- df_xen %>% filter(fov != "out")
df_cos <- df_cos %>% filter(filter != "no")

# Build CosMx Seurat object with rotated centroids
exp <- tiledb_scdataset$somas$RNA$X$members$counts$to_matrix(batch_mode = TRUE)[, rownames(df_cos)]
cosmx <- CreateSeuratObject(counts = as.matrix(exp))
cosmx <- AddMetaData(cosmx, df_cos)
coord <- cosmx@meta.data[, c("x_slide_um_shift", "y_slide_um_shift")]

coord_rotated <- coord %>%
  mutate(
    x_rotated = x_slide_um_shift * cos(-angle_rad) - y_slide_um_shift * sin(-angle_rad),
    y_rotated = x_slide_um_shift * sin(-angle_rad) + y_slide_um_shift * cos(-angle_rad)
  )

cents <- CreateCentroids(coord_rotated[, c("x_rotated", "y_rotated")])
cosmx[["fov"]] <- CreateFOV(coords = cents, type = "centroids")


# Save CosMx Seurat object
aligned_dir <- paste0(output_dir, "seurat_aligned")
dir.create(aligned_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(cosmx, paste0(aligned_dir, "/", sample, "_alignment_cosmx.rds"))

# Load and subset full Xenium Seurat object
xenium <- subset(xenium, cells = rownames(df_xen))
xenium@meta.data$fov <- df_xen$fov
saveRDS(xenium, paste0(aligned_dir, "/", sample, "_alignment_xenium.rds"))


