


library(Voyager)
library(stringr)
library(dplyr)
library(SpatialExperiment)
library(SpatialFeatureExperiment)


path_results <- "data/output/morani"
path_table <- "data/output/df/"

samples <- c("Colorectal", "Ovarian", "Lung", "Lymphoma", "Breast", "Bladder")


# Generate tiles of 50 um (rasterization)
df_tiles <- lapply(samples, function(sample) {
  
  df_info <- read.table(paste0(path_table, "df_", sample, ".txt"))
  df_xen <- read.table(paste0(path_table, "aligned_coord_", sample, "_xenium.txt"))
  df_cos <- read.table(paste0(path_table, "aligned_coord_", sample, "_cosmx.txt"))
  
  df_info[rownames(df_xen),c("x", "y")] <- df_xen[,c("x_shared", "y_shared")]
  df_info[rownames(df_cos),c("x", "y")] <- df_cos[,c("x_shared", "y_shared")]
  
  df_info <- df_info %>% filter(!is.na(ct))
  
  
  tile_size <- 50
  df_tiles <- df_info %>% mutate(x_centroid_10 = floor(x / tile_size ),
                                 y_centroid_10 = floor(y / tile_size )) %>% 
    mutate(fov = paste0(x_centroid_10, "_",y_centroid_10),
           x_tile = x_centroid_10 * tile_size  + tile_size/2,
           y_tile = y_centroid_10 * tile_size  + tile_size/2)
  
  df_tiles
})
df_tiles <- do.call(rbind, df_tiles)



df_sum <- df_tiles %>% 
  group_by(sample, platform, ct,fov) %>%
  summarise(n = n(), .groups = "drop") %>% 
  group_by(sample, platform, fov) %>%
  mutate(prop = n / sum(n)) %>%
  tidyr::pivot_wider(names_from = ct, values_from = c(n, prop), values_fill = list(n = 0, prop = 0))


df_prop <- left_join(df_sum, (df_tiles %>% select(sample, platform, fov, x_centroid_10, y_centroid_10) %>% unique()),
                     by = c("sample", "platform", "fov"))


# Compute Moran's I from tiles
combinations <- expand.grid(sample = samples, platform = c("cosmx", "xenium"))
results <- apply(combinations, 1, function(row) {
  tile_sub <- df_prop %>% filter(sample == row["sample"] & platform == row["platform"]) %>% ungroup() 
  
  sce = SingleCellExperiment(assays =list(SCT = tile_sub %>% select(n_B_cell:`n_T_cell:effector`) %>% t()))
  coord <- tile_sub %>% select(x=x_centroid_10,y=y_centroid_10)
  colData(sce) <- cbind(colData(sce), coord)
  spe <- toSpatialExperiment(sce, spatialCoordsNames = c("x", "y")) 
  sfe <- toSpatialFeatureExperiment(spe)
  
  
  colGraph(sfe, "knn5") <- findSpatialNeighbors(sfe, method = "knearneigh", 
                                                dist_type = "idw", k = 8, 
                                                style = "W")  
  
  sfe <- runMoransI(sfe, colGraphName = "knn5", exprs_values = "SCT")
  df <- as.data.frame(rowData(sfe))
  df$ct <- rownames(df)
  df$sample <- row["sample"] 
  df$platform <-  row["platform"]
  df
})

names(results) <- apply(combinations, 1, function(x) {
  paste0(x, collapse = "_")
})
results <- do.call(rbind, results)

write.table(results, paste0(path_results, "morani_ct.txt"))



