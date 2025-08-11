library(RANN)
library(sf)
library(concaveman)
library(arrow)


path_object <- "data/objects/"
path_results <- paste0(path_object, "XeniumV1V2_alignment/")

for (sample in c("Lung", "Lymphoma", "Colorectal")) {
  print(sample)
  
  dfV1 <- read_parquet(paste0(path_results, sample, "_V1_coords_aligned.parquet"))
  dfV1 <- dfV1 %>% tibble::column_to_rownames("cell_id")
  dfV2 <- read_parquet(paste0(path_results, sample, "_V2_coords_aligned.parquet"))
  dfV2 <- dfV2 %>% tibble::column_to_rownames("cell_id")
  
  dfV2 <- dfV2[,1:2]
  dfV1 <- dfV1[,1:2]
  
  min_distance <- 20
  nn <- nn2(data = dfV1, query = dfV2, k = 2)
  min_dists <- nn$nn.dists[,2]
  validV2 <- min_dists <= min_distance
  nn <- nn2(data = dfV2, query = dfV1, k = 2)
  min_dists <- nn$nn.dists[,2]
  validV1 <- min_dists <= min_distance
  
  
  
  dfV2$group <- validV2
  dfV1$group <- validV1

  
  
  if (sum(!validV1)  <= sum(!validV2)) {
    df_blue <- dfV2 %>% filter(group) %>% dplyr::select(-group)
  }  else {
    df_blue <- dfV1 %>% filter(group) %>% dplyr::select(-group)
  }
  
  points_sf <- st_as_sf(df_blue, coords = c("x", "y"), crs = 4326)
  boundary <- concaveman(points_sf)
  boundary_coords <- st_coordinates(boundary)[, 1:2]
  #boundary plot
  ggplot(boundary_coords, aes(X,Y)) + geom_point(size=0.1) + coord_equal()
  
  dfV1$inside  <- point.in.polygon(dfV1$x, dfV1$y, boundary_coords[,1], boundary_coords[,2])
  dfV2$inside  <- point.in.polygon(dfV2$x, dfV2$y, boundary_coords[,1], boundary_coords[,2])
  
  

  cells_V1 <- dfV1 %>% filter(inside == 1) %>% rownames()
  cells_V2 <- dfV2 %>% filter(inside == 1) %>% rownames()
  
  write.table(cells_V1, paste0(path_results, "/cells_filter_", sample, "_V1.txt"))
  write.table(cells_V2, paste0(path_results, "/cells_filter_", sample, "_V2.txt"))
}