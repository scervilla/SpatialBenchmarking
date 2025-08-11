library(Seurat)
library(dplyr)
library(spatstat.geom)
library(spatstat.explore)

path_objects <- "data/objects/seurat_aligned/"
path_output <- "data/output/"

sample <- "Ovarian"
object <- readRDS(paste0(path_objects, sample, "_xenium_in_hd.rds"))
ct <- readRDS(paste0(path_output, "SingleR/predictions_", sample, "_xenium_xen_hd_hpca.rds"))
ct_macro <- read.table(paste0(path_output, "subtype/", sample, "_subtype_xenium_markers.txt"))

coords <- GetTissueCoordinates(object@images[[1]][["centroids"]])
coords$ct <-  ct$pruned.labels
rownames(coords) <- coords$cell
coords[rownames(ct_macro), "ct"] <- ct_macro$seurat_clusters


hpca.se <- celldex::HumanPrimaryCellAtlasData()
combined_t <- colData(hpca.se) %>% as.data.frame() %>% filter(label.main == "T_cells") %>% pull(label.fine) %>% unique
coords$ct[coords$ct %in% combined_t] <- "T_cell"

coords <- coords[!is.na(coords$ct),] 

win <- owin(
  xrange = range(coords$x),
  yrange = range(coords$y)
)

ppo <- ppp(
  x = coords$x,
  y = coords$y,
  window = win,
  marks = factor(coords$ct)
)


radii <- seq(0, 1000, by = 5)

rkc_cancer <- do.call(rbind, lapply(0:11, function(x) {
  test <- Lcross(ppo, "Epithelial_cells", x, r = radii, correction = "isotropic")
  df_test <- data.frame(reference = "Epithelial_cells",
                        neighbor = x,
                        radius = test$r)
  ## iso-corrected minus theoretical
  df_test$score <- test$iso - test$theo
  
  # Compute AUC using trapezoidal rule
  auc <- sum(diff(df_test$radius) * 
               (head(df_test$score, -1) + tail(df_test$score, -1)) / 2)
  
  # Add AUC as a constant column (repeated for all rows of this neighbor type)
  df_test$AUC <- auc
  
  return(df_test)
}))
write.table(rkc_cancer, paste0(path_output,"Ripley/ripley_L_", sample, "_xenium_cancer.txt"))


rkc_tcells <- do.call(rbind, lapply(0:11, function(x) {
  test <- Lcross(ppo, "T_cell", x, r = radii, correction = "isotropic")
  df_test <- data.frame(reference = "T_cell",
                        neighbor = x,
                        radius = test$r)
  ## iso-corrected minus theoretical
  df_test$score <- test$iso - test$theo
  
  # Compute AUC using trapezoidal rule
  auc <- sum(diff(df_test$radius) * 
               (head(df_test$score, -1) + tail(df_test$score, -1)) / 2)
  
  # Add AUC as a constant column (repeated for all rows of this neighbor type)
  df_test$AUC <- auc
  
  return(df_test)
}))
write.table(rkc_cancer, paste0(path_output, "Ripley/ripley_L_", sample, "_xenium_tcell.txt"))