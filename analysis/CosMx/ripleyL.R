library(spatstat.geom)
library(spatstat.explore)
library(dplyr)
library(Seurat)

args = commandArgs(trailingOnly=TRUE)
sample <- args[1]


path_object <- "data/objects/seurat_aligned/"

platform <- "cosmx"
object <-  readRDS(paste0(path_object, sample, "_alignment_", platform, ".rds"))
ct <- readRDS(paste0("data/output/SingleR/predictions_", sample,"_", platform," _final.rds"))

coords <- GetTissueCoordinates(object[["fov"]][["centroids"]])
coords$ct <-  ct$pruned.labels 

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

rkc <- do.call(rbind, lapply(levels(ppo$marks), function(x) {
  test <- Lcross(ppo, x, x, r = radii, correction = "isotropic")
  df_test <- data.frame(reference = "",
                        neighbor = x,
                        radius = test$r)
  
  df_test$score <- test$iso - test$theo
  
  auc <- sum(diff(df_test$radius) * 
               (head(df_test$score, -1) + tail(df_test$score, -1)) / 2)
  
  df_test$AUC <- auc
  
  return(df_test)
}))

write.table(rkc, paste0("data/output/Ripley/ripley_L_", sample, "_xenium.txt"))
