library(dplyr)
library(Seurat)

path_objects <- "data/objects/"

samples <- c("Breast", "Lung", "Colorectal", "Lymphoma", "Ovarian")
df_neg = data.frame()
for (sample in samples) {
  print(sample)
  for (platform in c("V1", "V2")) {
    object <- readRDS(paste0(path_objects, "/seurat_aligned/", sample, "_", platform, "_aligned.rds"))
    df_neg <- rbind(df_neg, data.frame(type = "control",
                                       value = apply(object@assays$ControlProbe$counts, 1, mean),
                                       platform=platform,
                                       sample=sample))
    df_neg <- rbind(df_neg,data.frame(type = "gene",
                                      value = apply(object@assays$Xenium@layers$counts, 1, mean),
                                      platform=platform,
                                      sample=sample)) 
    gc()
  }
}

write.table(df_neg, "data/output/neg_rate/neg_mean_v1v2.txt")