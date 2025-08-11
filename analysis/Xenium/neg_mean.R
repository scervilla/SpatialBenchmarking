library(dplyr)

path_objects <- "data/objects/"

samples <- c("Bladder", "Breast", "Lung", "Colorectal", "Lymphoma", "Ovarian")
df_neg <- data.frame()

for (sample in samples) {
  print(sample)
  
  # CosMx
  cosmx <- readRDS(paste0(path_objects, sample, "_alignment_cosmx.rds"))
  gex <- apply(cosmx@assays$RNA@layers$counts, 2, sum) %>% sum
  exp_neg <- tiledb_scdataset$somas$negprobes$X$members$counts$to_matrix(batch_mode = TRUE)[, colnames(cosmx)]
  negativex <- apply(exp_neg, 2, sum) %>% sum
  
  print(paste0("CosMx: ", (negativex / (negativex + gex)) *
                 (nrow(cosmx@assays$RNA@layers$counts) / nrow(exp_neg)) * 100))
  
  df_neg <- rbind(df_neg, data.frame(type = "control",
                                     value = apply(exp_neg, 1, mean),
                                     platform = "cosmx",
                                     sample = sample))
  df_neg <- rbind(df_neg, data.frame(type = "gene",
                                     value = apply(cosmx@assays$RNA@layers$counts, 1, mean),
                                     platform = "cosmx",
                                     sample = sample))
  
  # Xenium
  xenium <- readRDS(paste0(path_objects, sample, "_alignment_xenium.rds"))
  gex <- apply(xenium@assays$Xenium@layers$counts, 2, sum) %>% sum
  negativex <- apply(xenium@assays$ControlProbe$counts, 2, sum) %>% sum
  
  print(paste0("Xenium: ", (negativex / (negativex + gex)) *
                 (nrow(xenium@assays$Xenium@layers$counts) /
                    nrow(xenium@assays$ControlProbe$counts)) * 100))
  
  df_neg <- rbind(df_neg, data.frame(type = "control",
                                     value = apply(xenium@assays$ControlProbe$counts, 1, mean),
                                     platform = "xenium",
                                     sample = sample))
  df_neg <- rbind(df_neg, data.frame(type = "gene",
                                     value = apply(xenium@assays$Xenium@layers$counts, 1, mean),
                                     platform = "xenium",
                                     sample = sample))
}

write.table(df_neg, "data/output/neg_rate/neg_mean_cosmx_xenium.txt")