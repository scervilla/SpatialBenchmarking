library(SingleR)
library(celldex)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
sample <- args[1]

hpca.se <- HumanPrimaryCellAtlasData()
wanted <- c("DC",   "Epithelial_cells",     "B_cell",              "Neutrophils" ,
            "T_cells" ,  "Monocyte" , "Endothelial_cells",  "Macrophage",          "NK_cell",
            "Fibroblasts")
hpca.se <- hpca.se[,hpca.se$label.main %in% wanted]
hpca.se$label <- colData(hpca.se) %>% as.data.frame() %>% 
  filter(label.main %in% wanted) %>%
  mutate(labels = label.main) %>% rowwise() %>% 
  mutate(labels = case_when(label.main == "T_cells" ~ label.fine,
                            label.main == "B_cell" ~ label.fine,
                            TRUE ~ label.main)) %>% pull(labels) 

main_path <- "data/objects/"
platform <- "cosmx"
hESCs <- readRDS(paste0(main_path,"/seurat_aligned/", sample, "_alignment_", platform, ".rds"))@assays$SCT@data
pred.hesc <- SingleR(test = hESCs, ref = hpca.se, assay.type.test=1,
                     labels = hpca.se$label)
saveRDS(pred.hesc, paste0("data/output/SingleR/predictions_", sample, "_", platform, "_final.rds"))

