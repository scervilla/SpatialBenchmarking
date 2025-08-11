library(Seurat)
library(decoupleR)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
sample <- args[1]

net <- decoupleR::get_progeny(organism = 'human', top = 500)

mat <- readRDS(paste0("data/objects/seurat/", sample, "_HD_sct.rds"))


acts <- decoupleR::run_mlm(mat = mat, 
                           net = net, 
                           .source = 'source', 
                           .target = 'target',
                           .mor = 'weight', 
                           minsize = 5)

acts <- acts %>% tidyr::pivot_wider(id_cols = 'source', 
                                    names_from = 'condition',
                                    values_from = 'score') %>%
  tibble::column_to_rownames(var = 'source') 

write.table(acts, paste0("data/outputs/progeny/progeny_", sample, "_HD.txt"))
