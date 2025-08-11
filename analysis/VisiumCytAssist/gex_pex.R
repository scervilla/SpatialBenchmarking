library(stringr)
library(dplyr)
library(circlize)
library(Seurat)
library(ComplexHeatmap)

path_object <- "data/objects/seurat_aligned/"
path_output <- "data/output/"

for (sample in c("Bladder", "Lung", "Colorectal", "Lymphoma", "Ovarian")) {
  print(sample)
  seurat <- readRDS(paste0(path_object, sample, "_cyt_xenium.rds"))
  
  #Normalize
  DefaultAssay(seurat) <- "Protein"
  seurat <- subset(seurat, features = (seurat@assays$Protein %>% rownames())[-c(15:18)])
  seurat <- NormalizeData(seurat, assay = "Protein", normalization.method = "CLR")
  seurat <- SCTransform(seurat, assay = "pseudo", new.assay.name = "SCTpseudo")
  
  proteins <- seurat@assays$Protein %>% rownames() %>% str_replace(pattern = "\\.1|\\.2", replacement = "")
  genes <- seurat@assays$RNA %>% rownames()
  xenium_genes <- seurat@assays$pseudo %>% rownames() %>% str_replace(pattern = "-pseudo", replacement = "")
  
  #PEX x GEX
  int_gen_prot <- intersect(proteins, genes)
  
  cor_df <- data.frame()
  for(gene in int_gen_prot) {
    cor_df <- rbind(cor_df, c(gene, cor(seurat@assays$Protein$data[paste0(gene,".1"),], 
                                        seurat@assays$SCT$data[gene,], method = "spearman")) )
  }
  colnames(cor_df) <- c("gene", "cor")
  write.table(cor_df, paste0(path_output, "gex_pex/" ,sample, "_cor_gex_pex.txt"))
  
  #PEX x Xenium
  int_gen_prot <- intersect(proteins, xenium_genes)
  
  cor_df <- data.frame()
  for(gene in int_gen_prot) {
    cor_df <- rbind(cor_df, c(gene, cor(seurat@assays$Protein$data[paste0(gene,".1"),], 
                                        seurat@assays$SCTpseudo$data[paste0(gene,"-pseudo"),], method = "spearman")) )
  }
  colnames(cor_df) <- c("gene", "cor")
  write.table(cor_df, paste0(path_output, "gex_pex/" ,sample, "_cor_xenium_pex.txt"))
  
  #PEX x CT
  int_gen_prot <- intersect(genes, proteins)
  cor_df <- data.frame()
  for(gene in int_gen_prot) {
    for (ct in colnames(cytassist@meta.data)[8:13]) {
      cor_df <- rbind(cor_df, c(gene, ct, cor(seurat@assays$Protein$data[paste0(gene,".1"),], 
                                              seurat@meta.data[,ct], method = "spearman")) )
    }
  }
  colnames(cor_df) <- c("gene","cell_type", "cor")
  cor_df[,3] <- as.numeric(cor_df[,3])
  write.table(cor_df, paste0(path_output, "gex_pex/" ,sample, "_cor_ct_pex.txt"))
}
