library(Seurat)
library(dplyr)

path_object <- "data/objects/"

args = commandArgs(trailingOnly=TRUE)
sample <- args[1]


df <- data.frame(row.names  = c("Bladder", "Colorectal", "Lymphoma", "Ovarian"),
                 xmin = c(0, 0, 1500, 1000),
                 xmax = c(6500, 5000, 8000, 7500),
                 ymin = c(0, 0, 2550, 400),
                 ymax = c(6500, 6500, 9050, 6900))

xenium <- readRDS(paste0(path_object, "seurat/", sample, "_seurat.rds"))

if (sample == "Lung") {
  coords <- xenium@images$fov@boundaries$centroids@coords %>% as.data.frame()
  rownames(coords) <- colnames(xenium)
  
  isin <- sapply(1:nrow(coords), function(i) {
    m1 <- 5.5
    m2 <- -1/m1
    y <- coords[i,"x"]
    x <- coords[i,"y"]
    xright <- (y+45000)/m1
    xleft <- (y+9500)/m1
    ymin <- x*m2+2200
    ymax <- x*m2+9000
    return(x >= xleft & x <= xright &
             y >= ymin & y <= ymax)
  })
  cells <- colnames(xenium)[isin]
  
} else {
  
  coords <- xenium@images$fov@boundaries$centroids@coords 
  rownames(coords) <- colnames(xenium)
  
  xmin = df[sample, "xmin"]
  ymax = df[sample, "ymax"]
  ymin = df[sample, "ymin"]
  xmax = df[sample, "xmax"]
  
  cells <- coords  %>% as.data.frame() %>% dplyr::filter( xmin <= x & xmax >= x &
                                                            ymin <= y & ymax >= y) %>% rownames()
}

xenium <- subset(xenium, cells=cells) 
xenium <- SCTransform(xenium, assay = "Xenium")


saveRDS(xenium, paste0(path_object, "seurat_aligned/", sample, "_xenium_in_hd.rds"))
