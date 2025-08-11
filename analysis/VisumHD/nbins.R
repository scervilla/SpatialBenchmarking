library(arrow)
library(dplyr)


args = commandArgs(trailingOnly=TRUE)
sample <- args[1]

xenium_folder <- c( "Bladder" = "output-XETG00051__0003605__S150192__20230922__150140",
                    "Lung" = "output-XETG00051__0003605__S150193__20230922__150140",
                    "Breast" = "output-XETG00051__0003605__S150194__20230922__150140",
                    "Colorectal" = "output-XETG00051__0003606__S150195__20230922__150140",
                    "Lymphoma" = "output-XETG00051__0003606__S150196__20230922__150140",
                    "Ovarian" = "output-XETG00051__0003606__S150197__20230922__150140")

path_data <- "data/raw/Xenium/"
path_output <- "data/output/nbins/"

transcripts <- read_parquet(paste0(path_data, xenium_folder[sample], "/transcripts.parquet"))
transcripts <- transcripts %>% filter(qv >= 20 & cell_id != "UNASSIGNED")

sizes <- c(2,8,16)

for (size in sizes) {
    print(c)
    transcripts_grouped <-  transcripts %>% 
      mutate(
        x_grid = floor(x_location / size),
        y_grid = floor(y_location / size),
        grid = paste0(x_grid, "_", y_grid)
      ) %>%
      distinct(cell_id, grid)
    
    n_cells_grid <- transcripts_grouped %>% group_by(grid) %>%
      count
    
    n_grid_cells <- transcripts_grouped  %>% group_by(cell_id) %>%
      count
    write.table(n_cells_grid, paste0(path_output, "n_cells_grid_",  size, "um_", sample, ".txt"))
    write.table(n_grid_cells, paste0(path_output, "n_grid_cells_", size, "um_", sample, ".txt"))
}
