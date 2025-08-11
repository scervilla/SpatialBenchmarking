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

path <- "data/raw/Xenium/"
path_output <- "data/output/nbins/"

transcripts <- read_parquet(paste0(path, xenium_folder[sample], "/transcripts.parquet"))
transcripts <- transcripts %>% filter(qv >= 20 & cell_id != "UNASSIGNED")

df <- read.table(paste0("data/output/df/df_", sample, ".txt"))  %>% tibble::rownames_to_column(var = "cell_id") %>% 
  filter(platform == "xenium" & !is.na(ct)) %>% 
  select(cell_id, ct)

transcripts <- inner_join(transcripts, df, by="cell_id")

sizes <- c(2,8,16)

for (size in sizes) {
  for (c in unique(df$ct)) {
    print(c)
    transcripts_grouped <-  transcripts %>% 
      filter(c == ct) %>% 
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
    write.table(n_cells_grid, paste0(path_output, "n_cells_grid_", c, "_", size, "um_", sample, ".txt"))
    write.table(n_grid_cells, paste0(path_output, "n_grid_cells_", c, "_", size, "um_", sample, ".txt"))
  }
}
