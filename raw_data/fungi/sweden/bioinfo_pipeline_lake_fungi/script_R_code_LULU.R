.libPaths(c("/cluster/home/pedromm/R/3_5_1",.libPaths()))

library(tidyverse)
library(openxlsx)
library(lulu)
path<-getwd()
library(tidyverse)
OTU_table_path <- list.files(path, pattern = ".otutable", full.names = TRUE)
OTU_table<-read.csv(OTU_table_path, sep="\t")
OTU_table2<-OTU_table %>% 
  separate(OTUId,c("OTUid"),";") %>%
  column_to_rownames("OTUid")

match_list_path<-list.files(path, pattern = "list.txt", full.names = TRUE)
match_list<-read.csv(match_list_path, sep="\t")
lulu_curation<-lulu(OTU_table2, match_list)
OTU_table_lulu_curated<-lulu_curation$curated_table %>% 
  rownames_to_column("OTUid")
write.xlsx(OTU_table_lulu_curated,"OTU_table_lulu_curated.xlsx")
saveRDS(lulu_curation, "lulu_analyses.R")
