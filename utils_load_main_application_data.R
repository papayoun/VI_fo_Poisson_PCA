Y_full <- read.table("data_sets/borneo/data_abundance.txt", 
                     sep = ";", header = TRUE)
X_full <- read.table("data_sets/borneo/data_soil_characteristics_after_pca.txt", 
                     sep = ";", header = TRUE) %>% 
  arrange(Site) %>% 
  dummy_cols("Sol", remove_selected_columns = TRUE)
kept_sites <- pull(X_full, Site)
Y_complete <- Y_full %>% 
  filter(Site %in% kept_sites) %>% 
  arrange(Site) %>% 
  select(-Site) %>% 
  as.matrix()