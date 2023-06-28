library(PLNmodels)

Y_full <- read.table("data_sets/borneo/data_abundance.txt", 
                     sep = ";", header = TRUE)
X_full <- read.table("data_sets/borneo/data_soil_characteristics_after_pca.txt", 
                     sep = ";", header = TRUE) %>% 
  arrange(Site) %>% 
  dummy_cols("Sol", remove_selected_columns = TRUE)
kept_sites <- pull(X_full, Site)
X <- select(X_full, -Site) %>% 
  as.matrix()
Y <- Y_full %>% 
  filter(Site %in% kept_sites) %>% 
  arrange(Site) %>% 
  select(-Site) %>% 
  as.matrix()

pln_data <- prepare_data(counts = Y_full %>% 
                           filter(Site %in% kept_sites) %>% 
                           arrange(Site) %>% 
                           column_to_rownames(var = "Site"), 
                         covariates = X_full %>% 
                           column_to_rownames(var = "Site"))

res_PLN <- PLNPCA(as.formula("Abundance ~ 1 + Sol_Dunaire + Sol_Grès +
                             Z1 + Z2 + Z3 + Z4"),
                  ranks = 1:3, data = pln_data)
plot(res_PLN)
myPCA_ICL <- getBestModel(res_PLN, "ICL") 
plot(myPCA_ICL)

sigma(myPCA_ICL) %>% cov2cor() %>%
  corrplot::corrplot(is.corr = TRUE, tl.cex = .1, 
                     hclust.method = "ward.D2", order = "hclust")

myPCA_ICL$model_par$Theta %>% 
  round(1) %>%  
  head()