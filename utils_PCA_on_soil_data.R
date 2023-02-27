rm(list = ls())

# Librairies --------------------------------------------------------------

library(tidyverse) # For data manipulation
library(FactoMineR) # For PCA
library(factoextra) # For PCA graphs

# Loading data ------------------------------------------------------------

soil_raw <- read.table("data_sets/borneo/data_soil_characteristics.txt",
                        sep = ";", header = TRUE)
abundance_raw <- read.table("data_sets/borneo/data_abundance.txt",
                             sep = ";", header = TRUE)

# Get transformed data for PCA ------------------------------------------

soil_for_pca <- select_at(soil_raw, 
                          vars(Site, Sol, contains("_Superficiel"))) %>% 
  na.omit() %>% 
  select(-Limon_Superficiel)


# Ajustement de l'ACP -----------------------------------------------------

PCA_result <- PCA(select_if(soil_for_pca, is.numeric), graph = FALSE)


# Variance des composantes ------------------------------------------------

fviz_eig(PCA_result) +
  labs(y = "Pourcentage de variance expliquée", 
       title = "Variance expliquée par les axes principaux")

# Cercle des corrélations -------------------------------------------------

fviz_pca_var(PCA_result, axes = c(1, 2))   
fviz_pca_var(PCA_result, axes = c(3, 4))  


# Représentation des individus --------------------------------------------

fviz_pca_ind(PCA_result, axes = c(1, 2), col.ind = soil_for_pca$Sol)  +
  labs(title = "Représentation des individus",
       color = "Type de sol", shape = "Type de sol")
fviz_pca_ind(PCA_result, axes = c(3, 4))  


# Extraction d'un tableau avec toute les nouvelles variables --------------


PCA_result$ind$coord %>% 
  as_tibble() %>% 
  select(1:4) %>% 
  mutate(Site = soil_for_pca$Site,
         Sol = soil_for_pca$Sol) %>% 
  rename_at(vars(contains("Dim.")), 
            function(nom) str_replace(nom, "Dim.", "Z")) %>% 
  write.table(file = "data_sets/borneo/data_soil_characteristics_after_pca.txt",
              sep = ";", col.names = TRUE, row.names = FALSE)

