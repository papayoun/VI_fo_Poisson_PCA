rm(list = ls())


# librairies ---------------------------------------------------------------

library(tidyverse)
library(abind) # To gather results together
library(parallel) # For some parallel computations
library(pacman) # For progress bar
library(torch)
library(fastDummies)
library(ggcorrplot) # For correlation plots
library(ellipse) # For ellipse plots

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1 - cormat) / 2)
  hc <- hclust(dd, method = "ward.D2")
  cormat[hc$order, hc$order]
}

# Graphics theme ----------------------------------------------------------

my_theme <- function(base_size = 8){
  theme_bw()  %+replace%
    theme(
      panel.border = element_rect(colour = "black", 
                                  fill = rgb(0, 0, 0, 0)),
      # plot.background = element_rect(fill = "white"),# bg around panel
      legend.background = element_blank(), 
      text = element_text(family = "LM Roman 10", size = base_size),
      axis.title = element_text(size = rel(1)),
      legend.text = element_text(size = rel(1)),
      legend.title = element_text(size = rel(1)),
      axis.text = element_text(size = rel(1)),
      strip.background = element_rect(fill = "lightgoldenrod1",
                                      color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size = rel(1)),
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5))
}

theme_set(my_theme(base_size = 16))


# Data --------------------------------------------------------------------

source("utils_load_main_application_data.R")

# Residual correlation on model without Sol -------------------------------

results_wo_sol <- load("result_VI_full_CAVI_on_borneo_without_sol.RData") %>% 
  get()
results_complete <- load("result_VI_full_CAVI_on_borneo.RData") %>% 
  get()

estimated_cor_wo_sol <- (t(results_wo_sol$params$Lambda$M) %*% 
                    results_wo_sol$params$Lambda$M) %>%
  {. + diag(results_wo_sol$params$Sigma$A / 
              results_wo_sol$params$Sigma$B)} %>% 
  cov2cor() %>% 
  reorder_cormat()

chosen_species_indexes <- which(abs(estimated_cor_wo_sol)  > .6, arr.ind = TRUE) %>% 
  as.data.frame() %>% 
  filter(row != col) %>% # Remove diagonal 
  pull(row) %>% 
  unique()
chosen_species <- colnames(Y_complete)[chosen_species_indexes]

example_cor_wo_sol <- estimated_cor_wo_sol[chosen_species_indexes, 
                                    chosen_species_indexes] %>%
  matrix(nrow = nrow(.), ncol = ncol(.),
         dimnames = list(str_replace_all(chosen_species, "_", " "),
                         str_replace_all(chosen_species, "_", " "))) %>% 
  reorder_cormat()

ggcorrplot(example_cor_wo_sol[,rev(1:nrow(example_cor_wo_sol))], 
           show.diag = TRUE, 
           legend.title = "Residual correlation") +
  theme(text = element_text(face = "italic"),
        legend.title = element_text(face = "plain"))


ggcorrplot(estimated_cor_wo_sol[, rev(1:nrow(estimated_cor_wo_sol))],
           legend.title = "Residual correlation",
           title = "Complete residual correlation matrix") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

# Plot covariates parameters posterior ellipses ---------------------------

covariate_indexes <- c(1, 4)
conf_ellipses_wo_sol <- map_dfr(1:length(colnames(Y_complete)),
                               function(j){
                                 ellipse(x = results_wo_sol$params$Beta$Cov[covariate_indexes, covariate_indexes, j], 
                                         center = results_wo_sol$params$Beta$M[covariate_indexes, j]) %>% 
                                   as.data.frame() %>% 
                                   mutate(Species = colnames(Y_complete)[j],
                                          C_x = results_wo_sol$params$Beta$M[covariate_indexes[1], j],
                                          C_y = results_wo_sol$params$Beta$M[covariate_indexes[2], j])
                               })

filter(conf_ellipses_wo_sol, Species %in% chosen_species) %>% 
  mutate(Species = str_replace_all(Species, "_", " ")) %>% 
  ggplot(aes(x = x, y = y, group = Species)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(aes(x = C_x, y = C_y, color = Species)) +
  geom_polygon(aes(fill = Species), alpha = .5) +
  labs(x = "Reaction to available cations",
       y = "Reaction to available NH4",
       title = "95% credible ellipses") +
  theme(legend.text = element_text(face = "italic"))


# Figure in the complete results case -------------------------------------


estimated_cor_complete <- (t(results_complete$params$Lambda$M) %*% 
                           results_complete$params$Lambda$M) %>%
  {. + diag(results_complete$params$Sigma$A / 
              results_complete$params$Sigma$B)} %>% 
  cov2cor() %>% 
  reorder_cormat()

example_cor_complete <- estimated_cor_complete[chosen_species_indexes, 
                                           chosen_species_indexes] %>%
  matrix(nrow = nrow(.), ncol = ncol(.),
         dimnames = list(str_replace_all(chosen_species, "_", " "),
                         str_replace_all(chosen_species, "_", " "))) %>% 
  {.[rownames(example_cor_wo_sol),
     colnames(example_cor_wo_sol)]}

ggcorrplot(example_cor_complete[,rev(1:nrow(example_cor_complete))], 
           show.diag = TRUE, 
           legend.title = "Residual correlation, with soil typology") +
  theme(text = element_text(face = "italic"),
        legend.title = element_text(face = "plain"))

ggcorrplot(estimated_cor_complete[, rev(1:nrow(estimated_cor_complete))],
           legend.title = "Residual correlation",
           title = "Complete residual correlation matrix, with soil typology") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))


conf_ellipses_complete <- map_dfr(1:length(colnames(Y_complete)),
                                  function(j){
                                    ellipse(x = results_complete$params$Beta$Cov[covariate_indexes, covariate_indexes, j], 
                                            center = results_complete$params$Beta$M[covariate_indexes, j]) %>% 
                                      as.data.frame() %>% 
                                      mutate(Species = colnames(Y_complete)[j],
                                             C_x = results_complete$params$Beta$M[covariate_indexes[1], j],
                                             C_y = results_complete$params$Beta$M[covariate_indexes[2], j])
                                  })


filter(conf_ellipses_complete, Species %in% chosen_species) %>% 
  mutate(Species = str_replace_all(Species, "_", " ")) %>% 
  ggplot(aes(x = x, y = y, group = Species)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(aes(x = C_x, y = C_y, color = Species)) +
  geom_polygon(aes(fill = Species), alpha = .5) +
  labs(x = "Reaction to available cations",
       y = "Reaction to available NH4",
       title = "95% posterior ellipses, with soil typology") +
  theme(legend.text = element_text(face = "italic"))

