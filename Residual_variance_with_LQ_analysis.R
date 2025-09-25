######### Author : Adeyemi Akinade
######### Place : Clemson University, SC
######### lab : Gopalan's lab
######### Date : September 9, 2025
######### Project : Comparative Epigenetics of Ageing
######### This script was written to determine the residual variance analysis of all cpg sites from DNAm vs aging

### Loading the slope data (age_in_years)
blood_regression_result <- readRDS("/Users/aakinad/Downloads/PHD FOLDER/blood_regression_result.rds")
head(blood_regression_result)

library(dplyr)
library(tidyr)
library(purrr)
# Suppose your data frame is blood_regression_result
# we can do residual variance directly
resid_var <- blood_regression_result %>%
  mutate(resid_var = map_dbl(residuals, ~ var(.x, na.rm = TRUE))) %>%
  select(species, cpg_site, resid_var)

blood_residual_variance_wt_LQ <- merge(species_with_mass_max_lifespan, resid_var, by = "species")

#### After merging the LQ and maximum lifespan information, the next is standard correlation approach

#  blood_residual_variance_wt_LQ has columns: species, CpG, var_residuals, LQ
cpg_sites <- unique(blood_residual_variance_wt_LQ$cpg_site)  
resid_var_cor_results <- data.frame()  # initialize empty dataframe

for (CpG_name in cpg_sites) {
  message(Sys.time(), " - Processing CpG: ", CpG_name)
  
  # Subset data for this CpG and complete cases
  residual_trait_cpg <- blood_residual_variance_wt_LQ %>%
    filter(cpg_site == CpG_name, !is.na(resid_var), !is.na(LQ))
  
 
    test <- cor.test(residual_trait_cpg$resid_var, residual_trait_cpg$LQ)
    
    resid_var_cor_results <- rbind(resid_var_cor_results, data.frame(
      CpG_site = CpG_name,
      cor = test$estimate,
      p_value = test$p.value
    ))
  } 

resid_var_cor_results$p_adjust <- p.adjust(resid_var_cor_results$p_value, method = "bonferroni")

###### write the residual variance correlation with LQ as a table
saveRDS(resid_var_cor_results, file = "blood_residual_variance_vs_LQ_standard_correlation.rds")

library(dplyr)
library(ggplot2)

ggplot(resid_var_cor_results, aes(x = cor)) +
  geom_histogram(aes(y = ..density..), bins = 40, 
                 fill = "skyblue", color = "white") +   # histogram
  geom_density(color = "red", size = 0.5) +            # density curve
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Distribution of correlations between CpG residual variance and LQ",
    x = "Correlation coefficient",
    y = "Density"
  ) +
  theme_minimal(base_size = 14)


# Add a flag for raw significance
resid_var_cor_results <- resid_var_cor_results %>%
  mutate(sig_flag = ifelse(p_value < 0.05, "significant", "NS"))

# Plot density with fill by significance
ggplot(resid_var_cor_results, aes(x = cor, fill = sig_flag)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("significant" = "tomato", "NS" = "skyblue")) +
  labs(
    title = "Distribution of residual variance at CpG sites correlations with LQ",
    x = "Correlation coefficient",
    y = "Density",
    fill = "Significance"
  ) +
  theme_minimal(base_size = 14)






#################### EVOLUTIONARY COVARIANCE AND CORRELATION BETWEEN residuals variance vs LQ
# Make sure species names match tree
library(devtools)
library(usethis)
library(phytools)
library(mvMORPH)

## loading the tree data
blood_tree <- read.tree("~/Downloads/blood_species_34.nwk")

## correcting for the unmatched species in tree_species and DNAm_rates_data
blood_residual_variance_wt_LQ$species <- gsub(" ", "_", blood_residual_variance_wt_LQ$species)
blood_residual_variance_wt_LQ$species <- gsub("Canis_lupus_familiaris", "Canis_lupus", blood_residual_variance_wt_LQ$species)
blood_residual_variance_wt_LQ$species <- gsub("Sus_scrofa", "Sus_scrofa_domesticus", blood_residual_variance_wt_LQ$species)
blood_residual_variance_wt_LQ$species <- gsub("", "Sus_scrofa_domesticus", blood_residual_variance_wt_LQ$species)


### There was chanllenges and error with the species match between blood tree and trait data
common_species <- setdiff(unique( blood_residual_variance_wt_LQ$species), blood_tree$tip.label)


cpg_sites <- unique(blood_residual_variance_wt_LQ$cpg_site)  # all CpGs
resid_var_results <- data.frame()  # initialize results

for(CpG_name in cpg_sites) {
  
  message(Sys.time(), " - Processing CpG: ", CpG_name)
  
  # Subset data for this CpG (already per species)
  resid_var_trait_cpg_species <- blood_residual_variance_wt_LQ %>%
    filter(cpg_site == CpG_name, !is.na(resid_var), !is.na(LQ))
  
  n_species <- nrow(resid_var_trait_cpg_species)
  if(n_species < 2) {
    message("Skipping CpG ", CpG_name, ": less than 2 species with data")
    next
  }
  
  # Trait matrix (use resid_var OR slope as CpG trait depending on your analysis)
  traits_tmp <- as.matrix(data.frame(
    CpG = resid_var_trait_cpg_species$resid_var,   # change to slope if needed
    LQ  = resid_var_trait_cpg_species$LQ
  ))
  rownames(traits_tmp) <- resid_var_trait_cpg_species$species
  
  # Fit mvBM
  fit <- tryCatch({
    mvBM(tree = blood_tree, data = traits_tmp, model = "BM1")
  }, error = function(e){
    message("CpG failed: ", CpG_name, " - ", e$message)
    return(NULL)
  })
  
  if(is.null(fit)) next
  
  # Covariance + correlation
  cor_matrix <- cov2cor(fit$sigma)
  
  resid_var_results <- rbind(resid_var_results, data.frame(
    CpG_site = CpG_name,
    cov      = fit$sigma[1,2],
    cor      = cor_matrix[1,2]
  ))
}

# Example: save results dataframe
saveRDS(results, file = "blood_results_DNAm_rate_LQ_mvBM.rds")




