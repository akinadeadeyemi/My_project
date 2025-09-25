######### Author : Adeyemi Akinade
######### Place : Clemson University, SC
######### lab : Gopalan's lab
######### Date : September 9, 2025
######### Project : Comparative Epigenetics of Ageing
######### This script was written to determine the correlation between DNAm rate vs LQ

############### MULTIVARIATE TRAIT EVOLUTION ######################
# Here I first determined the Longevity quotient (LQ). 
# The predicted lifespand was determined using the equation proposed by 
# Prothero & Jürgens for all mammals: 𝑦=  5.3 × 𝑚^0.174  
# All species mass (as adult weight in AnAge database) was measured in g but converted to kg.


# STEP 1: calculating LQ for all 34 species of the blood DNAm species
###  Loading the data
library(phytools)
library(foreach)
library(parallel)
library(doParallel)

library(readxl)
library(dplyr)
blood_species_with_mass_from_AnAge <- read_excel("/Users/aakinad/Downloads/blood_species_with_mass_from_AnAge.xls")
blood_metadata_wt_life_span <- readRDS("/Users/aakinad/Downloads/PHD FOLDER/blood_metadata_wt_life_span.rds")
blood_regression_result <- readRDS("/Users/aakinad/Downloads/PHD FOLDER/blood_regression_result.rds")
blood_tree <- read.tree("~/Downloads/blood_species_34.nwk")

## selecting the distinct species and their max. lifespan
species_lifespan <- blood_metadata_wt_life_span %>%
  distinct(species, max_lifespan, .keep_all = FALSE)

# column bind  blood_species_with_mass_from_AnAge and species_lifespan
species_with_mass_max_lifespan <- merge(species_lifespan,blood_species_with_mass_from_AnAge, by = "species")

# Calculate predicted lifespan using Prothero & Jürgens equation (1987))
species_with_mass_max_lifespan$predicted_lifespan <- 5.3 * (species_with_mass_max_lifespan$body_mass_kg ^ 0.174)

# Calculate LQ = observed / predicted
species_with_mass_max_lifespan$LQ <- species_with_mass_max_lifespan$max_lifespan / species_with_mass_max_lifespan$predicted_lifespan

# merge the LQ with the blood_metadata_wt_life_span

blood_metadata_wt_life_span_and_LQ <- merge(species_with_mass_max_lifespan, blood_metadata_wt_life_span,  by = "species")

# extracting the trait data for blood DNAm
blood_trait_data <- merge(species_with_mass_max_lifespan, blood_regression_result, by = "species")



################ USING FOR LOOP FOR MULTIVARIATE TRAIT EVOLUTION ######################
library(dplyr)

cpg_sites <- unique(blood_trait_data$cpg_site)  # all CpGs
cor_results <- data.frame()  # initialize empty dataframe

for(CpG_name in cpg_sites) {
  
  message(Sys.time(), " - Processing CpG: ", CpG_name)
  
  # Subset data for this CpG and complete cases
  trait_cpg <- blood_trait_data %>%
    filter(cpg_site == CpG_name) %>%
    filter(!is.na(slope) & !is.na(LQ) & !is.na(std_error))
  
  if(nrow(trait_cpg) > 1) {  # Need at least 2 points to compute correlation
    test <- cor.test(trait_cpg$slope, trait_cpg$LQ)
    
    cor_results <- rbind(cor_results, data.frame(
      CpG_site = CpG_name,
      cor = test$estimate,
      p_value = test$p.value
    ))
  } else {
    message("Skipping CpG ", CpG_name, ": not enough data points")
  }
}

### getting the adjusted p-values
resid_var_cor_results$p.adjust <- p.adjust(resid_var_cor_results$p_value, method = "bonferroni")

# Example: save results dataframe

saveRDS(cor_results, file = "blood_results_DNAm_rate_LQ_normal_correlation.rds")

##### comapring evolutionary correlations with standard correlations
blood_results_DNAm_rate_LQ_mvBM <- readRDS("/Users/aakinad/Downloads/PHD FOLDER/Evolution of Epigenetic Aging across mammalian species /Blood_tissues/blood_results_DNAm_rate_LQ_mvBM.rds")

evo_std_comparison <- data.frame(
  CpG = cor_results$CpG_site,
  cor_standard = cor_results$cor,
  cor_evolutionary = blood_results_DNAm_rate_LQ_mvBM$cor
)

library(ggplot2)

ggplot(evo_std_comparison, aes(x = cor_standard, y = cor_evolutionary)) +
  geom_point(alpha = 0.6, color = "blue") +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(
    title = "Comparison of Standard vs Evolutionary Correlations",
    x = "Standard correlation (within-species)",
    y = "Evolutionary correlation (across species)"
  ) +
  theme_minimal(base_size = 14)

cor(evo_std_comparison$cor_standard, evo_std_comparison$cor_evolutionary, use = "complete.obs")

library(ggplot2)
ggplot(evo_std_comparison, aes(x = cor_standard, y = cor_evolutionary)) +
  geom_point(alpha = 0.05, color = "blue") +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(
    title = "Comparison of Standard vs Evolutionary Correlations",
    x = "Standard correlation (across species)",
    y = "Evolutionary correlation (across species)"
  ) +
  annotate("text",
           x = min(evo_std_comparison$cor_standard, na.rm = TRUE),
           y = max(evo_std_comparison$cor_evolutionary, na.rm = TRUE),
           label = "corr = 0.875",
           hjust = 0, vjust = 1.2, size = 5, color = "black") + theme(
             axis.title.x = element_text(size = 20, face = "bold"),
             axis.title.y = element_text(size = 20, face = "bold"),
             axis.text = element_text(size = 20)
           )



library(ggplot2)
library(dplyr)
library(tidyr)


# Reshape into long format
evo_long <- evo_std_comparison %>%
  pivot_longer(
    cols = c(cor_standard, cor_evolutionary),
    names_to = "Correlation_Type",
    values_to = "Correlation"
  )

# Plot histogram + density outline
ggplot(evo_long, aes(x = Correlation, fill = Correlation_Type, color = Correlation_Type)) +
  geom_histogram(alpha = 0.4, position = "identity", bins = 40) +
  geom_density(size = 1.2) +
  labs(
    title = "Distribution of  Evolutionary and Standard Correlations",
    x = "Correlation coefficient",
    y = "Count / Density",
    fill = "Correlation type",
    color = "Correlation type"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

library(ggplot2)

ggplot(evo_std_comparison) +
  geom_density(aes(x = cor_standard, color = "Standard"), size = 1.2) +
  geom_density(aes(x = cor_evolutionary, color = "Evolutionary"), size = 1.2) +
  scale_fill_manual(values = c("Standard" = "blue", "Evolutionary" = "green")) +
  scale_color_manual(values = c("Standard" = "blue", "Evolutionary" = "green")) +
  labs(
    title = "Distribution of Evolutionary and Standard Correlations at all CpG sites",
    x = "Correlation",
    y = "Density"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 20)
  )


