######### Author : Adeyemi Akinade
######### Place : Clemson University, SC
######### lab : Gopalan's lab
######### Date : August 13, 2025
######### Project : Comparative Epigenetics of Ageing
######### This script was written to determine the multivariate trait evolution 

############### MULTIVARIATE TRAIT EVOLUTION ######################
# Here I first determined the Longevity quotient (LQ). 
# The predicted lifespand was determined using the equation proposed by 
# Prothero & Jürgens for all mammals: 𝑦=  5.3 × 𝑚^0.174  
# All species mass (as adult weight in AnAge database) was measured in g but converted to kg.


# STEP 1: calculating LQ for all 34 species of the blood DNAm species
###  Loading the data
library(readxl)
library(dplyr)
blood_species_with_mass_from_AnAge <- read_excel("Downloads/blood_species_with_mass_from_AnAge.xls")
blood_metadata_wt_life_span <- readRDS("~/Downloads/blood_metadata_wt_life_span.rds")
blood_regression_result <- readRDS("~/Downloads/blood_regression_result.rds")
blood_tree <- read.tree("~/Downloads/blood_species_34.nwk")

## selecting the distinct species and their max. lifespan
species_lifespan <- blood_metadata_wt_life_span %>%
  distinct(species, max_lifespan, .keep_all = FALSE)

# column bind  blood_species_with_mass_from_AnAge and species_lifespan
species_with_mass_max_lifespan <- merge(species_lifespan,blood_species_with_mass_from_AnAge, by = "species")

# Calculate predicted lifespan using Prothero & Jürgens equation
species_with_mass_max_lifespan$predicted_lifespan <- 5.3 * (species_with_mass_max_lifespan$body_mass_kg ^ 0.174)

# Calculate LQ = observed / predicted
species_with_mass_max_lifespan$LQ <- species_with_mass_max_lifespan$max_lifespan / species_with_mass_max_lifespan$predicted_lifespan

# merge the LQ with the blood_metadata_wt_life_span

blood_metadata_wt_life_span_and_LQ <- merge(species_with_mass_max_lifespan, blood_metadata_wt_life_span,  by = "species")

# extracting the trait data for blood DNAm
blood_trait_data <- merge(species_with_mass_max_lifespan, blood_regression_result, by = "species")


### Running parallel Analysis
library(phytools)
library(foreach)
library(parallel)
library(doParallel)

################ USING FOR LOOP FOR MULTIVARIATE TRAIT EVOLUTION ######################
library(phytools)
## correcting for the unmatched species in tree_species and DNAm_rates_data
blood_trait_data$species <- gsub("Canis_lupus_familiaris", "Canis_lupus", blood_trait_data$species)
blood_trait_data$species <- gsub("Sus_scrofa_domesticus_domesticus", "Sus_scrofa_domesticus", blood_trait_data$species)


# Make sure species names match tree
library(devtools)
library(usethis)
library(mvMORPH)

blood_trait_data$species <- gsub(" ", "_", blood_trait_data$species)

cpg_sites <- unique(blood_trait_data$cpg_site)  # all CpGs
results <- data.frame()  # initialize empty dataframe


for(CpG_name in cpg_sites) {
  
  message(Sys.time(), " - Processing CpG: ", CpG_name)
  
  # Subset data for this CpG and complete cases
  trait_cpg <- blood_trait_data %>%
    filter(cpg_site == CpG_name) %>%
    filter(!is.na(slope) & !is.na(LQ) & !is.na(std_error))
  
  # Aggregate multiple individuals per species
  trait_cpg_species <- trait_cpg %>%
    group_by(species) %>%
    summarise(
      slope = mean(slope, na.rm = TRUE),
      LQ = mean(LQ, na.rm = TRUE),
      SE = sqrt(mean(std_error^2, na.rm = TRUE))
    ) %>%
    ungroup()
  
  n_species <- nrow(trait_cpg_species)
  if(n_species < 2) {
    message("Skipping CpG ", CpG_name, ": less than 2 species with data")
    next
  }
  
  # Trait matrix
  traits_tmp <- data.frame(
    CpG = trait_cpg_species$slope,
    LQ  = trait_cpg_species$LQ
  )
  rownames(traits_tmp) <- trait_cpg_species$species
  
  # Fit mvBM (Brownian motion)
  fit <- tryCatch({
    mvBM(tree = blood_tree, data = traits_tmp, model = "BM1")
  }, error = function(e){
    message("CpG failed: ", CpG_name, " - ", e$message)
    return(NULL)
  })
  
  if(is.null(fit)) next
  
  # Covariance and correlation
  cor_matrix <- cov2cor(fit$sigma)
  
  results <- rbind(results, data.frame(
    CpG_site = CpG_name,
    cov = fit$sigma[1,2],
    cor = cor_matrix[1,2]
  ))
}

# Example: save results dataframe
saveRDS(results, file = "blood_results_DNAm_rate_LQ_mvBM.rds")


library(phytools)

# Step 1: Get numeric vector of LQ with species names
LQ_vec <- blood_trait_data %>%
  distinct(species, LQ) %>% 
  { setNames(.$LQ, .$species) }   # creates named numeric vector

# Step 2: Keep only species present in both tree and vector
common_species <- intersect(blood_tree$tip.label, names(LQ_vec))
LQ_vec <- LQ_vec[common_species]
blood_tree <- drop.tip(blood_tree, setdiff(blood_tree$tip.label, common_species))

# Step 3: Reorder vector to match tree tip order
LQ_vec <- LQ_vec[blood_tree$tip.label]

# Step 4: Plot contMap
LQ_plot <- contMap(blood_tree, LQ_vec, plot = FALSE)
plot(LQ_plot)

## BAR PLOT OF LQ
# Define color palette
custom_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF00FF", "#838B8B", "#0000FF", "#050505", "#00CD00", "#FF1493",
  "#FFC1C1", "#B03060", "#2F4F4F", "#551A8B", "#FF7F00", "#A65628", "#00F5FF", "#F781BF", "#CD8500", "#66C2A5",
  "#7CFC00", "#0000CD", "#8B2323", "#6B8E23", "#87CEFF", "#A4D3EE", "#CD2990", "#FC8D62", "#458B74", "#CDBA96",
  "#8B3626", "#006400", "#FFD700", "#8B6508", "#00008B", "#9370DB", "#96CDCD", "#8B4500", "#2F4F4F", "#87CEFA"
)
ggplot(LQ, aes(x = LQ, fill = species)) +
  geom_histogram() +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() + #theme(legend.position = "none") +
  labs(title = "Distribution of LQ by Species",
       x = "Longevity Quotient (LQ)",
       y = "Count")


########### How many cpg in phylogenetic signal show up in LQ/CpG evolutionary covariation #############
phylo_cpg <- blood_phylo_lambda_sig %>% select(cpg_site)

phylo_LQ_cov <- results %>% filter(CpG_site %in% phylo_cpg$cpg_site) 





