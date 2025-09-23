######### Author : Adeyemi Akinade
######### Place : Clemson University, SC
######### lab : Gopalan's lab
######### Date : March 19, 2025
######### Project : Comparative Epigenetics of Ageing
######### This script was written to determine the linear regression between age and beta values of DNA methylation sites with phylogenetic frame work


######################## PHYLOGENETIC LINEAR REGRESSION IN DNA METHYLATION BLOOD SAMPLES ONLY ##########################################
# Incorporating phylogeny into the model ----------------------------------
### Loading the necessary libraries
library(dplyr)
library(phylolm)
library(ape)

## working on the loading datasets
blood_species_tree <- read.tree("~/Phd_data/blood_species_tree.newick")
blood_metadata_wt_life_span <- readRDS("~/Phd_data/blood_metadata_wt_life_span.rds")
blood_metadata_wt_life_span_betas <- readRDS("~/Phd_data/blood_metadata_wt_life_span_betas.rds")

## processing the tree data
# Remove the OTL suffix to match metadata species
clean_labels <- gsub("_ott\\d+", "", blood_species_tree$tip.label)

# Replace underscores with spaces (if your metadata uses spaces)
clean_labels <- gsub("_", " ", clean_labels)

# Assign cleaned labels back to the tree
blood_species_tree$tip.label <- clean_labels
blood_species_tree$tip.label 

######### Author : Adeyemi Akinade
######### Place : Clemson University, SC
######### lab : Gopalan's lab
######### Date : June 16, 2025
######### Project : Comparative Epigenetics of Ageing
######### This script was written to determine the linear regression between age and beta values of DNA methylation sites in tissues


######################## PHYLOGENETIC REGRESSION MODEL IN SKIN SAMPLES ONLY ##########################################

library(dplyr)
library(phylolm)
library(ape)

## working on the loading skin species tree
skin_species_tree <- read.tree("~/Phd_data/skin_species.nwk")
View(skin_species_tree$tip.label)
is.ultrametric(skin_species_tree)   # Should return TRUE
is.binary(skin_species_tree)        # Should return TRUE
plot(skin_species_tree, cex = 0.7)


### Loading other datasets
skin_metadata_wt_life_span <- readRDS("~/Phd_data/skin_metadata_wt_life_span.rds")
skin_metadata_wt_life_span_betas <- readRDS("~/Phd_data/skin_metadata_wt_life_span_betas.rds")
skin_regression_wt_lifespan_result <- readRDS("~/Phd_data/outputs/skin_regression_wt_lifespan_result.rds")

### Extract the betip.label### Extract the beta values  (DNAm rates) from the age vs DNAm regression
DNAm_rates_with_age <- skin_regression_wt_lifespan_result [ , c("species", "cpg_site", "slope")]
DNAm_rates_std_err <- skin_regression_wt_lifespan_result [ , c("species", "cpg_site", "std_error")]

### Reshape the slope data 
DNAm_rates_with_age_wide <- DNAm_rates_with_age %>%
  tidyr::pivot_wider(
    names_from = species,
    values_from = slope,
    id_cols = cpg_site 
  ) %>%
  tibble::column_to_rownames(var = "cpg_site")

View(DNAm_rates_with_age_wide)

### Reshape the slope data 
DNAm_rates_std_err_wide <- DNAm_rates_std_err %>%
  tidyr::pivot_wider(
    names_from = species,
    values_from = std_error,
    id_cols = cpg_site 
  ) %>%
  tibble::column_to_rownames(var = "cpg_site")

View(DNAm_rates_std_err_wide)

## Write file as RDS
saveRDS(DNAm_rates_with_age_wide, file = "DNAm_rates_with_age_wide.rds")
saveRDS(DNAm_rates_std_err_wide, file = "DNAm_rates_std_err_wide.rds")

library(dplyr)
colnames(DNAm_rates_with_age_wide)
### adding undersore to fit the phylogenetic tree style
colnames(DNAm_rates_with_age_wide) <- gsub(" ", "_", colnames(DNAm_rates_with_age_wide))
colnames(DNAm_rates_std_err_wide) <- gsub(" ", "_", colnames(DNAm_rates_std_err_wide))


### Working on the phylosignal


##Setting the difference between the tree_species and DNAm_rates_with_age_wide datasets.
setdiff(skin_species_tree$tip.label, colnames(DNAm_rates_with_age_wide))
setdiff(skin_species, colnames(DNAm_rates_std_err_wide))
## correcting for the unmatched species in tree_species and DNAm_rates_with_age_wide
colnames(DNAm_rates_with_age_wide)[colnames(DNAm_rates_with_age_wide) == "Cephalorhynchus_hectori_hectori"] <- "Cephalorhynchus_hectori"
colnames(DNAm_rates_std_err_wide)[colnames(DNAm_rates_std_err_wide) == "Cephalorhynchus_hectori_hectori"] <- "Cephalorhynchus_hectori"

# Ensure species order matches tree
slope_matrix <- slope_matrix[, skin_species_tree$tip.label]
se_matrix <- se_matrix[, skin_species_tree$tip.label]

colnames(se_matrix)
skin_species_tree$tip.label

## load phytools
library(phytools)

# Ensure species order matches tree
species_order <- skin_species_tree$tip.label

lambda_with_se <- function(tree, slope_matrix, se_matrix) {
  library(phytools)
  
  species_order <- tree$tip.label
  slope_matrix <- slope_matrix[, species_order]
  se_matrix <- se_matrix[, species_order]
  
  results <- data.frame(
    CpG_ID = rownames(slope_matrix),
    Lambda = NA_real_,
    P_value = NA_real_,
    LogLik = NA_real_,
    Sigma2 = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(slope_matrix))) {
    slope_vec <- as.numeric(slope_matrix[i, ])
    se_vec <- as.numeric(se_matrix[i, ])
    names(slope_vec) <- species_order
    names(se_vec) <- species_order
    
    cpg_id <- rownames(slope_matrix)[i]
    message(sprintf("Processing CpG: %s (%d of %d)", cpg_id, i, nrow(slope_matrix)))
    
    res <- tryCatch({
      out <- phylosig(tree, x = slope_vec, method = "lambda", test = TRUE, se = se_vec)
      
      # Show structure of output for inspection
      # Uncomment this to inspect the output
      # print(str(out))
      
      c(
        Lambda  = out$lambda,
        P_value = out$P,
        LogLik  = out$logL,
        Sigma2  = out$sig2
      )
    }, error = function(e) {
      message(sprintf("Error at CpG %s: %s", cpg_id, conditionMessage(e)))
      return(c(Lambda = NA, P_value = NA, LogLik = NA, Sigma2 = NA))
    })
    
    # Check if result has correct length before assigning
    if (length(res) == 4) {
      results[i, 2:5] <- res
    } else {
      message(sprintf("Warning: Unexpected result length at CpG %s", cpg_id))
    }
  }
  
  return(results)
}


skin_phylo_lambda_results <- lambda_with_se(
  tree = skin_species_tree,
  slope_matrix = slope_matrix,
  se_matrix = se_matrix
)

saveRDS(skin_phylo_lamda_result, file = "~/Phd_data/outputs/skin_phylo_lamda_result.rds")

