######### Author : Adeyemi Akinade
######### Place : Clemson University, SC
######### lab : Gopalan's lab
######### Date : June 16, 2025
######### Project : Comparative Epigenetics of Ageing
######### This script was written to determine the linear regression between age and beta values of DNA methylation sites in tissues


######################## PHYLOGENETIC REGRESSION MODEL IN BLOOD SAMPLES ONLY ##########################################

library(dplyr)
library(phylolm)
library(ape)

## working on the loading skin species tree
blood_species_tree <- read.tree("~/Phd_data/blood_species.nwk")
blood_species_tree$tip.label
plot(blood_species_tree, cex = 0.8)
### Loading the datasets
blood_regression_wt_lifespan_result <- readRDS("~/Phd_data/outputs/blood_regression_test_wt_lifespan_result.rds")

### Extract the betip.label### Extract the beta values  (DNAm rates) from the age vs DNAm regression
blood_regression_slope <- blood_regression_wt_lifespan_result [ , c("species", "cpg_site", "slope")]
blood_regression_std_error <- blood_regression_wt_lifespan_result [ , c("species", "cpg_site", "std_error")]

### Reshape the slope data 
blood_regression_slope_wide <- blood_regression_slope %>%
  tidyr::pivot_wider(
    names_from = species,
    values_from = slope,
    id_cols = cpg_site 
  ) %>%
  tibble::column_to_rownames(var = "cpg_site")

View(blood_regression_slope_wide)

### Reshape the slope data 
blood_regression_std_error_wide <- blood_regression_std_error %>%
  tidyr::pivot_wider(
    names_from = species,
    values_from = std_error,
    id_cols = cpg_site 
  ) %>%
  tibble::column_to_rownames(var = "cpg_site")

View(blood_regression_std_error_wide)


library(dplyr)
colnames(blood_regression_std_error_wide)
### adding undersore to fit the phylogenetic tree style
colnames(blood_regression_slope_wide) <- gsub(" ", "_", colnames(blood_regression_slope_wide))
colnames(blood_regression_std_error_wide) <- gsub(" ", "_", colnames(blood_regression_std_error_wide))


### Working on the phylosignal


##Setting the difference between the tree_species and DNAm_rates_with_age_wide datasets.
# setdiff(blood_species_tree$tip.label, colnames(blood_regression_slope_wide))
# setdiff(blood_species_tree$tip.label, colnames(blood_regression_std_error_wide))

## correcting for the unmatched species in tree_species and DNAm_rates_with_age_wide
#colnames(DNAm_rates_with_age_wide)[colnames(DNAm_rates_with_age_wide) == "Cephalorhynchus_hectori_hectori"] <- "Cephalorhynchus_hectori"
#colnames(DNAm_rates_std_err_wide)[colnames(DNAm_rates_std_err_wide) == "Cephalorhynchus_hectori_hectori"] <- "Cephalorhynchus_hectori"

## Write file as RDS
# saveRDS(blood_regression_slope_wide, file = "blood_regression_slope_wide.rds")
# saveRDS(blood_regression_std_error_wide, file = "blood_regression_std_error_wide.rds")


# Ensure species order matches tree
slope_matrix <- blood_regression_slope_wide[, blood_species_tree$tip.label]
se_matrix <- blood_regression_std_error_wide[, blood_species_tree$tip.label]

colnames(se_matrix)
blood_species_tree$tip.label

## load phytools
library(phytools)

# Ensure species order matches tree
species_order <- blood_species_tree$tip.label

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


blood_phylo_lambda_results <- lambda_with_se(
  tree = blood_species_tree,
  slope_matrix = slope_matrix,
  se_matrix = se_matrix
)

saveRDS(blood_phylo_lamda_result, file = "~/Phd_data/outputs/blood_phylo_lamda_result.rds")

