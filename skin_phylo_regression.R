######### Author : Adeyemi Akinade
######### Place : Clemson University, SC
######### lab : Gopalan's lab
######### Date : June 16, 2025
######### Project : Comparative Epigenetics of Ageing
######### This script was written to determine the linear regression between age and beta values of DNA methylation sites in tissues


######################## PHYLOGENETIC REGRESSION MODEL IN SKIN SAMPLES ONLY ##########################################
cat("Loading the required libraries...")
library(dplyr)
library(phylolm)
library(ape)

## working on the loading skin species tree
skin_species_tree <- read.tree("~/Phd_data/skin_species.nwk")
plot(skin_species_tree, cex =0.7)
### Loading other datasets
skin_metadata_wt_life_span <- readRDS("~/Phd_data/skin_metadata_wt_life_span.rds")
skin_metadata_wt_life_span_betas <- readRDS("~/Phd_data/skin_metadata_wt_life_span_betas.rds")
skin_regression_wt_lifespan_result <- readRDS("~/Phd_data/outputs/skin_regression_wt_lifespan_result.rds")

### selecting for the significant sites in CpG sites


### Loading the matrix of slope and standard error
skin_regression_slope <- readRDS("~/Phd_data/skin_regression_slopes.rds")
View(skin_regression_slope)
skin_regression_std_error <- readRDS("~/Phd_data/skin_regression_std_error.rds")

### Ensure species order matches tree
slope_matrix <- skin_regression_slope[, skin_species_tree$tip.label]
se_matrix <- skin_regression_std_error[, skin_species_tree$tip.label]

colnames(se_matrix)
skin_species_tree$tip.label
setdiff(skin_species_tree$tip.label, colnames(skin_regression_slope))


## load phytools
library(phytools)

# Ensure species order matches tree
species_order <- skin_species_tree$tip.label

cat("Starting the phylogenetics analysis...")
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

#### starting the analysis
skin_phylo_lambda_results <- lambda_with_se(
  tree = skin_species_tree,
  slope_matrix = skin_regression_slope,
  se_matrix = skin_regression_std_error
)


saveRDS(skin_phylo_lamda_result, file = "~/Phd_data/outputs/skin_phylo_lamda_result.rds")
cat("Finished with the analysis and storage of the output...")



### trying a small example

skin_rates_50 <- skin_regression_slope[1:50, ]
saveRDS(skin_rates_50, file = "~/Phd_data/skin_rates_50.rds")
skin_std_err_50 <- skin_regression_std_error[1:50, ]
saveRDS(skin_std_err_50, file = "~/Phd_data/skin_std_err_50.rds")
### match species in tree and slope/std_error data
### Ensure species order matches tree
slope_matrix <- skin_rates_50[, skin_species_tree$tip.label]
se_matrix <- skin_std_err_50[, skin_species_tree$tip.label]


### runing the function
skin_phylo_test <- lambda_with_se(
  tree = skin_species_tree,
  slope_matrix = skin_rates_50,
  se_matrix = skin_std_err_50
)

