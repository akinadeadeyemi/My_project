######### Author : Adeyemi Akinade
######### Place : Clemson University, SC
######### lab : Gopalan's lab
######### Date : September 25, 2025
######### Project : Comparative Epigenetics of Ageing
######### This script was written to determine the phylogenetic signal of the DNAm rates at CpG sites in blood tissue sample
######### To do this, the slope of linear model between DNAm level and relative age was used


######################## PHYLOGENETIC REGRESSION MODEL IN BLOOD SAMPLES ONLY ##########################################
library(dplyr)
library(phytools)
library(ape)

## working on the loading skin species tree
blood_tree <- read.tree("~/Phd_data/blood_species_34.nwk")
# blood_tree$tip.label
# plot(blood_tree, cex = 0.8, edge.color = "blue")

### Loading the datasets
# blood_regression_wt_relative_age <- readRDS("~/Phd_data/outputs/blood_regression_wt_relative_age.rds")

### Extract the betip.label### Extract the beta values  (DNAm rates) from the age vs DNAm regression
# blood_regression_wt_relative_age_slope <- blood_regression_wt_relative_age [ , c("species", "cpg_site", "slope")]
# blood_regression_wt_relative_age_std_error <- blood_regression_wt_relative_age [ , c("species", "cpg_site", "std_error")]

### Reshape the slope data 
#blood_regression_wt_relative_age_slope_wide <- blood_regression_wt_relative_age_slope %>%
#  tidyr::pivot_wider(
 #   names_from = species,
  #  values_from = slope,
   # id_cols = cpg_site 
 # ) %>%
#  tibble::column_to_rownames(var = "cpg_site")

# View(blood_regression_wt_relative_age_slope_wide)

### Reshape the slope data 
# blood_regression_wt_relative_age_std_error_wide <- blood_regression_wt_relative_age_std_error %>%
#  tidyr::pivot_wider(
 #   names_from = species,
  #  values_from = std_error,
#    id_cols = cpg_site 
#  ) %>%
 # tibble::column_to_rownames(var = "cpg_site")
# View(blood_regression_wt_relative_age_std_error_wide)


### adding undersore to fit the phylogenetic tree style
# colnames(blood_regression_wt_relative_age_slope_wide) <- gsub(" ", "_", colnames(blood_regression_wt_relative_age_slope_wide))
# colnames(blood_regression_wt_relative_age_std_error_wide) <- gsub(" ", "_", colnames(blood_regression_wt_relative_age_std_error_wide))

#### set  a difference between tree names and data column names
# setdiff( colnames(blood_regression_wt_relative_age_slope_wide), blood_tree$tip.label)
# setdiff(blood_tree$tip.label, colnames(blood_regression_wt_relative_age_std_error_wide))


##### Changing the names on the slope and standard error dataset to match the 
# colnames(blood_regression_wt_relative_age_slope_wide)[colnames(blood_regression_wt_relative_age_slope_wide) == "us_scrofa_domesticus"] <- "Sus_scrofa_domesticus"
# colnames(blood_regression_wt_relative_age_std_error_wide)[colnames(blood_regression_wt_relative_age_std_error_wide) == "Canis_lupus_familiaris"] <- "Canis_lupus"

# saveRDS(blood_regression_wt_relative_age_slope_wide, file = "~/Phd_data/outputs/blood_regression_wt_relative_age_slope_wide.rds")
# saveRDS(blood_regression_wt_relative_age_std_error_wide, file = "~/Phd_data/outputs/blood_regression_wt_relative_age_std_error_wide.rds")



################# RUNNING PHYLOSIGNAL ANALYSIS ########################################
##                                                                                  ###
##                                                                                  ###
#######################################################################################
## loadiing the datasets
blood_regression_wt_relative_age_slope_wide <- readRDS("~/Phd_data/outputs/blood_regression_wt_relative_age_slope_wide.rds")
slope_dat <- blood_regression_wt_relative_age_slope_wide

blood_regression_wt_relative_age_std_error_wide <- readRDS("~/Phd_data/outputs/blood_regression_wt_relative_age_std_error_wide.rds")
se_dat <- blood_regression_wt_relative_age_std_error_wide


#### running phylosig

library(parallel)
run_phylosig <- function(i, tree, slope_dat, se_dat) {
  if (i %% 500 == 0) cat("Processing CpG site", i, "of", nrow(slope_dat), "\n")
  
  slopes <- as.numeric(slope_dat[i, ])
  std_errs <- as.numeric(se_dat[i, ])
  
  names(slopes) <- colnames(slope_dat)
  names(std_errs) <- colnames(se_dat)
  
  slopes <- slopes[tree$tip.label]
  std_errs <- std_errs[tree$tip.label]
  
  # Try-catch block
  out <- tryCatch({
    x <- phylosig(tree, x = slopes, se = std_errs, method = "lambda", test = TRUE)
    c(lambda = x$lambda,
      sig2   = x$sig2,
      p_val  = if (!is.null(x$P)) x$P else NA_real_)
  }, error = function(e) {
    c(lambda = NA_real_, sig2 = NA_real_, p_val = NA_real_)
  })
  
  return(out)
}


# detect available cores on your HPC node
ncores <- detectCores() - 1  

# run in parallel
blood_phylo_relative_age_results <- mclapply(1:nrow(slope_dat), run_phylosig,
                                             tree = blood_tree,
                                             slope_dat = slope_dat,
                                             se_dat = se_dat,
                                             mc.cores = ncores)

# convert to data.frame
blood_phylo_relative_age_results <- as.data.frame(do.call(rbind, blood_phylo_relative_age_results))
rownames(blood_phylo_relative_age_results) <- rownames(slope_dat)

saveRDS(blood_phylo_relative_age_results, "~/Phd_data/outputs/blood_regression_relative_age_phylosig_result.rds")

