######### Author : Adeyemi Akinade
######### Place : Clemson University, SC
######### lab : Gopalan's lab
######### Date : March 19, 2025
######### Project : Comparative Epigenetics of Ageing
######### This script was written to determine the linear regression between age and beta values of DNA methylation sites in tissues


######################## LINEAR REGRESSION IN DNA METHYLATION BLOOD SAMPLES ONLY ##########################################

## saving the data as .rds
cat("loading the DNA methylation data...", "\n")
blood_metadata_wt_life_span_betas <- readRDS("~/Phd_data/blood_metadata_wt_life_span_betas.rds")
cat("Finished loading the DNA mthylation data...", "\n")

cat("loading the species metadata...", "\n")
blood_metadata_wt_life_span <- readRDS("~/Phd_data/blood_metadata_wt_life_span.rds")
cat("Finished loading the species metadata...", "\n")



########## conducting the linear regression in age and beta values across the cpg sites ############
library(dplyr)
blood_regression_test_wt_max_lifespan <- function(data, metadata) {
  # Get unique species from the metadata
  unique_species <- unique(metadata$species)
  
  # Initialize a list to store results for each species
  species_results <- list()
  
  for (sp in unique_species) {
    cat("Processing species:", sp, "\n")
    
    # Subset metadata for the current species
    species_metadata <- metadata[metadata$species == sp, ]
    if (nrow(species_metadata) == 0) {
      warning(paste("No metadata found for species:", sp))
      next
    }
    
    # Extract sample IDs and corresponding age data
    sampled_ids <- species_metadata$geo_accession
    matched_indices <- match(sampled_ids, colnames(data))
    
    if (all(is.na(matched_indices))) {
      warning(paste("No matching samples found for species:", sp))
      next
    }
    
    # Subset data
    sampled_data <- data[, matched_indices, drop = FALSE]
    
    if (nrow(sampled_data) == 0) {
      warning(paste("No data available for species:", sp))
      next
    }
    
    # Extract lifespan data
    max_lifespan_data <- as.numeric(as.character(species_metadata$max_lifespan))
    
    # Ensure cpg_site is correctly assigned
    cpg_sites <- if (!is.null(rownames(sampled_data)) && length(rownames(sampled_data)) > 0) {
      rownames(sampled_data)
    } else {
      seq_len(nrow(sampled_data))
    }
    
    # Initialize dataframe
    species_out <- data.frame(
      species = rep(sp, length(cpg_sites)),
      cpg_site = cpg_sites,
      intercept = NA,
      slope = NA,
      p_value = NA,
      std_error = NA,
      conf_int_low = NA,
      conf_int_upp = NA,
      residuals = I(vector("list", length(cpg_sites))) # Initialize list column for residuals
    )
    
    # Perform regression for each CpG site
    for (i in seq_len(nrow(sampled_data))) {
      x <- as.numeric(sampled_data[i, ])
      y <- max_lifespan_data
      
      # Perform linear regression
      model <- lm(y ~ x)
      model_summary <- summary(model)
      
      # Extract regression results
      species_out$intercept[i] <- coef(model)[1]
      species_out$slope[i] <- coef(model)[2]
      species_out$p_value[i] <- coef(summary(model))[2, 4]  # p-value of the slope
      species_out$std_error[i] <- coef(summary(model))[2, 2]  # Standard error
      
      # Confidence intervals
      conf_int <- confint(model, level = 0.95)
      species_out$conf_int_low[i] <- conf_int[2, 1]
      species_out$conf_int_upp[i] <- conf_int[2, 2]
      
      # Store residuals
      species_out$residuals[[i]] <- residuals(model)
    }
    
    # Store results in the list
    species_results[[sp]] <- species_out
  }
  
  # Combine results into a single dataframe
  combined_results <- dplyr::bind_rows(species_results)
  
  return(combined_results)
}

### conducting the correlation analysis
blood_regression_wt_max_lifespan_result <- blood_regression_test_wt_max_lifespan(data = blood_metadata_wt_life_span_betas,
                                                                       metadata = blood_metadata_wt_life_span)

saveRDS(blood_regression_wt_max_lifespan_result, file = "~/Phd_data/outputs/blood_regression_wt_max_lifespan_result.rds")

