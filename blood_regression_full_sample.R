######### Author : Adeyemi Akinade
######### Place : Clemson University, SC
######### lab : Gopalan's lab
######### Date : March 19, 2025
######### Project : Comparative Epigenetics of Ageing
######### This script was written to determine the linear regression between age and beta values of DNA methylation sites in tissues


######################## LINEAR REGRESSION IN DNA METHYLATION BLOOD SAMPLES ONLY ##########################################
library(readr)
blood_skin_metadata_with_order <- read_csv("blood_skin_metadata_with_order.csv")
View(blood_skin_metadata_with_order)

## Extracting the blood metadata
blood_metadata <- blood_skin_metadata_with_order %>% filter(tissue == "Blood") 
blood_metadata <- blood_metadata %>% filter(!age_years =="NA")

# Filter species where all tissue types have more than 10 samples
filtered_species <- blood_metadata %>%
  group_by(species) %>%        # Group by species and tissue type
  tally() %>%                   # Count the number of rows in each group
  filter(n >= 10) %>%                  # Keep species where count > 10
  ungroup()                           # Remove grouping after filtering

View(filtered_species)
unique(filtered_species$species)


# Extract the species names that meet the condition
species_with_10 <- filtered_species$species

# Subset the metadata to include only those species with blood and skin count > 10
blood_metadata_final <- blood_metadata %>%
  filter(species %in% species_with_10)
View (blood_metadata_final)

#selecting unique species in skin_metadata_final
final_blood_species <- data.frame(unique(blood_metadata_final$species)) 
colnames (final_blood_species) <- "species"
View(final_blood_species)

# loading the file
library(readxl)
species_metadata_with_max_lifespan <- read_excel("species_metadata_with_max_lifespan.xlsx")
View(species_metadata_with_max_lifespan)

## changing column names
library(dplyr)
species_metadata_with_max_lifespan <- species_metadata_with_max_lifespan %>% rename(max_lifespan = 'Maximum Lifespan (yrs)')
species_metadata_with_max_lifespan <- species_metadata_with_max_lifespan %>% rename(species = 'Species Latin Name')


# merging species that are common in skin_metadata_final and species_with_lifespan data to add lifespan column.
species_with_lifespan_blood <- merge(species_metadata_with_max_lifespan, final_blood_species, by = "species")
(species_with_lifespan_blood)

# merging the species with life span with the skin_metadata
blood_metadata_wt_life_span <- merge(blood_metadata_final, species_with_lifespan_blood, by = "species")
View(blood_metadata_wt_life_span)

### calculating the current age of individuals as a % of max lifespan of the species.
blood_metadata_wt_life_span$lifespan_percent <- (blood_metadata_wt_life_span$age_years/blood_metadata_wt_life_span$max_lifespan) * 100


## Extracting the sample ids of skin_metadata_wt_lifespan from the DNA methylation data
blood_metadata_wt_life_span_betas <- merged_betas_in_sample[, colnames(merged_betas_in_sample) %in% blood_metadata_wt_life_span$geo_accession]
blood_metadata_wt_life_span <- blood_metadata_wt_life_span[blood_metadata_wt_life_span$geo_accession %in% colnames(blood_metadata_wt_life_span_betas), ]

## saving the data as .rds
saveRDS(blood_metadata_wt_life_span_betas, file = "~/Phd_data/blood_metadata_wt_life_span_betas.rds")
saveRDS(blood_metadata_wt_life_span, file = "~/Phd_data/blood_metadata_wt_life_span.rds")


View(blood_metadata_wt_life_span_betas)

########## conducting the linear regression in age and beta values across the cpg sites ############
library(dplyr)
blood_regression_test <- function(data, metadata) {
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
    
    # Extract age data
    age_data <- as.numeric(as.character(species_metadata$age_years))
    
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
      y <- age_data
      
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
blood_regression_result <- blood_regression_test(data = blood_metadata_wt_life_span_betas,
                                                                         metadata = blood_metadata_wt_life_span)

saveRDS(blood_regression_result, file = "~/Phd_data/outputs/blood_regression_result.rds")

View(blood_regression_result)
