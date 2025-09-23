######### Author : Adeyemi Akinade
######### Place : Clemson University, SC
######### lab : Gopalan's lab
######### Date : March 19, 2025
######### Project : Comparative Epigenetics of Ageing



# This script is meant to carry out correlation test from blood samples. The codes will select samples that have at least 10 samples 
# in each species present in the skin_metadata1. datasets and in the DNA methylation datasets (named: merged_betas.RDS).

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


##### conducting correlation analysis by creating a function
library(dplyr)

blood_correlation_test <- function(data, metadata) {
  # Get unique species from the metadata
  unique_species <- unique(metadata$species)
  
  # Initialize a list to store results for each species
  species_results <- list()
  
  for (sp in unique_species) {
    cat("Processing species:", sp, "\n")
    
    # Subset metadata for the current species
    species_metadata <- metadata[metadata$species == sp, ]
    
    # Use all available samples for this species
    n_samples <- nrow(species_metadata)
    
    # Initialize dataframe for storing correlations and p-values
    n_features <- nrow(data)
    
    species_out <- data.frame(
      species = sp,  # Add species column
      cpg_site = rownames(data), 
      correlation = NA, 
      p_value = NA, 
      conf_int_low = NA,
      conf_int_upp = NA,
      std_err = NA
    )
    
    # Extract sample IDs and corresponding age data
    sampled_ids <- species_metadata$geo_accession
    age_data <- as.numeric(as.character(species_metadata$age_years))  # Ensure numeric
    
    # Subset the data for the sampled IDs
    matched_indices <- match(sampled_ids, colnames(data))
    
    if (any(is.na(matched_indices))) {
      warning("Some sampled IDs do not match column names in the data.")
    }
    
    sampled_data <- data[, matched_indices, drop = FALSE]
    
    # Calculate correlation and p-values for each row (feature)
    for (i in seq_len(n_features)) {
      x <- as.numeric(sampled_data[i, ])
      y <- age_data
      
      # Perform the correlation test
      cor_test_result <- cor.test(x, y, use = "complete.obs")
      
      species_out$correlation[i] <- cor_test_result$estimate
      species_out$p_value[i] <- cor_test_result$p.value
      species_out$conf_int_low[i] <- cor_test_result$conf.int[1]
      species_out$conf_int_upp[i] <- cor_test_result$conf.int[2]
      species_out$std_err[i] <- sqrt((1 - cor_test_result$estimate^2) / (n_samples - 2))
    }
    
    # Store the results for the current species in the list
    species_results[[sp]] <- species_out
  }
  
  # Combine all species results into one dataframe
  blood_combined_results <- bind_rows(species_results)
  
  return(blood_combined_results)
}

### conducting the correlation analysis
blood_correlation_result <- blood_correlation_test(data = blood_metadata_wt_life_span_betas ,
                      metadata = blood_metadata_wt_life_span )

saveRDS(blood_correlation_result, file = "~/Phd_data/outputs/blood_correlation_full_sample.rds")




##### conducting correlation analysis using DNA methylation site and percentage (%) lifespan  by creating a function
library(dplyr)

blood_correlation_wt_lifespan_test <- function(data, metadata) {
  # Get unique species from the metadata
  unique_species <- unique(metadata$species)
  
  # Initialize a list to store results for each species
  species_results <- list()
  
  for (sp in unique_species) {
    cat("Processing species:", sp, "\n")
    
    # Subset metadata for the current species
    species_metadata <- metadata[metadata$species == sp, ]
    
    # Use all available samples for this species
    n_samples <- nrow(species_metadata)
    
    # Initialize dataframe for storing correlations and p-values
    n_features <- nrow(data)
    
    species_out <- data.frame(
      species = sp,  # Add species column
      cpg_site = rownames(data), 
      correlation = NA, 
      p_value = NA, 
      conf_int_low = NA,
      conf_int_upp = NA,
      std_err = NA
    )
    
    # Extract sample IDs and corresponding age data
    sampled_ids <- species_metadata$geo_accession
    lifespan_data <- as.numeric(as.character(species_metadata$lifespan_percent))  # Ensure numeric
    
    # Subset the data for the sampled IDs
    matched_indices <- match(sampled_ids, colnames(data))
    
    if (any(is.na(matched_indices))) {
      warning("Some sampled IDs do not match column names in the data.")
    }
    
    sampled_data <- data[, matched_indices, drop = FALSE]
    
    # Calculate correlation and p-values for each row (feature)
    for (i in seq_len(n_features)) {
      x <- as.numeric(sampled_data[i, ])
      y <- lifespan_data
      
      # Perform the correlation test
      cor_test_result <- cor.test(x, y, use = "complete.obs")
      
      species_out$correlation[i] <- cor_test_result$estimate
      species_out$p_value[i] <- cor_test_result$p.value
      species_out$conf_int_low[i] <- cor_test_result$conf.int[1]
      species_out$conf_int_upp[i] <- cor_test_result$conf.int[2]
      species_out$std_err[i] <- sqrt((1 - cor_test_result$estimate^2) / (n_samples - 2))
    }
    
    # Store the results for the current species in the list
    species_results[[sp]] <- species_out
  }
  
  # Combine all species results into one dataframe
  blood_combined_results <- bind_rows(species_results)
  
  return(blood_combined_results)
}

### conducting the correlation analysis
blood_correlation_wt_lifespan_result <- blood_correlation_wt_lifespan_test(data = blood_metadata_wt_life_span_betas ,
                                                   metadata = blood_metadata_wt_life_span )

saveRDS(blood_correlation_wt_lifespan_result, file = "~/Phd_data/outputs/blood_correlation_wt_lifespan_full_sample.rds")




