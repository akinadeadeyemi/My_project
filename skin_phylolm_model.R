######### Author : Adeyemi Akinade
######### Place : Clemson University, SC
######### lab : Gopalan's lab
######### Date : March 19, 2025
######### Project : Comparative Epigenetics of Ageing
######### This script was written to determine the linear regression between age and beta values of DNA methylation sites with phylogenetic frame work


######################## PHYLOGENETIC LINEAR REGRESSION IN DNA METHYLATION SKIN SAMPLES ONLY ##########################################
# Incorporating phylogeny into the model ----------------------------------
### Loading the necessary libraries
library(dplyr)
library(phylolm)
library(ape)

## working on the loading datasets
skin_species_tree <- read.tree("~/Phd_data/skin_species_tree.newick")
skin_metadata_wt_life_span <- readRDS("~/Phd_data/skin_metadata_wt_life_span.rds")
skin_metadata_wt_life_span_betas <- readRDS("~/Phd_data/skin_metadata_wt_life_span_betas.rds")

## processing the tree data
# Remove the OTL suffix to match metadata species
clean_labels <- gsub("_ott\\d+", "", skin_species_tree$tip.label)

# Replace underscores with spaces (if your metadata uses spaces)
clean_labels <- gsub("_", " ", clean_labels)

# Assign cleaned labels back to the tree
skin_species_tree$tip.label <- clean_labels
skin_species_tree$tip.label 


skin_regression_phylogenetic_simple <- function(data, metadata, phy_tree) {
  
  cat("Starting Phylogenetic Linear Models (PGLM) across species using 'age_years'...\n")
  
  # Matching species between metadata and the tree
  common_species <- intersect(metadata$species, phy_tree$tip.label)
  
  if (length(common_species) < 3) {
    stop("Not enough common species (need at least 3) between your data and the phylogenetic tree for PGLM. Please check your inputs.")
  }
  
  # Filtering the metadata and pruning the tree to include only these common species.
  filtered_metadata <- metadata[metadata$species %in% common_species, ]
  pruned_tree <- ape::keep.tip(phy_tree, common_species)
  
  # Ensure the tree tip labels are correctly ordered for later matching
  # (though direct matching is handled by row names later)
  # For robustness, let's also ensure 'filtered_metadata' species are ordered consistently if needed,
  # but the aggregation step below will handle it based on unique species names.
  
  # 2. Aggregate 'age_years' data to the species level
  # Since you have multiple individuals per species, you need to decide how to get one 'age_years' per species.
  # Common choices: mean age, max age. Let's use MEAN age for this example.
  # If you want max age, change 'mean' to 'max'.
  
  species_age_values <- filtered_metadata %>%
    group_by(species) %>%
    summarise(
      # IMPORTANT: Choose your aggregation method here.
      # Using mean age for the species
      species_age_years = mean(as.numeric(as.character(age_years)), na.rm = TRUE)
      # Or, if you want max age for the species (e.g., as a proxy for lifespan for your samples):
      # species_age_years = max(as.numeric(as.character(age_years)), na.rm = TRUE)
    ) %>%
    ungroup() %>%
    # Convert to a named vector, matching the pruned_tree tip labels order
    tibble::deframe() # Converts a two-column tibble to a named vector
  
  # Ensure the species_age_values vector is ordered according to the pruned_tree$tip.label
  # This is crucial for matching with the methylation data later if not using direct named vector matching.
  # However, deframe() creates a named vector, which is robust to order as long as names match.
  # Let's explicitly order it to be safe for phylolm's internal checks.
  species_age_values <- species_age_values[pruned_tree$tip.label]
  
  
  # 3. Aggregate methylation data per species (this part is largely the same as before)
  species_methylation_list <- list()
  for (s_name in common_species) {
    s_accessions <- filtered_metadata$geo_accession[filtered_metadata$species == s_name]
    s_accessions_in_data <- s_accessions[s_accessions %in% colnames(data)]
    
    if (length(s_accessions_in_data) > 0) {
      # Aggregate (mean) methylation for all individuals of this species for each CpG site
      species_methylation_list[[s_name]] <- rowMeans(data[, s_accessions_in_data, drop = FALSE], na.rm = TRUE)
    } else {
      species_methylation_list[[s_name]] <- rep(NA, nrow(data))
    }
  }
  species_methylation_data <- do.call(cbind, species_methylation_list)
  rownames(species_methylation_data) <- rownames(data)
  species_methylation_data <- species_methylation_data[, pruned_tree$tip.label, drop = FALSE] # Ensure order matches tree
  
  cpg_sites <- rownames(species_methylation_data)
  
  # 4. Set up a place to store our PGLM results
  pglm_results <- data.frame(
    cpg_site = character(),
    intercept = numeric(),
    slope = numeric(),
    p_value = numeric(),
    std_error = numeric(),
    conf_int_low = numeric(),
    conf_int_upp = numeric(),
    lambda = numeric()
  )
  
  # 5. Run PGLM for each CpG site
  for (i in seq_len(nrow(species_methylation_data))) {
    cpg_site_name <- cpg_sites[i]
    cat("  Running PGLM for CpG site:", cpg_site_name, "\n")
    
    # Create a temporary data frame for this specific CpG site's regression
    temp_data_for_pglm <- data.frame(
      # NOW using the species-aggregated 'age_years'
      age_years = species_age_values,
      # And the species-aggregated 'methylation' for this CpG site
      methylation = as.numeric(species_methylation_data[i, ])
    )
    # The row names are already set by the aggregation and explicit ordering
    rownames(temp_data_for_pglm) <- pruned_tree$tip.label
    
    
    # Remove any rows with missing data (NA) for a clean regression
    complete_cases_pglm <- complete.cases(temp_data_for_pglm)
    temp_data_clean_pglm <- temp_data_for_pglm[complete_cases_pglm, , drop = FALSE]
    
    # Prune the tree again to only include species with complete data for this CpG site
    current_pruned_tree <- ape::keep.tip(pruned_tree, rownames(temp_data_clean_pglm))
    
    if (nrow(temp_data_clean_pglm) < 3) {
      warning(paste("Not enough complete data points for CpG site:", cpg_site_name, " for PGLM. Skipping."))
      next
    }
    
    tryCatch({
      pglm_model <- phylolm(age_years ~ methylation, data = temp_data_clean_pglm, phy = current_pruned_tree, model = "lambda")
      pglm_summary <- summary(pglm_model)
      
      conf_int <- confint(pglm_model, level = 0.95)
      
      pglm_results <- dplyr::bind_rows(pglm_results,
                                       data.frame(
                                         cpg_site = cpg_site_name,
                                         intercept = coef(pglm_model)[1],
                                         slope = coef(pglm_model)[2],
                                         p_value = coef(pglm_summary)[2, 4],
                                         std_error = coef(pglm_summary)[2, 2],
                                         conf_int_low = conf_int["methylation", 1],
                                         conf_int_upp = conf_int["methylation", 2],
                                         lambda = pglm_model$optpar
                                       )
      )
    }, error = function(e) {
      warning(paste("Error performing PGLM for CpG site:", cpg_site_name, "-", e$message, "Results will be NA."))
      pglm_results <- dplyr::bind_rows(pglm_results,
                                       data.frame(
                                         cpg_site = cpg_site_name,
                                         intercept = NA, slope = NA, p_value = NA, std_error = NA,
                                         conf_int_low = NA, conf_int_upp = NA, lambda = NA
                                       )
      )
    })
  }
  
  return(pglm_results)
}



skin_phylolm_result <- skin_regression_phylogenetic_simple(data = skin_metadata_wt_life_span_betas,
                                                             metadata = skin_metadata_wt_life_span,
                                                             phy_tree = skin_species_tree)
saveRDS(skin_phylolm_result, file = "~/Phd_data/outputs/skin_phylolm_result.rds")


