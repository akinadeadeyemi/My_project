######### Author : Adeyemi Akinade
######### Place : Clemson University, SC
######### lab : Gopalan's lab
######### Date : September 22, 2025
######### Project : Comparative Epigenetics of Ageing
######### This script was written to determine the linear regression between age and beta values of 
######### DNAm levels in skin tissue in each species


######################## LINEAR REGRESSION IN DNA METHYLATION SKIN SAMPLES (AGE IN YEARS) ##########################################

## saving the data as .rds
cat("loading the DNA methylation data...", "\n")
skin_metadata_wt_life_span_betas <- readRDS("~/Phd_data/skin_metadata_wt_life_span_betas.rds")
cat("Finished loading the DNA mthylation data...", "\n")

cat("loading the species metadata...", "\n")
skin_metadata_wt_life_span <- readRDS("~/Phd_data/skin_metadata_wt_life_span.rds")
cat("Finished loading the species metadata...", "\n")

# Loading the necessary libraries
library(data.table)
library(dplyr)

skin_regression_stream_residuals <- function(data, metadata, output_file = "skin_regression_age_in_years_wt_residuals_ID.csv") {
  library(data.table)
  
  # Initialize CSV
  fwrite(
    data.table(
      species = character(),
      cpg_site = character(),
      sample_id = character(),
      residual = numeric(),
      intercept = numeric(),
      slope = numeric(),
      p_value = numeric(),
      std_error = numeric(),
      conf_low = numeric(),
      conf_upp = numeric()
    ),
    file = output_file
  )
  
  for (sp in unique(metadata$species)) {
    cat("Processing species:", sp, "\n")
    species_metadata <- metadata[metadata$species == sp, ]
    sampled_ids <- species_metadata$geo_accession
    matched_indices <- match(sampled_ids, colnames(data))
    sampled_data <- data[, matched_indices, drop = FALSE]
    age_years_data <- as.numeric(as.character(species_metadata$age_years))
    cpg_sites <- rownames(sampled_data)
    
    for (i in seq_len(nrow(sampled_data))) {
      model <- lm(age_years_data ~ as.numeric(sampled_data[i, ]))
      coefs <- coef(model)
      summ <- summary(model)
      conf <- confint(model)
      
      residuals_dt <- data.table(
        species = sp,
        cpg_site = cpg_sites[i],
        sample_id = sampled_ids,
        residual = residuals(model),
        intercept = coefs[1],
        slope = coefs[2],
        p_value = summ$coefficients[2, 4],
        std_error = summ$coefficients[2, 2],
        conf_low = conf[2, 1],
        conf_upp = conf[2, 2]
      )
      
      fwrite(residuals_dt, file = output_file, append = TRUE)
    }
  }
  
  cat("Finished processing all species. Residuals saved to:", output_file, "\n")
}


### calling the function
skin_regression_stream_residuals(
  data = skin_metadata_wt_life_span_betas,
  metadata = skin_metadata_wt_life_span,
  output_file = "~/Phd_data/outputs/skin_regression_age_in_years_wt_residuals_ID.csv"
)



