########## Author: Adeyemi Akinade
########## Project: Comparative Epigenetics of Aging in mammalian species
########## Date: 22 February, 2025
########## Place: Clemson University, SC, USA.

## loading libraries
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)

# selecting the blood and skin samples
# blood_skin_metadata_with_order <- read_csv("blood_skin_metadata_with_order.csv")

# selecting for only blood samples  
# skin_metadata <- blood_skin_metadata_with_order %>% filter(tissue == "Skin")

# remove all NAs in the age column
# blood_metadata <- blood_metadata %>% filter(!age_years == "NA") 

# load the data
# write.csv(skin_metadata1, file = "skin_metadata1.csv", row.names = FALSE)
# skin_metadata1 <- read_csv("~/Phd_data/skin_metadata1.csv")

### downsampling to N = 30 for each species
### Subsample N in the meta data
# set.seed(123)
# sample_species <- function(metadata, n_samples) {
  #Get unique species from the metadata
 # unique_species <- unique(metadata$species)

# Create an empty list to store sampled data
#  sampled_list <- list()

 # for (species in unique_species) {
#   cat("Processing species:", species, "\n")

 #Subset metadata for the current species
#  species_metadata <- metadata[metadata$species == species, ]

# Check if the species has enough rows to sample
# if (nrow(species_metadata) >= n_samples) {
 #Sample exactly `n_samples` rows
#  sampled_data <- species_metadata[sample(1:nrow(species_metadata), 
 #                                        size = n_samples, 
  #                                      replace = FALSE), ]

# Store the result in the list
#   sampled_list[[species]] <- sampled_data
#  } else {
#   cat("Skipping species:", species, " - not enough samples\n")
#  }
# }

# Combine all sampled data into one dataframe
#final_sampled_data <- do.call(rbind, sampled_list)

#return(final_sampled_data)
#}

#skin_downsampled_N30 <- sample_species(metadata = skin_metadata1, n_samples = 30)
#head(skin_downsampled_N30)


### select samples that are present in the merged betas
#skin_downsampled_betas_N30 <- merged_betas[, colnames(merged_betas) %in% skin_downsampled_N30$geo_accession]
#write.csv(skin_downsampled_N30, file = "skin_downsampled_N30.csv", row.names = FALSE)
#write.csv(skin_downsampled_betas_N30, file = "skin_downsampled_betas_N30.csv", row.names = FALSE)

#### getting the number of 


##### Reading and loading files 
skin_downsampled_betas_N30 <- read_csv("~/Phd_data/skin_downsampled_betas_N30.csv")
skin_downsampled_N30 <- read_csv("~/Phd_data/skin_downsampled_N30.csv")

correlation_test <- function(data, metadata) {
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
    
    # Initialize matrices for storing correlations and p-values
    n_features <- nrow(data)
    
    species_out <- data.frame(cpg_site = n_feature, correlation = NA, p_value = NA, conf_int_low = NA,conf_int_upp = NA )
    
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
      species_out$std_err <- (1- species_out$correlation) / (n_samples - 2)
    }
    
    # Store the results for the current species in the list
    species_results[[sp]] <- list()
  }
  
  return(species_results)
}

# Set seed for reproducibility
set.seed(123)

# Run the updated function
skin_downsampled_N30_corr <- correlation_test(
  data = skin_downsampled_betas_N30,  # The beta values of the CpG sites
  metadata = skin_downsampled_N30  # The metadata data frame 
)

cat("Saving the result to the folder")

saveRDS(skin_downsampled_N30_corr , file = "~/Phd_data/outputs/skin_downsampled_N30_corr.rds")

cat("Results have been saved to the current working directory.\n")


skin_downsampled_N30_corr <- readRDS("~/Phd_data/outputs/skin_downsampled_N30_corr.rds")
View(skin_downsampled_N30_corr)
## Extracting the data
### Extract the species in the list 
library(tibble)   # For rownames_to_column()
library(dplyr)    # For data manipulation

# Read metadata and get unique species names
skin_downsampled_N30 <- read.csv("~/Phd_data/skin_downsampled_N30.csv")

unique_species <- unique(skin_downsampled_N30$species)

# Initialize a named vector to store the number of significant sites for each species
num_sig_site <- setNames(rep(NA, length(unique_species)), unique_species)

# Read the correlation results (assuming this file contains the results for all species)
skin_corr <- readRDS("~/Phd_data/outputs/skin_downsampled_N30_corr.rds")

# Loop over each species (using a different variable name for clarity)
for (sp in unique_species) {
  
  cat("Processing species:", sp, "\n")
  
  # Subset correlation data for the current species.
  # Assuming your RDS list is structured by species, e.g., with names matching species
  if (!sp %in% names(skin_corr)) {
    warning(paste("Species", sp, "not found in the correlation data."))
    next
  }
  
  # Extract the data for the current species (assuming it is a data frame)
  sp_data <- as.data.frame(skin_corr[[sp]])
  
  # Optionally, if rownames represent CpG sites, convert them to a column:
  sp_data <- rownames_to_column(sp_data, var = "cpg_site")
  
  # Calculate the number of significant CpG sites with p_value <= 0.01
  Num_of_sig_cpg <- sum(sp_data$p_values <= 0.01, na.rm = TRUE)
  
  # Store the result in the vector using the species name as index
  num_sig_site[sp] <- Num_of_sig_cpg
}

# Convert to data frame if needed:
num_sig_site_df <- data.frame(
  species = names(num_sig_site),
  num_sig_sites = as.numeric(num_sig_site)
)
#add total cpg sites column as a number 
num_sig_site_df$total_cpg_sites <- "37554"
num_sig_site_df$total_cpg_sites <- as.numeric(as.character(num_sig_site_df$total_cpg_sites))
str(num_sig_site_df)

## deriving the percentage
num_sig_site_df$percentage <- (num_sig_site_df$num_sig_sites/num_sig_site_df$total_cpg_sites) * 100 
head(num_sig_site_df)
num_sig_site_df <- num_sig_site_df %>%
  mutate(percentage = round(percentage, 2))

### the plot
custom_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF00FF","#C71585", "#9ACD32", "#8B8B00","#BF3EFF",
  "#FF7F00", "#A65628", "#00F5FF", "#F781BF", "#999999", "#66C2A5", "#7CFC00","#FFE7BA","#CAE1FF",
  "#FC8D62", "#458B74", "#CDBA96", "#8B3626", "#006400", "#FFD700", "#FFFF33", "#8B5A2B", "#2F4F4F"
)

ggplot(num_sig_site_df, aes(x = species, y = num_sig_site)) +
  geom_bar(stat = "identity", fill= "#8B5A2B" ) +  #scale_fill_manual(values = custom_colors) +
  geom_text((aes(label = percentage)), vjust = 0.5, hjust = -0.1, size = 2, check_overlap = TRUE) +
  labs(
    title = "Number and Percentage of Significant CpG Sites at p = 0.01 in skin samples",
    x = "Mammalian Species",
    y = "Number of Significant CpG Sites"
  ) + coord_flip() + scale_y_continuous(limits = c(0, 20000), breaks = seq(0, 20000, by = 1000)) +
  theme(panel.grid = element_blank(),
        legend.title = element_text(family = "Roboto", color = "black", size = 12, face = "bold"),
        panel.border = element_rect(fill = "transparent", color = "black", size = 0.5)) +
  theme(axis.text.x = element_text(color = 1,
                                   size = 7,hjust = 0.6, angle = 75, vjust = 0.6),
        axis.text.y = element_text(color = 1,
                                   size = 10,
                                   hjust = 1)) 
View(num_sig_site_df)



