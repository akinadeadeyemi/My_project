########## Author: Adeyemi Akinade
########## Project: Comparative Epigenetics of Aging in mammalian species
########## Date: 22 February, 2025
########## Place: Clemson University, SC, USA.
########## Aim : To detect the false negative rate and the power of our sample size

# Here in this script, the analysis is based on determining the false negative rate using the final_results_with_pvalue as the 
# true estimate of the  correlation and then compared with the downsampled data for both blood and skin tissues.

# loading the data
final_results_with_pvalue <- read_csv("final_results_with_pvalue.csv")
View(final_results_with_pvalue)
blood_true_corr <- final_results_with_pvalue %>% filter(tissue=="Blood")
skin_true_corr <- final_results_with_pvalue %>% filter(tissue=="Skin")

## selecting only pvalues <= 0.05
blood_true_not_sig <- blood_true_corr %>% filter(species == "Homo sapiens") %>% filter(p_value < 0.05)
View(blood_true_not_sig)
## Calculating the fdr in blood samples only
## FDR in blood_sample_N = 12
library(tibble)   # For rownames_to_column()
library(dplyr)    # For data manipulation


# Read the correlation results (assuming this file contains the results for all species)
blood_sample_12_corr <- readRDS("~/Phd_data/outputs/blood_samples_12_correlation.rds")

# Initialize an empty list to store results for all species
############### Homo sapiens #########################

# Extract correlation and p-value data for the species
rownames(blood_sample_12_corr[["Homo sapiens"]][["correlations"]]) <- rownames(blood_metadata_betas)
sp_data_corr <- as.data.frame(blood_sample_12_corr[["Homo sapiens"]][["correlations"]])
rownames(blood_sample_12_corr[["Homo sapiens"]][["p_values"]]) <- rownames(blood_metadata_betas)
sp_data_pval <- as.data.frame(blood_sample_12_corr[["Homo sapiens"]][["p_values"]])

library(dplyr)
library(tibble)

# Extract p-values data frame for the current species
sp_data_pval <- as.data.frame(blood_sample_12_corr[["Homo sapiens"]][["p_values"]])

# Ensure blood_true_corr rownames are a column for joining
#blood_true_sig <- blood_true_sig %>%
# rownames_to_column(var = "Beta_ID")

# Initialize a result storage
homo_neg_fractions <- data.frame(Replicate = colnames(sp_data_pval), False_negative = NA, percent = NA)

# Loop over each column in sp_data_pval
for (i in seq_along(colnames(sp_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sp_data_pval)[i]
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sp_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) > 0.05) %>%  # Filter p > 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with blood_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, blood_true_not_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  #  Calculate the fraction of matches
  total_rows <- nrow(sp_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  # Store the result
  homo_neg_fractions$False_negative[i] <- fraction_matching
  homo_neg_fractions$percent[i] <- percent
}

# Print the final results
print(homo_neg_fractions)
view(homo_neg_fractions)
# Save results to a CSV file
write.csv(matching_fractions, "~/Phd_data/outputs/matching_fractions.csv", row.names = FALSE)
homo_neg_fractions$species <- "Homo sapiens"
homo_neg_fractions$sample_size <- "12"
homo_neg_fractions$tissue <- "blood"
homo_neg_fractions$power <- 1 - (homo_neg_fractions$False_negative)

##################### False negative for blood sample N = 36 in. homo sapiens  ##################
# Extract correlation and p-value data for the species
rownames(blood_sample_36_corr[["Homo sapiens"]][["correlations"]]) <- rownames(blood_metadata_betas)
sp_data_corr <- as.data.frame(blood_sample_36_corr[["Homo sapiens"]][["correlations"]])
rownames(blood_sample_36_corr[["Homo sapiens"]][["p_values"]]) <- rownames(blood_metadata_betas)
sp_data_pval <- as.data.frame(blood_sample_36_corr[["Homo sapiens"]][["p_values"]])

library(dplyr)
library(tibble)

# Extract p-values data frame for the current species
sp_data_pval <- as.data.frame(blood_sample_36_corr[["Homo sapiens"]][["p_values"]])

# Ensure blood_true_corr rownames are a column for joining
#blood_true_sig <- blood_true_sig %>%
# rownames_to_column(var = "Beta_ID")

# Initialize a result storage
homo_neg_fractions_36 <- data.frame(Replicate = colnames(sp_data_pval), False_negative = NA, percent = NA)

# Loop over each column in sp_data_pval
for (i in seq_along(colnames(sp_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sp_data_pval)[i]
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sp_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) > 0.05) %>%  # Filter p > 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with blood_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, blood_true_not_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  #  Calculate the fraction of matches
  total_rows <- nrow(sp_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  # Store the result
  homo_neg_fractions_36$False_negative[i] <- fraction_matching
  homo_neg_fractions_36$percent[i] <- percent
}

# Print the final results
print(homo_neg_fractions_36)
view(homo_neg_fractions_36)
# Save results to a CSV file
write.csv(matching_fractions_36, "~/Phd_data/outputs/matching_fractions_36.csv", row.names = FALSE)
homo_neg_fractions_36$species <- "Homo sapiens"
homo_neg_fractions_36$sample_size <- "36"
homo_neg_fractions_36$tissue <- "blood"
homo_neg_fractions_36$power <- 1 - (homo_neg_fractions_36$False_negative)

###################### False negative rate blood sample N = 60. ################

# Extract correlation and p-value data for the species
rownames(blood_sample_60_corr[["Homo sapiens"]][["correlations"]]) <- rownames(blood_metadata_betas)
sp_data_corr <- as.data.frame(blood_sample_60_corr[["Homo sapiens"]][["correlations"]])
rownames(blood_sample_60_corr[["Homo sapiens"]][["p_values"]]) <- rownames(blood_metadata_betas)
sp_data_pval <- as.data.frame(blood_sample_60_corr[["Homo sapiens"]][["p_values"]])

library(dplyr)
library(tibble)

# Extract p-values data frame for the current species
sp_data_pval <- as.data.frame(blood_sample_60_corr[["Homo sapiens"]][["p_values"]])

# Ensure blood_true_corr rownames are a column for joining
#blood_true_sig <- blood_true_sig %>%
# rownames_to_column(var = "Beta_ID")

# Initialize a result storage
homo_neg_fractions_60 <- data.frame(Replicate = colnames(sp_data_pval), False_negative = NA, percent = NA)

# Loop over each column in sp_data_pval
for (i in seq_along(colnames(sp_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sp_data_pval)[i]
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sp_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) > 0.05) %>%  # Filter p > 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with blood_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, blood_true_not_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  #  Calculate the fraction of matches
  total_rows <- nrow(sp_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  # Store the result
  homo_neg_fractions_60$False_negative[i] <- fraction_matching
  homo_neg_fractions_60$percent[i] <- percent
}

# Print the final results
print(homo_neg_fractions_60)
view(homo_neg_fractions_60)
# Save results to a CSV file
write.csv(matching_fractions_60, "~/Phd_data/outputs/matching_fractions_60.csv", row.names = FALSE)
homo_neg_fractions_60$species <- "Homo sapiens"
homo_neg_fractions_60$sample_size <- "60"
homo_neg_fractions_60$tissue <- "blood"
homo_neg_fractions_60$power <- 1 - (homo_neg_fractions_60$False_negative)




######################### False negative rate blood samples N = 90 ###############################
# Extract correlation and p-value data for the species
rownames(blood_sample_90_corr[["Homo sapiens"]][["correlations"]]) <- rownames(blood_metadata_betas)
sp_data_corr <- as.data.frame(blood_sample_90_corr[["Homo sapiens"]][["correlations"]])
rownames(blood_sample_90_corr[["Homo sapiens"]][["p_values"]]) <- rownames(blood_metadata_betas)
sp_data_pval <- as.data.frame(blood_sample_90_corr[["Homo sapiens"]][["p_values"]])

library(dplyr)
library(tibble)

# Extract p-values data frame for the current species
sp_data_pval <- as.data.frame(blood_sample_90_corr[["Homo sapiens"]][["p_values"]])

# Ensure blood_true_corr rownames are a column for joining
#blood_true_sig <- blood_true_sig %>%
# rownames_to_column(var = "Beta_ID")

# Initialize a result storage
homo_neg_fractions_90 <- data.frame(Replicate = colnames(sp_data_pval), False_negative = NA, percent = NA)

# Loop over each column in sp_data_pval
for (i in seq_along(colnames(sp_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sp_data_pval)[i]
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sp_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) > 0.05) %>%  # Filter p > 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with blood_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, blood_true_not_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  #  Calculate the fraction of matches
  total_rows <- nrow(sp_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  # Store the result
  homo_neg_fractions_90$False_negative[i] <- fraction_matching
  homo_neg_fractions_90$percent[i] <- percent
}

# Print the final results
print(homo_neg_fractions_90)
view(homo_neg_fractions_90)
# Save results to a CSV file
write.csv(matching_fractions_90, "~/Phd_data/outputs/matching_fractions_90.csv", row.names = FALSE)
homo_neg_fractions_90$species <- "Homo sapiens"
homo_neg_fractions_90$sample_size <- "90"
homo_neg_fractions_90$tissue <- "blood"
homo_neg_fractions_90$power <- 1 - (homo_neg_fractions_90$False_negative)

###### plot the false negative rate and power
## Homo sapiens
homo_neg_combined_fdr <- rbind(homo_neg_fractions, homo_neg_fractions_36, homo_neg_fractions_60, homo_neg_fractions_90)
view(homo_neg_combined_fdr)

## False negative distribution plot
ggplot(homo_neg_combined_fdr, aes(x = factor(sample_size), y = percent, fill = sample_size)) + 
  geom_violin(alpha = 0.7) + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, alpha = 0.5) + 
  theme(panel.background = element_rect(color = "black", size = 0.2)) + 
  labs(x = "Sample Sizes", 
       y = "False Negative Rate (in %)", 
       title = "The False Negative Rate using downsampled DNAm blood samples from Homo sapiens") +
  theme_minimal() + 
  theme(legend.position = "none")  # Hide legend if not needed

## Power distribution plot
ggplot(homo_neg_combined_fdr, aes(x= power,  fill=sample_size)) +
  geom_density(alpha=0.2) + theme_bw() + scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme(panel.background = element_rect(color = "black", # Color of the border
                                        size = 0.2)) + 
  labs(x = "power (1 - β)", 
       title = "The power of different sample sizes of DNAm blood samples from Homo sapiens")

## box plot of power
ggplot(homo_neg_combined_fdr, aes(x = sample_size, y=power,  fill=sample_size)) +
  geom_violin(alpha = 0.7) + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, alpha = 0.5) + 
  scale_x_discrete() +
  theme(panel.background = element_rect(color = "black", size = 0.2)) + 
  labs(title = "Violin plot of sample size power from DNA methylation from Homo sapiens blood sample", x = "Sample size", 
       y = "Power (1 - β)") +
  theme_minimal() + theme(panel.background = element_rect(color = "black", # Color of the border
                                                        size = 0.2))

########################### canis lupus familiaris N = 12 #################################  
canis_true_not_sig <- blood_true_corr %>% filter(species == "Canis lupus familiaris") %>% filter(p_value < 0.05)
View(canis_true_not_sig)

# Extract correlation and p-value data for the species
rownames(blood_sample_12_corr[["Canis lupus familiaris"]][["correlations"]]) <- rownames(blood_metadata_betas)
sp_data_corr <- as.data.frame(blood_sample_12_corr[["Canis lupus familiaris"]][["correlations"]])
rownames(blood_sample_12_corr[["Canis lupus familiaris"]][["p_values"]]) <- rownames(blood_metadata_betas)
sp_data_pval <- as.data.frame(blood_sample_12_corr[["Canis lupus familiaris"]][["p_values"]])

library(dplyr)
library(tibble)

# Extract p-values data frame for the current species
sp_data_pval <- as.data.frame(blood_sample_12_corr[["Canis lupus familiaris"]][["p_values"]])

# Ensure blood_true_corr rownames are a column for joining
#blood_true_sig <- blood_true_sig %>%
# rownames_to_column(var = "Beta_ID")

# Initialize a result storage
canis_neg_fractions <- data.frame(Replicate = colnames(sp_data_pval), False_negative = NA, percent = NA)

# Loop over each column in sp_data_pval
for (i in seq_along(colnames(sp_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sp_data_pval)[i]
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sp_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) > 0.05) %>%  # Filter p < 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with blood_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, canis_true_not_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  #  Calculate the fraction of matches
  total_rows <- nrow(sp_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  # Store the result
  canis_neg_fractions$False_negative[i] <- fraction_matching
  canis_neg_fractions$percent[i] <- percent
}

# Print the final results
print(canis_neg_fractions)
view(canis_neg_fractions)
canis_neg_fractions$species <- "Canis lupus familiaris"
canis_neg_fractions$sample_size <- "12"
canis_neg_fractions$tissue <- "blood"
canis_neg_fractions$power <- 1 - (canis_neg_fractions$False_negative)

########################### canis lupus familiaris N = 36 #################################  
canis_true_not_sig <- blood_true_corr %>% filter(species == "Canis lupus familiaris") %>% filter(p_value < 0.05)
View(canis_true_not_sig)

# Extract correlation and p-value data for the species
rownames(blood_sample_36_corr[["Canis lupus familiaris"]][["correlations"]]) <- rownames(blood_metadata_betas)
sp_data_corr <- as.data.frame(blood_sample_36_corr[["Canis lupus familiaris"]][["correlations"]])
rownames(blood_sample_36_corr[["Canis lupus familiaris"]][["p_values"]]) <- rownames(blood_metadata_betas)
sp_data_pval <- as.data.frame(blood_sample_36_corr[["Canis lupus familiaris"]][["p_values"]])

library(dplyr)
library(tibble)

# Extract p-values data frame for the current species
sp_data_pval <- as.data.frame(blood_sample_36_corr[["Canis lupus familiaris"]][["p_values"]])

# Ensure blood_true_corr rownames are a column for joining
#blood_true_sig <- blood_true_sig %>%
# rownames_to_column(var = "Beta_ID")

# Initialize a result storage
canis_neg_fractions_36 <- data.frame(Replicate = colnames(sp_data_pval), False_negative = NA, percent = NA)

# Loop over each column in sp_data_pval
for (i in seq_along(colnames(sp_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sp_data_pval)[i]
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sp_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) > 0.05) %>%  # Filter p < 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with blood_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, canis_true_not_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  #  Calculate the fraction of matches
  total_rows <- nrow(sp_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  # Store the result
  canis_neg_fractions_36$False_negative[i] <- fraction_matching
  canis_neg_fractions_36$percent[i] <- percent
}

# Print the final results
print(canis_neg_fractions_36)
view(canis_neg_fractions_36)
canis_neg_fractions_36$species <- "Canis lupus familiaris"
canis_neg_fractions_36$sample_size <- "36"
canis_neg_fractions_36$tissue <- "blood"
canis_neg_fractions_36$power <- 1 - (canis_neg_fractions_36$False_negative)


########################### canis lupus familiaris N = 60 #################################  
canis_true_not_sig <- blood_true_corr %>% filter(species == "Canis lupus familiaris") %>% filter(p_value < 0.05)
View(canis_true_not_sig)

# Extract correlation and p-value data for the species
rownames(blood_sample_60_corr[["Canis lupus familiaris"]][["correlations"]]) <- rownames(blood_metadata_betas)
sp_data_corr <- as.data.frame(blood_sample_60_corr[["Canis lupus familiaris"]][["correlations"]])
rownames(blood_sample_60_corr[["Canis lupus familiaris"]][["p_values"]]) <- rownames(blood_metadata_betas)
sp_data_pval <- as.data.frame(blood_sample_60_corr[["Canis lupus familiaris"]][["p_values"]])

library(dplyr)
library(tibble)

# Extract p-values data frame for the current species
sp_data_pval <- as.data.frame(blood_sample_60_corr[["Canis lupus familiaris"]][["p_values"]])

# Ensure blood_true_corr rownames are a column for joining
#blood_true_sig <- blood_true_sig %>%
# rownames_to_column(var = "Beta_ID")

# Initialize a result storage
canis_neg_fractions_60 <- data.frame(Replicate = colnames(sp_data_pval), False_negative = NA, percent = NA)

# Loop over each column in sp_data_pval
for (i in seq_along(colnames(sp_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sp_data_pval)[i]
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sp_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) > 0.05) %>%  # Filter p < 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with blood_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, canis_true_not_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  #  Calculate the fraction of matches
  total_rows <- nrow(sp_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  # Store the result
  canis_neg_fractions_60$False_negative[i] <- fraction_matching
  canis_neg_fractions_60$percent[i] <- percent
}

# Print the final results
print(canis_neg_fractions_60)
view(canis_neg_fractions_60)
canis_neg_fractions_60$species <- "Canis lupus familiaris"
canis_neg_fractions_60$sample_size <- "60"
canis_neg_fractions_60$tissue <- "blood"
canis_neg_fractions_60$power <- 1 - (canis_neg_fractions_60$False_negative)

########################### canis lupus familiaris N = 90 #################################  
canis_true_not_sig <- blood_true_corr %>% filter(species == "Canis lupus familiaris") %>% filter(p_value < 0.05)
View(canis_true_not_sig)

# Extract correlation and p-value data for the species
rownames(blood_sample_90_corr[["Canis lupus familiaris"]][["correlations"]]) <- rownames(blood_metadata_betas)
sp_data_corr <- as.data.frame(blood_sample_90_corr[["Canis lupus familiaris"]][["correlations"]])
rownames(blood_sample_90_corr[["Canis lupus familiaris"]][["p_values"]]) <- rownames(blood_metadata_betas)
sp_data_pval <- as.data.frame(blood_sample_90_corr[["Canis lupus familiaris"]][["p_values"]])

library(dplyr)
library(tibble)

# Extract p-values data frame for the current species
sp_data_pval <- as.data.frame(blood_sample_90_corr[["Canis lupus familiaris"]][["p_values"]])

# Ensure blood_true_corr rownames are a column for joining
#blood_true_sig <- blood_true_sig %>%
# rownames_to_column(var = "Beta_ID")

# Initialize a result storage
canis_neg_fractions_90 <- data.frame(Replicate = colnames(sp_data_pval), False_negative = NA, percent = NA)

# Loop over each column in sp_data_pval
for (i in seq_along(colnames(sp_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sp_data_pval)[i]
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sp_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) > 0.05) %>%  # Filter p < 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with blood_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, canis_true_not_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  #  Calculate the fraction of matches
  total_rows <- nrow(sp_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  # Store the result
  canis_neg_fractions_90$False_negative[i] <- fraction_matching
  canis_neg_fractions_90$percent[i] <- percent
}

# Print the final results
print(canis_neg_fractions_90)
view(canis_neg_fractions_90)
canis_neg_fractions_90$species <- "Canis lupus familiaris"
canis_neg_fractions_90$sample_size <- "90"
canis_neg_fractions_90$tissue <- "blood"
canis_neg_fractions_90$power <- 1 - (canis_neg_fractions_90$False_negative)

canis_neg_combined_fdr <- rbind(canis_neg_fractions, canis_neg_fractions_36, canis_neg_fractions_60, canis_neg_fractions_90)
view(canis_neg_combined_fdr)

## False negative distribution plot
ggplot(canis_neg_combined_fdr, aes(x = factor(sample_size), y = percent, fill = sample_size)) + 
  geom_violin(alpha = 0.7) + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, alpha = 0.5) + 
  scale_x_discrete() +
  theme(panel.background = element_rect(color = "black", size = 0.2)) + 
  labs(x = "Sample Sizes", 
       y = "False Negative Rate (in %)", 
       title = "The False Negative Rate using Downsampled DNAm Blood Samples from Canis lupus familiaris") +
  theme_minimal() + 
  theme(legend.position = "none")  # Hide legend if not needed

## Power distribution plot
ggplot(canis_neg_combined_fdr, aes(x= power,  fill=sample_size)) +
  geom_density(alpha=0.2) + theme_bw() + scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme(panel.background = element_rect(color = "black", # Color of the border
                                        size = 0.2)) + 
  labs(x = "power (1 - β)", 
       title = "The power of different sample sizes of DNAm blood samples from Canis lupus familiaris")

## box plot of power
ggplot(canis_neg_combined_fdr, aes(x = factor(sample_size), y = power, fill = sample_size)) + 
  geom_violin(alpha = 0.7) + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, alpha = 0.5) + 
  scale_x_discrete() +
  theme(panel.background = element_rect(color = "black", size = 0.2)) + 
  labs(x = "Sample Sizes", 
       y = "False Negative Rate (in %)", 
       title = "The False Negative Rate using Downsampled DNAm Blood Samples from Canis lupus familiaris") +
  theme_minimal() + 
  theme(legend.position = "none")  +
  labs(title = "Boxplot of sample size power from DNA methylation from Canis lupus familiaris blood sample", 
       x = "Sample size", y = "Power(1 - β)") +
  theme_minimal() + theme(panel.background = element_rect(color = "black", # Color of the border
                                                          size = 0.2))









############################################## FALSE NEGATIVE RATE FOR SKIN SAMPLES ###################################################
# loading the data
final_results_with_pvalue <- read_csv("final_results_with_pvalue.csv")
View(final_results_with_pvalue)
blood_true_corr <- final_results_with_pvalue %>% filter(tissue=="Blood")
skin_true_corr <- final_results_with_pvalue %>% filter(tissue=="Skin")

## selecting only pvalues <= 0.05
skin_true_not_sig <- skin_true_corr %>% filter(species == "Homo sapiens") %>% filter(p_value < 0.05)
View(skin_true_not_sig)

## Calculating the fdr in blood samples only
## FDR in blood_sample_N = 12
library(tibble)   # For rownames_to_column()
library(dplyr)    # For data manipulation


# Read the correlation results (assuming this file contains the results for all species)
skin_sample_14_corr <- readRDS("~/Phd_data/outputs/skin_sample_14_corr_test.rds")

# Initialize an empty list to store results for all species
############### Homo sapiens N = 14 #########################

# Extract correlation and p-value data for the species
rownames(skin_sample_14_corr[["Homo sapiens"]][["correlations"]]) <- rownames(blood_metadata_betas)
sk_data_corr <- as.data.frame(skin_sample_14_corr[["Homo sapiens"]][["correlations"]])
rownames(skin_sample_14_corr[["Homo sapiens"]][["p_values"]]) <- rownames(blood_metadata_betas)
sk_data_pval <- as.data.frame(skin_sample_14_corr[["Homo sapiens"]][["p_values"]])

library(dplyr)
library(tibble)

# Extract p-values data frame for the current species
sk_data_pval <- as.data.frame(skin_sample_14_corr[["Homo sapiens"]][["p_values"]])

# Ensure blood_true_corr rownames are a column for joining
#blood_true_sig <- blood_true_sig %>%
# rownames_to_column(var = "Beta_ID")

# Initialize a result storage
homo_skin_neg_14 <- data.frame(Replicate = colnames(sk_data_pval), False_negative = NA, percent = NA)

# Loop over each column in sp_data_pval
for (i in seq_along(colnames(sk_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sk_data_pval)[i]
  print(current_column)
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sk_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) > 0.05) %>%  # Filter p > 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with blood_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, skin_true_not_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  
  #  Calculate the fraction of matches
  total_rows <- nrow(sk_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  
  # Store the result
  homo_skin_neg_14$False_negative[i] <- fraction_matching
  homo_skin_neg_14$percent[i] <- percent
  
  cat("calculation is completed...")
}

# Print the final results
print(homo_skin_neg_14)
view(homo_skin_neg_14)
# Save results to a CSV file
homo_skin_neg_14$species <- "Homo sapiens"
homo_skin_neg_14$sample_size <- "14"
homo_skin_neg_14$tissue <- "skin"
homo_skin_neg_14$power <- 1 - (homo_skin_neg_14$False_negative) 


############### Homo sapiens N = 32 #########################
# Read the correlation results (assuming this file contains the results for all species)
skin_sample_32_corr <- readRDS("~/Phd_data/outputs/skin_sample_32_correlation.rds")

# Extract correlation and p-value data for the species
rownames(skin_sample_32_corr[["Homo sapiens"]][["correlations"]]) <- rownames(blood_metadata_betas)
sk_data_corr <- as.data.frame(skin_sample_32_corr[["Homo sapiens"]][["correlations"]])
rownames(skin_sample_32_corr[["Homo sapiens"]][["p_values"]]) <- rownames(blood_metadata_betas)
sk_data_pval <- as.data.frame(skin_sample_32_corr[["Homo sapiens"]][["p_values"]])

library(dplyr)
library(tibble)

# Extract p-values data frame for the current species
sk_data_pval <- as.data.frame(skin_sample_32_corr[["Homo sapiens"]][["p_values"]])

# Ensure blood_true_corr rownames are a column for joining
#blood_true_sig <- blood_true_sig %>%
# rownames_to_column(var = "Beta_ID")

# Initialize a result storage
homo_skin_neg_32 <- data.frame(Replicate = colnames(sk_data_pval), False_negative = NA, percent = NA)

# Loop over each column in sp_data_pval
for (i in seq_along(colnames(sk_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sk_data_pval)[i]
  print(current_column)
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sk_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) > 0.05) %>%  # Filter p > 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with blood_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, skin_true_not_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  
  #  Calculate the fraction of matches
  total_rows <- nrow(sk_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  
  # Store the result
  homo_skin_neg_32$False_negative[i] <- fraction_matching
  homo_skin_neg_32$percent[i] <- percent
  
  cat("calculation is completed...")
}

# Print the final results
print(homo_skin_neg_32)
view(homo_skin_neg_32)
# Save results to a CSV file
homo_skin_neg_32$species <- "Homo sapiens"
homo_skin_neg_32$sample_size <- "32"
homo_skin_neg_32$tissue <- "skin"
homo_skin_neg_32$power <- 1 - (homo_skin_neg_32$False_negative) 



############### Homo sapiens N = 50 #########################
# Read the correlation results (assuming this file contains the results for all species)
skin_sample_50_corr <- readRDS("~/Phd_data/outputs/skin_sample_50_correlation.rds")

# Extract correlation and p-value data for the species
rownames(skin_sample_50_corr[["Homo sapiens"]][["correlations"]]) <- rownames(blood_metadata_betas)
sk_data_corr <- as.data.frame(skin_sample_50_corr[["Homo sapiens"]][["correlations"]])
rownames(skin_sample_50_corr[["Homo sapiens"]][["p_values"]]) <- rownames(blood_metadata_betas)
sk_data_pval <- as.data.frame(skin_sample_50_corr[["Homo sapiens"]][["p_values"]])

library(dplyr)
library(tibble)

# Extract p-values data frame for the current species
sk_data_pval <- as.data.frame(skin_sample_50_corr[["Homo sapiens"]][["p_values"]])

# Ensure blood_true_corr rownames are a column for joining
#blood_true_sig <- blood_true_sig %>%
# rownames_to_column(var = "Beta_ID")

# Initialize a result storage
homo_skin_neg_50 <- data.frame(Replicate = colnames(sk_data_pval), False_negative = NA, percent = NA)

# Loop over each column in sp_data_pval
for (i in seq_along(colnames(sk_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sk_data_pval)[i]
  print(current_column)
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sk_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) > 0.05) %>%  # Filter p > 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with blood_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, skin_true_not_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  
  #  Calculate the fraction of matches
  total_rows <- nrow(sk_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  
  # Store the result
  homo_skin_neg_50$False_negative[i] <- fraction_matching
  homo_skin_neg_50$percent[i] <- percent
  
  cat("calculation is completed...")
}

# Print the final results
print(homo_skin_neg_50)
view(homo_skin_neg_50)
# Save results to a CSV file
homo_skin_neg_50$species <- "Homo sapiens"
homo_skin_neg_50$sample_size <- "50"
homo_skin_neg_50$tissue <- "skin"
homo_skin_neg_50$power <- 1 - (homo_skin_neg_50$False_negative) 


############### Homo sapiens N = 70 #########################
# Read the correlation results (assuming this file contains the results for all species)
skin_sample_70_corr <- readRDS("~/Phd_data/outputs/skin_sample_70_correlation.rds")

# Extract correlation and p-value data for the species
rownames(skin_sample_70_corr[["Homo sapiens"]][["correlations"]]) <- rownames(blood_metadata_betas)
sk_data_corr <- as.data.frame(skin_sample_70_corr[["Homo sapiens"]][["correlations"]])
rownames(skin_sample_70_corr[["Homo sapiens"]][["p_values"]]) <- rownames(blood_metadata_betas)
sk_data_pval <- as.data.frame(skin_sample_70_corr[["Homo sapiens"]][["p_values"]])

library(dplyr)
library(tibble)

# Extract p-values data frame for the current species
sk_data_pval <- as.data.frame(skin_sample_70_corr[["Homo sapiens"]][["p_values"]])

# Ensure blood_true_corr rownames are a column for joining
#blood_true_sig <- blood_true_sig %>%
# rownames_to_column(var = "Beta_ID")

# Initialize a result storage
homo_skin_neg_70 <- data.frame(Replicate = colnames(sk_data_pval), False_negative = NA, percent = NA)

# Loop over each column in sp_data_pval
for (i in seq_along(colnames(sk_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sk_data_pval)[i]
  print(current_column)
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sk_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) > 0.05) %>%  # Filter p > 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with blood_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, skin_true_not_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  
  #  Calculate the fraction of matches
  total_rows <- nrow(sk_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  
  # Store the result
  homo_skin_neg_70$False_negative[i] <- fraction_matching
  homo_skin_neg_70$percent[i] <- percent
  
  cat("calculation is completed...")
}

# Print the final results
print(homo_skin_neg_70)
view(homo_skin_neg_70)
# Save results to a CSV file
homo_skin_neg_70$species <- "Homo sapiens"
homo_skin_neg_70$sample_size <- "70"
homo_skin_neg_70$tissue <- "skin"
homo_skin_neg_70$power <- 1 - (homo_skin_neg_70$False_negative) 


### combining the dataframes
homo_skin_neg_combined <- rbind(homo_skin_neg_14, homo_skin_neg_32, homo_skin_neg_50, homo_skin_neg_70)
view(homo_skin_neg_combined)

## False negative distribution plot
ggplot(homo_skin_neg_combined, aes(x = factor(sample_size), y = percent, fill = sample_size)) + 
  geom_violin(alpha = 0.7) + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, alpha = 0.5) + 
  scale_x_discrete() +
  theme(panel.background = element_rect(color = "black", size = 0.2)) + 
  labs(x = "Sample Sizes", 
       y = "False Negative Rate (in %)", 
       title = "The False Negative Rate using downsampled DNAm skin samples from Homo sapiens") +
  theme_minimal() + 
  theme(legend.position = "none")  # Hide legend if not needed

## Power distribution plot
ggplot(homo_skin_neg_combined, aes(x= power,  fill=sample_size)) +
  geom_density(alpha=0.2) + theme_bw() + scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme(panel.background = element_rect(color = "black", # Color of the border
                                        size = 0.2)) + 
  labs(x = "power (1 - β)", 
       title = "The power of different sample sizes of DNAm blood samples from Canis lupus familiaris")

## box plot of power
ggplot(homo_skin_neg_combined, aes(x = factor(sample_size), y = power, fill = sample_size)) + 
  geom_violin(alpha = 0.7) + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, alpha = 0.5) + 
  scale_x_discrete() +
  theme(panel.background = element_rect(color = "black", size = 0.2)) + 
  labs(x = "Sample Sizes", 
       y = "False Negative Rate (in %)", 
       title = "The False Negative Rate using Downsampled DNAm Blood Samples from Canis lupus familiaris") +
  theme_minimal() + 
  theme(legend.position = "none")  +
  labs(title = "Violin plot of sample size power from DNA methylation from Homo sapiens skin sample", 
       x = "Sample size", y = "Power(1 - β)") +
  theme_minimal() + theme(panel.background = element_rect(color = "black", # Color of the border
                                                          size = 0.2))




############################ Tursiops aduncus skin samples N = 14 ##########################
## selecting only pvalues <= 0.05
tursiops_skin_true_not_sig <- skin_true_corr %>% filter(species == "Tursiops aduncus") %>% filter(p_value < 0.05)
View(tursiops_skin_true_not_sig)



# Read the correlation results (assuming this file contains the results for all species)
skin_sample_14_corr <- readRDS("~/Phd_data/outputs/skin_sample_14_correlation.rds")

# Extract correlation and p-value data for the species
rownames(skin_sample_14_corr[["Tursiops aduncus"]][["correlations"]]) <- rownames(blood_metadata_betas)
sk_data_corr <- as.data.frame(skin_sample_14_corr[["Tursiops aduncus"]][["correlations"]])
rownames(skin_sample_14_corr[["Tursiops aduncus"]][["p_values"]]) <- rownames(blood_metadata_betas)
sk_data_pval <- as.data.frame(skin_sample_14_corr[["Tursiops aduncus"]][["p_values"]])

library(dplyr)
library(tibble)

# Extract p-values data frame for the current species
sk_data_pval <- as.data.frame(skin_sample_14_corr[["Tursiops aduncus"]][["p_values"]])

# Ensure blood_true_corr rownames are a column for joining
#blood_true_sig <- blood_true_sig %>%
# rownames_to_column(var = "Beta_ID")

# Initialize a result storage
Tursiops_skin_neg_14 <- data.frame(Replicate = colnames(sk_data_pval), False_negative = NA, percent = NA)

# Loop over each column in sp_data_pval
for (i in seq_along(colnames(sk_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sk_data_pval)[i]
  print(current_column)
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sk_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) > 0.05) %>%  # Filter p > 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with blood_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, tursiops_skin_true_not_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  
  #  Calculate the fraction of matches
  total_rows <- nrow(sk_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  
  # Store the result
  Tursiops_skin_neg_14$False_negative[i] <- fraction_matching
  Tursiops_skin_neg_14$percent[i] <- percent
  
  cat("calculation is completed...")
}

# Print the final results
print(Tursiops_skin_neg_14)
view(Tursiops_skin_neg_14)
# Save results to a CSV file
Tursiops_skin_neg_14$species <- "Tursiops aduncus"
Tursiops_skin_neg_14$sample_size <- "14"
Tursiops_skin_neg_14$tissue <- "skin"
Tursiops_skin_neg_14$power <- 1 - (Tursiops_skin_neg_14$False_negative) 



############### Tursiops aduncus N = 32 #########################
# Read the correlation results (assuming this file contains the results for all species)
skin_sample_32_corr <- readRDS("~/Phd_data/outputs/skin_sample_32_correlation.rds")

# Extract correlation and p-value data for the species
rownames(skin_sample_32_corr[["Tursiops aduncus"]][["correlations"]]) <- rownames(blood_metadata_betas)
sk_data_corr <- as.data.frame(skin_sample_32_corr[["Tursiops aduncus"]][["correlations"]])
rownames(skin_sample_32_corr[["Tursiops aduncus"]][["p_values"]]) <- rownames(blood_metadata_betas)
sk_data_pval <- as.data.frame(skin_sample_32_corr[["Tursiops aduncus"]][["p_values"]])

library(dplyr)
library(tibble)

# Extract p-values data frame for the current species
sk_data_pval <- as.data.frame(skin_sample_32_corr[["Tursiops aduncus"]][["p_values"]])

# Ensure blood_true_corr rownames are a column for joining
#blood_true_sig <- blood_true_sig %>%
# rownames_to_column(var = "Beta_ID")

# Initialize a result storage
Tursiops_skin_neg_32 <- data.frame(Replicate = colnames(sk_data_pval), False_negative = NA, percent = NA)

# Loop over each column in sp_data_pval
for (i in seq_along(colnames(sk_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sk_data_pval)[i]
  print(current_column)
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sk_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) > 0.05) %>%  # Filter p > 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with blood_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, tursiops_skin_true_not_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  
  #  Calculate the fraction of matches
  total_rows <- nrow(sk_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  
  # Store the result
  Tursiops_skin_neg_32$False_negative[i] <- fraction_matching
  Tursiops_skin_neg_32$percent[i] <- percent
  
  cat("calculation is completed...")
}

# Print the final results
print(Tursiops_skin_neg_32)
view(Tursiops_skin_neg_32)
# Save results to a CSV file
Tursiops_skin_neg_32$species <- "Tursiops aduncus"
Tursiops_skin_neg_32$sample_size <- "32"
Tursiops_skin_neg_32$tissue <- "skin"
Tursiops_skin_neg_32$power <- 1 - (Tursiops_skin_neg_32$False_negative) 



############### Tursiops aduncus N = 50 #########################
# Read the correlation results (assuming this file contains the results for all species)
skin_sample_50_corr <- readRDS("~/Phd_data/outputs/skin_sample_50_correlation.rds")

# Extract correlation and p-value data for the species
rownames(skin_sample_50_corr[["Tursiops aduncus"]][["correlations"]]) <- rownames(blood_metadata_betas)
sk_data_corr <- as.data.frame(skin_sample_50_corr[["Tursiops aduncus"]][["correlations"]])
rownames(skin_sample_50_corr[["Tursiops aduncus"]][["p_values"]]) <- rownames(blood_metadata_betas)
sk_data_pval <- as.data.frame(skin_sample_50_corr[["Tursiops aduncus"]][["p_values"]])

library(dplyr)
library(tibble)

# Extract p-values data frame for the current species
sk_data_pval <- as.data.frame(skin_sample_50_corr[["Tursiops aduncus"]][["p_values"]])

# Ensure blood_true_corr rownames are a column for joining
#blood_true_sig <- blood_true_sig %>%
# rownames_to_column(var = "Beta_ID")

# Initialize a result storage
Tursiops_skin_neg_50 <- data.frame(Replicate = colnames(sk_data_pval), False_negative = NA, percent = NA)

# Loop over each column in sp_data_pval
for (i in seq_along(colnames(sk_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sk_data_pval)[i]
  print(current_column)
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sk_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) > 0.05) %>%  # Filter p > 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with blood_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, tursiops_skin_true_not_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  
  #  Calculate the fraction of matches
  total_rows <- nrow(sk_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  
  # Store the result
  Tursiops_skin_neg_50$False_negative[i] <- fraction_matching
  Tursiops_skin_neg_50$percent[i] <- percent
  
  cat("calculation is completed...")
}

# Print the final results
print(Tursiops_skin_neg_50)
view(Tursiops_skin_neg_50)
# Save results to a CSV file
Tursiops_skin_neg_50$species <- "Tursiops aduncus"
Tursiops_skin_neg_50$sample_size <- "50"
Tursiops_skin_neg_50$tissue <- "skin"
Tursiops_skin_neg_50$power <- 1 - (Tursiops_skin_neg_50$False_negative) 


############### Tursiops aduncus N = 70 #########################
# Read the correlation results (assuming this file contains the results for all species)
skin_sample_70_corr <- readRDS("~/Phd_data/outputs/skin_sample_70_correlation.rds")

# Extract correlation and p-value data for the species
rownames(skin_sample_70_corr[["Tursiops aduncus"]][["correlations"]]) <- rownames(blood_metadata_betas)
sk_data_corr <- as.data.frame(skin_sample_70_corr[["Tursiops aduncus"]][["correlations"]])
rownames(skin_sample_70_corr[["Tursiops aduncus"]][["p_values"]]) <- rownames(blood_metadata_betas)
sk_data_pval <- as.data.frame(skin_sample_70_corr[["Tursiops aduncus"]][["p_values"]])

library(dplyr)
library(tibble)

# Extract p-values data frame for the current species
sk_data_pval <- as.data.frame(skin_sample_70_corr[["Tursiops aduncus"]][["p_values"]])

# Ensure blood_true_corr rownames are a column for joining
#blood_true_sig <- blood_true_sig %>%
# rownames_to_column(var = "Beta_ID")

# Initialize a result storage
Tursiops_skin_neg_70 <- data.frame(Replicate = colnames(sk_data_pval), False_negative = NA, percent = NA)

# Loop over each column in sp_data_pval
for (i in seq_along(colnames(sk_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sk_data_pval)[i]
  print(current_column)
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sk_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) > 0.05) %>%  # Filter p > 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with blood_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, tursiops_skin_true_not_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  
  #  Calculate the fraction of matches
  total_rows <- nrow(sk_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  
  # Store the result
  Tursiops_skin_neg_70$False_negative[i] <- fraction_matching
  Tursiops_skin_neg_70$percent[i] <- percent
  
  cat("calculation is completed...")
}

# Print the final results
print(Tursiops_skin_neg_70)
view(Tursiops_skin_neg_70)
# Save results to a CSV file
Tursiops_skin_neg_70$species <- "Tursiops aduncus"
Tursiops_skin_neg_70$sample_size <- "70"
Tursiops_skin_neg_70$tissue <- "skin"
Tursiops_skin_neg_70$power <- 1 - (Tursiops_skin_neg_70$False_negative) 


### combining the dataframes
tursiops_skin_neg_combined <- rbind(Tursiops_skin_neg_14, Tursiops_skin_neg_32, Tursiops_skin_neg_50, Tursiops_skin_neg_70)
view(tursiops_skin_neg_combined)

## False negative distribution plot


ggplot(tursiops_skin_neg_combined, aes(x = factor(sample_size), y = percent, fill = sample_size)) + 
  geom_violin(alpha = 0.7) + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, alpha = 0.5) + 
  scale_x_discrete() +
  theme(panel.background = element_rect(color = "black", size = 0.2)) + 
  labs(x = "Sample Sizes", 
       y = "False Negative Rate (in %)", 
       title = "The False Negative Rate using downsampled DNAm skin samples from Tursiops aduncus") +
  theme_minimal() + 
  theme(legend.position = "none")  # Hide legend if not needed

## Power distribution plot
ggplot(trusiops_skin_neg_combined, aes(x= power,  fill=sample_size)) +
  geom_density(alpha=0.2) + theme_bw() + scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme(panel.background = element_rect(color = "black", # Color of the border
                                        size = 0.2)) + 
  labs(x = "power (1 - β)", 
       title = "The power of different sample sizes of DNAm blood samples from Tursiops aduncus")

## box plot of power
ggplot(tursiops_skin_neg_combined, aes(x = factor(sample_size), y = power, fill = sample_size)) + 
  geom_violin(alpha = 0.7) + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1.5, alpha = 0.5) + 
  scale_x_discrete() +
  theme(panel.background = element_rect(color = "black", size = 0.2)) + 
  theme_minimal() + 
  theme(legend.position = "none")  +
  labs(title = "Violin plot of sample size power from DNA methylation from Tursiops aduncus skin sample", 
       x = "Sample size", y = "Power(1 - β)") +
  theme_minimal() + theme(panel.background = element_rect(color = "black", # Color of the border
                                                          size = 0.2))




























