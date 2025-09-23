########## Author: Adeyemi Akinade
########## Project: Comparative Epigenetics of Aging in mammalian species
########## Date: 22 February, 2025
########## Place: Clemson University, SC, USA.
########## Aim : To detect the faslse positive rate and false negative rate and the power of our sample size

# Here in this script, the analysis is based on determining the false positive rate using the final_results_with_pvalue as the 
# true estimate of the  correlation and then compared with the downsampled data for both blood and skin tissues.

# loading the data
final_results_with_pvalue <- read_csv("final_results_with_pvalue.csv")
View(final_results_with_pvalue)
blood_true_corr <- final_results_with_pvalue %>% filter(tissue=="Blood")
skin_true_corr <- final_results_with_pvalue %>% filter(tissue=="Skin")

## selecting only pvalues > 0.05
blood_true_sig <- blood_true_corr %>% filter(species == "Homo sapiens") %>% filter(p_value > 0.05)
View(blood_true_sig)
## Calculating the fdr in blood samples only
## FDR in blood_sample_N = 12
library(tibble)   # For rownames_to_column()
library(dplyr)    # For data manipulation

# Read metadata and get unique species names
blood_samples_12 <- read.csv("~/Phd_data/blood_samples_12.csv")
unique_species <- unique(blood_samples_12$species)

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
  matching_fractions <- data.frame(Replicate = colnames(sp_data_pval), False_positive = NA, percent = NA)
  
  # Loop over each column in sp_data_pval
  for (i in seq_along(colnames(sp_data_pval))) {
    
    # Extract the current column
    current_column <- colnames(sp_data_pval)[i]
    
    #  Filter for significant p-values in the current column
    filtered_pval <- sp_data_pval %>%
      select(current_column) %>%  # Select only the current column
      filter(!!sym(current_column) < 0.05) %>%  # Filter p < 0.05
      rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
    
    #  Perform an inner join with blood_true_corr to find matching CpG sites
    matching_rows <- inner_join(filtered_pval, blood_true_sig, by = "Beta_ID")
    
    # Count the number of matching row names
    num_matching <- nrow(matching_rows)
    print(num_matching)
    #  Calculate the fraction of matches
    total_rows <- nrow(sp_data_pval)
    print(total_rows)
    fraction_matching <- (num_matching / total_rows)
    percent <- (num_matching / total_rows) * 100
    # Store the result
    matching_fractions$False_positive[i] <- fraction_matching
    matching_fractions$percent[i] <- percent
  }
  
  # Print the final results
  print(matching_fractions)
  view(matching_fractions)
  # Save results to a CSV file
  write.csv(matching_fractions, "~/Phd_data/outputs/matching_fractions.csv", row.names = FALSE)
  matching_fractions$species <- "Homo sapiens"
  matching_fractions$sample_size <- "12"
  
  
  
  
  
########################### canis lupus familiaris N = 12 #################################  
 canis_true_sig <- blood_true_corr %>% filter(species == "Canis lupus familiaris") %>% filter(p_value > 0.05)
  View(canis_true_sig)
  
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
  canis_fractions <- data.frame(Replicate = colnames(sp_data_pval), False_positive = NA, percent = NA)
  
  # Loop over each column in sp_data_pval
  for (i in seq_along(colnames(sp_data_pval))) {
    
    # Extract the current column
    current_column <- colnames(sp_data_pval)[i]
    
    #  Filter for significant p-values in the current column
    filtered_pval <- sp_data_pval %>%
      select(current_column) %>%  # Select only the current column
      filter(!!sym(current_column) < 0.05) %>%  # Filter p < 0.05
      rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
    
    #  Perform an inner join with blood_true_corr to find matching CpG sites
    matching_rows <- inner_join(filtered_pval, canis_true_sig, by = "Beta_ID")
    
    # Count the number of matching row names
    num_matching <- nrow(matching_rows)
    print(num_matching)
    #  Calculate the fraction of matches
    total_rows <- nrow(sp_data_pval)
    print(total_rows)
    fraction_matching <- (num_matching / total_rows)
    percent <- (num_matching / total_rows) * 100
    # Store the result
    canis_fractions$False_positive[i] <- fraction_matching
    canis_fractions$percent[i] <- percent
  }
  
  # Print the final results
  print(canis_fractions)
  view(canis_fractions)
 canis_fractions$species <- "Canis lupus familiaris"
 canis_fractions$sample_size <- "12"
  # Save results to a CSV file
  write.csv(matching_fractions, "~/Phd_data/outputs/matching_fractions.csv", row.names = FALSE)

  
  # Initialize an empty list to store results for all species
  ############### Homo sapiens #########################
## load the data
  blood_sample_36_corr <- readRDS("~/Phd_data/outputs/blood_samples_36_correlation.rds")
  

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
  homo_36_fractions <- data.frame(Replicate = colnames(sp_data_pval), False_positive = NA, percent = NA)
  
  # Loop over each column in sp_data_pval
  for (i in seq_along(colnames(sp_data_pval))) {
    
    # Extract the current column
    current_column <- colnames(sp_data_pval)[i]
    
    #  Filter for significant p-values in the current column
    filtered_pval <- sp_data_pval %>%
      select(current_column) %>%  # Select only the current column
      filter(!!sym(current_column) < 0.05) %>%  # Filter p < 0.05
      rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
    
    #  Perform an inner join with blood_true_corr to find matching CpG sites
    matching_rows <- inner_join(filtered_pval, blood_true_sig, by = "Beta_ID")
    
    # Count the number of matching row names
    num_matching <- nrow(matching_rows)
    print(num_matching)
    #  Calculate the fraction of matches
    total_rows <- nrow(sp_data_pval)
    print(total_rows)
    fraction_matching <- (num_matching / total_rows)
    percent <- (num_matching / total_rows) * 100
    # Store the result
    homo_36_fractions$False_positive[i] <- fraction_matching
    homo_36_fractions$percent[i] <- percent
  }
  
  # Print the final results
  print(homo_36_fractions)
  view(homo_36_fractions)
  # Save results to a CSV file
  # write.csv(matching_fractions, "~/Phd_data/outputs/matching_fractions.csv", row.names = FALSE)
  homo_36_fractions$species <- "Homo sapiens"
  homo_36_fractions$sample_size <- "36"
  
############################## Canis lupus N =36 #################################
  canis_true_sig <- blood_true_corr %>% filter(species == "Canis lupus familiaris") %>% filter(p_value > 0.05)
  View(canis_true_sig)
  
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
  canis_36_fractions <- data.frame(Replicate = colnames(sp_data_pval), False_positive = NA, percent = NA)
  
  # Loop over each column in sp_data_pval
  for (i in seq_along(colnames(sp_data_pval))) {
    
    # Extract the current column
    current_column <- colnames(sp_data_pval)[i]
    
    #  Filter for significant p-values in the current column
    filtered_pval <- sp_data_pval %>%
      select(current_column) %>%  # Select only the current column
      filter(!!sym(current_column) < 0.05) %>%  # Filter p < 0.05
      rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
    
    #  Perform an inner join with blood_true_corr to find matching CpG sites
    matching_rows <- inner_join(filtered_pval, canis_true_sig, by = "Beta_ID")
    
    # Count the number of matching row names
    num_matching <- nrow(matching_rows)
    print(num_matching)
    #  Calculate the fraction of matches
    total_rows <- nrow(sp_data_pval)
    print(total_rows)
    fraction_matching <- (num_matching / total_rows)
    percent <- (num_matching / total_rows) * 100
    # Store the result
    canis_36_fractions$False_positive[i] <- fraction_matching
    canis_36_fractions$percent[i] <- percent
  }
  
  # Print the final results
  print(canis_36_fractions)
  view(canis_36_fractions)
  canis_36_fractions$species <- "Canis lupus familiaris"
  canis_36_fractions$sample_size <- "36"
  # Save results to a CSV file
  write.csv(matching_fractions, "~/Phd_data/outputs/matching_fractions.csv", row.names = FALSE)

  ############################## Canis lupus N = 60 #################################
  canis_true_sig <- blood_true_corr %>% filter(species == "Canis lupus familiaris") %>% filter(p_value > 0.05)
  View(canis_true_sig)
  
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
  canis_60_fractions <- data.frame(Replicate = colnames(sp_data_pval), False_positive = NA, percent = NA)
  
  # Loop over each column in sp_data_pval
  for (i in seq_along(colnames(sp_data_pval))) {
    
    # Extract the current column
    current_column <- colnames(sp_data_pval)[i]
    
    #  Filter for significant p-values in the current column
    filtered_pval <- sp_data_pval %>%
      select(current_column) %>%  # Select only the current column
      filter(!!sym(current_column) < 0.05) %>%  # Filter p < 0.05
      rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
    
    #  Perform an inner join with blood_true_corr to find matching CpG sites
    matching_rows <- inner_join(filtered_pval, canis_true_sig, by = "Beta_ID")
    
    # Count the number of matching row names
    num_matching <- nrow(matching_rows)
    print(num_matching)
    #  Calculate the fraction of matches
    total_rows <- nrow(sp_data_pval)
    print(total_rows)
    fraction_matching <- (num_matching / total_rows)
    percent <- (num_matching / total_rows) * 100
    # Store the result
    canis_60_fractions$False_positive[i] <- fraction_matching
    canis_60_fractions$percent[i] <- percent
  }
  
  # Print the final results
  print(canis_60_fractions)
  view(canis_60_fractions)
  canis_60_fractions$species <- "Canis lupus familiaris"
  canis_60_fractions$sample_size <- "60"
  # Save results to a CSV file
  write.csv(matching_fractions, "~/Phd_data/outputs/matching_fractions.csv", row.names = FALSE)  
  
  
  
######################## Homo sapiens N = 60 ###############################
  ## load the data
  blood_sample_60_corr <- readRDS("~/Phd_data/outputs/blood_samples_60_correlation.rds")
  
  
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
  homo_60_fractions <- data.frame(Replicate = colnames(sp_data_pval), False_positive = NA, percent = NA)
  
  # Loop over each column in sp_data_pval
  for (i in seq_along(colnames(sp_data_pval))) {
    
    # Extract the current column
    current_column <- colnames(sp_data_pval)[i]
    
    #  Filter for significant p-values in the current column
    filtered_pval <- sp_data_pval %>%
      select(current_column) %>%  # Select only the current column
      filter(!!sym(current_column) < 0.05) %>%  # Filter p < 0.05
      rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
    
    #  Perform an inner join with blood_true_corr to find matching CpG sites
    matching_rows <- inner_join(filtered_pval, blood_true_sig, by = "Beta_ID")
    
    # Count the number of matching row names
    num_matching <- nrow(matching_rows)
    print(num_matching)
    #  Calculate the fraction of matches
    total_rows <- nrow(sp_data_pval)
    print(total_rows)
    fraction_matching <- (num_matching / total_rows)
    percent <- (num_matching / total_rows) * 100
    # Store the result
    homo_60_fractions$False_positive[i] <- fraction_matching
    homo_60_fractions$percent[i] <- percent
  }
  
  # Print the final results
  print(homo_60_fractions)
  view(homo_60_fractions)
  # Save results to a CSV file
  # write.csv(matching_fractions, "~/Phd_data/outputs/matching_fractions.csv", row.names = FALSE)
  homo_60_fractions$species <- "Homo sapiens"
  homo_60_fractions$sample_size <- "60"

  
######################## Homo sapiens N = 90 ###############################
  ## load the data
  blood_sample_90_corr <- readRDS("~/Phd_data/outputs/blood_samples_90_correlation.rds")
  
  
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
  homo_90_fractions <- data.frame(Replicate = colnames(sp_data_pval), False_positive = NA, percent = NA)
  
  # Loop over each column in sp_data_pval
  for (i in seq_along(colnames(sp_data_pval))) {
    
    # Extract the current column
    current_column <- colnames(sp_data_pval)[i]
    
    #  Filter for significant p-values in the current column
    filtered_pval <- sp_data_pval %>%
      select(current_column) %>%  # Select only the current column
      filter(!!sym(current_column) < 0.05) %>%  # Filter p < 0.05
      rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
    
    #  Perform an inner join with blood_true_corr to find matching CpG sites
    matching_rows <- inner_join(filtered_pval, blood_true_sig, by = "Beta_ID")
    
    # Count the number of matching row names
    num_matching <- nrow(matching_rows)
    print(num_matching)
    #  Calculate the fraction of matches
    total_rows <- nrow(sp_data_pval)
    print(total_rows)
    fraction_matching <- (num_matching / total_rows)
    percent <- (num_matching / total_rows) * 100
    # Store the result
    homo_90_fractions$False_positive[i] <- fraction_matching
    homo_90_fractions$percent[i] <- percent
  }
  
  # Print the final results
  print(homo_90_fractions)
  view(homo_90_fractions)
  # Save results to a CSV file
  # write.csv(matching_fractions, "~/Phd_data/outputs/matching_fractions.csv", row.names = FALSE)
  homo_90_fractions$species <- "Homo sapiens"
  homo_90_fractions$sample_size <- "90"    
  
  
  
  ############################## Canis lupus N = 90 #################################
  canis_true_sig <- blood_true_corr %>% filter(species == "Canis lupus familiaris") %>% filter(p_value > 0.05)
  View(canis_true_sig)
  
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
  canis_90_fractions <- data.frame(Replicate = colnames(sp_data_pval), False_positive = NA, percent = NA)
  
  # Loop over each column in sp_data_pval
  for (i in seq_along(colnames(sp_data_pval))) {
    
    # Extract the current column
    current_column <- colnames(sp_data_pval)[i]
    
    #  Filter for significant p-values in the current column
    filtered_pval <- sp_data_pval %>%
      select(current_column) %>%  # Select only the current column
      filter(!!sym(current_column) < 0.05) %>%  # Filter p < 0.05
      rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
    
    #  Perform an inner join with blood_true_corr to find matching CpG sites
    matching_rows <- inner_join(filtered_pval, canis_true_sig, by = "Beta_ID")
    
    # Count the number of matching row names
    num_matching <- nrow(matching_rows)
    print(num_matching)
    #  Calculate the fraction of matches
    total_rows <- nrow(sp_data_pval)
    print(total_rows)
    fraction_matching <- (num_matching / total_rows)
    percent <- (num_matching / total_rows) * 100
    # Store the result
    canis_90_fractions$False_positive[i] <- fraction_matching
    canis_90_fractions$percent[i] <- percent
  }
  
  # Print the final results
  print(canis_90_fractions)
  view(canis_90_fractions)
  canis_90_fractions$species <- "Canis lupus familiaris"
  canis_90_fractions$sample_size <- "90"
  # Save results to a CSV file
  write.csv(matching_fractions, "~/Phd_data/outputs/matching_fractions.csv", row.names = FALSE)  
  
######### combine all the false positive rate across the samples sizes in homo_sapiens #############
## homo sapiens
combined_fdr <- rbind(matching_fractions,homo_36_fractions, homo_60_fractions, homo_90_fractions)  
View(combined_fdr)  

## violin plot  of FDR 
ggplot(combined_fdr, aes(x =factor(sample_size), y=percent)) + geom_boxplot() +
   geom_point(size = 2, alpha= 0.8) +
  labs(title = "False positive rate DNA methylation from Homo sapiens blood sample",
       x = "Sample size", y = "False positive rate (%)") +
  theme_minimal() + theme(panel.background = element_rect(color = "black", # Color of the border
                                                          size = 0.2))

ggplot(combined_fdr, aes(x=percent,  fill=sample_size)) +
    geom_density(alpha=0.2) + theme_bw() + scale_x_continuous(breaks = seq(0, 5, by = 0.5)) +
  scale_y_continuous(breaks = seq(0, 2, by = 0.5)) +
  theme(panel.background = element_rect(color = "black", # Color of the border
                                        size = 0.2)) + 
  labs(x = "false positive rate (in %)", 
       title = "The false positive rate using downsampled DNAm blood samples from Homo sapiens")

## with facet wrap
ggplot(combined_fdr, aes(x=percent,  fill= sample_size)) +
  geom_density(alpha=0.2) + theme_bw() + scale_x_continuous(breaks = seq(0, 5, by = 0.5)) +
  scale_y_continuous(breaks = seq(0, 2, by = 0.5)) +
  theme(panel.background = element_rect(color = "black", # Color of the border
                                        size = 0.2)) + facet_wrap(~sample_size) + #theme_tufte()
  labs(x = "false positive rate (in %)", 
       title = "The false positive rate using downsampled DNAm blood samples from Homo sapiens")




## canis
canis_combined_fdr <- rbind(canis_fractions, canis_36_fractions, canis_60_fractions, canis_90_fractions)
view(canis_combined_fdr)

ggplot(canis_combined_fdr, aes(x=percent,  fill=sample_size)) +
  geom_density(alpha=0.2) + theme_bw() + scale_x_continuous(breaks = seq(0, 6, by = 0.5)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme(panel.background = element_rect(color = "black", # Color of the border
                                        size = 0.2)) + 
  labs(x = "false positive rate (in %)", 
  title = "The false positive rate using downsampled DNAm blood samples from Canis lupus familiaris")


ggplot(canis_combined_fdr, aes(x=percent,  fill=sample_size)) +
  geom_density(alpha=0.2) + theme_bw() + scale_x_continuous(breaks = seq(0, 6, by = 0.5)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme(panel.background = element_rect(color = "black", # Color of the border
                                        size = 0.2)) + facet_wrap(~ sample_size) +
  labs(x = "false positive rate (in %)", 
       title = "The false positive rate using downsampled DNAm blood samples from Canis lupus familiaris")


################ N =32 
# Read the correlation results (assuming this file contains the results for all species)
skin_sample_32_corr <- readRDS("~/Phd_data/outputs/skin_sample_32_correlation.rds")

# Initialize an empty list to store results for all species
############### Homo sapiens #########################

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
homo_skin_fdr32 <- data.frame(Replicate = colnames(sk_data_pval), False_positive = NA, percent = NA)

# Loop over each column in sp_data_pval
for (i in seq_along(colnames(sk_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sk_data_pval)[i]
  print(current_column)
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sk_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) < 0.05) %>%  # Filter p < 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with blood_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, skin_true_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  
  #  Calculate the fraction of matches
  total_rows <- nrow(sk_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  
  # Store the result
  homo_skin_fdr32$False_positive[i] <- fraction_matching
  homo_skin_fdr32$percent[i] <- percent
  
  cat("calculation is completed... ")
}

# Print the final results
print(homo_skin_fdr32)
view(homo_skin_fdr32)
# Save results to a CSV file
homo_skin_fdr32$species <- "Homo sapiens"
homo_skin_fdr32$sample_size <- "32"
homo_skin_fdr32$tissue <- "skin"

############# N =50

# Read the correlation results (assuming this file contains the results for all species)
skin_sample_50_corr <- readRDS("~/Phd_data/outputs/skin_sample_50_correlation.rds")

# Initialize an empty list to store results for all species
############### Homo sapiens #########################

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
homo_skin_fdr50 <- data.frame(Replicate = colnames(sk_data_pval), False_positive = NA, percent = NA)

# Loop over each column in sp_data_pval
for (i in seq_along(colnames(sk_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sk_data_pval)[i]
  print(current_column)
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sk_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) < 0.05) %>%  # Filter p < 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with blood_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, skin_true_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  
  #  Calculate the fraction of matches
  total_rows <- nrow(sk_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  
  # Store the result
  homo_skin_fdr50$False_positive[i] <- fraction_matching
  homo_skin_fdr50$percent[i] <- percent
  
  cat("calculation is completed... ")
}

# Print the final results
print(homo_skin_fdr50)
view(homo_skin_fdr50)
# Save results to a CSV file
homo_skin_fdr50$species <- "Homo sapiens"
homo_skin_fdr50$sample_size <- "50"
homo_skin_fdr50$tissue <- "skin"


##################### N = 70 ###################
# Read the correlation results (assuming this file contains the results for all species)
skin_sample_70_corr <- readRDS("~/Phd_data/outputs/skin_sample_70_correlation.rds")

# Initialize an empty list to store results for all species
############### Homo sapiens #########################

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
homo_skin_fdr70 <- data.frame(Replicate = colnames(sk_data_pval), False_positive = NA, percent = NA)

# Loop over each column in sp_data_pval
for (i in seq_along(colnames(sk_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sk_data_pval)[i]
  print(current_column)
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sk_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) < 0.05) %>%  # Filter p < 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with blood_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, skin_true_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  
  #  Calculate the fraction of matches
  total_rows <- nrow(sk_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  
  # Store the result
  homo_skin_fdr70$False_positive[i] <- fraction_matching
  homo_skin_fdr70$percent[i] <- percent
  
  cat("calculation is completed... ")
}

# Print the final results
print(homo_skin_fdr70)
view(homo_skin_fdr70)
# Save results to a CSV file
homo_skin_fdr70$species <- "Homo sapiens"
homo_skin_fdr70$sample_size <- "70"
homo_skin_fdr70$tissue <- "skin"

homo_skin_combined_fdr <- rbind(homo_skin_fdr14, homo_skin_fdr32, homo_skin_fdr50, homo_skin_fdr70)
view(homo_skin_combined_fdr)

ggplot(homo_skin_combined_fdr, aes(x=percent,  fill=sample_size)) +
  geom_density(alpha=0.2) + theme_bw() + scale_x_continuous(breaks = seq(0, 11, by = 0.5)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme(panel.background = element_rect(color = "black", # Color of the border
                                        size = 0.2)) + 
  labs(x = "false positive rate (in %)", 
       title = "The false positive rate using downsampled DNAm skin samples from Homo sapiens")


ggplot(homo_skin_combined_fdr, aes(x=percent,  fill=sample_size)) +
  geom_density(alpha=0.2) + theme_bw() + scale_x_continuous(breaks = seq(0, 11, by = 2)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme(panel.background = element_rect(color = "black", # Color of the border
                                        size = 0.2)) + facet_wrap(~ sample_size) +
  labs(x = "false positive rate (in %)", 
       title = "The false positive rate using downsampled DNAm skin samples from Homo sapiens")









########################### Tursiops aduncus skin samples N = 12 #################################  
Tursiops_true_sig <- skin_true_corr %>% filter(species == "Tursiops aduncus") %>% filter(p_value > 0.05)
View(Tursiops_true_sig)

# Extract correlation and p-value data for the species
rownames(skin_sample_14_corr[["Tursiops aduncus"]][["correlations"]]) <- rownames(skin_metadata_betas)
sp_data_corr <- as.data.frame(skin_sample_14_corr[["Tursiops aduncus"]][["correlations"]])
rownames(skin_sample_14_corr[["Tursiops aduncus"]][["p_values"]]) <- rownames(skin_metadata_betas)
sp_data_pval <- as.data.frame(skin_sample_14_corr[["Tursiops aduncus"]][["p_values"]])

library(dplyr)
library(tibble)

# Extract p-values data frame for the current species
sk_data_pval <- as.data.frame(skin_sample_14_corr[["Tursiops aduncus"]][["p_values"]])

# Ensure skin_true_corr rownames are a column for joining
#skin_true_sig <- skin_true_sig %>%
# rownames_to_column(var = "Beta_ID")

# Initialize a result storage
Tursiops_fractions <- data.frame(Replicate = colnames(sp_data_pval), False_positive = NA, percent = NA)

# Loop over each column in sp_data_pval
for (i in seq_along(colnames(sk_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sk_data_pval)[i]
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sk_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) < 0.05) %>%  # Filter p < 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with skin_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, Tursiops_true_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  #  Calculate the fraction of matches
  total_rows <- nrow(sk_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  # Store the result
  Tursiops_fractions$False_positive[i] <- fraction_matching
  Tursiops_fractions$percent[i] <- percent
}

# Print the final results
print(Tursiops_fractions)
view(Tursiops_fractions)
Tursiops_fractions$species <- "Tursiops aduncus"
Tursiops_fractions$sample_size <- "14"
Tursiops_fractions$tissue <- "skin"



##################### Tursiops aduncus skin samples N = 32 ############################
# Extract correlation and p-value data for the skecies
rownames(skin_sample_32_corr[["Tursiops aduncus"]][["correlations"]]) <- rownames(skin_metadata_betas)
sk_data_corr <- as.data.frame(skin_sample_32_corr[["Tursiops aduncus"]][["correlations"]])
rownames(skin_sample_32_corr[["Tursiops aduncus"]][["p_values"]]) <- rownames(skin_metadata_betas)
sk_data_pval <- as.data.frame(skin_sample_32_corr[["Tursiops aduncus"]][["p_values"]])

library(dplyr)
library(tibble)

# Extract p-values data frame for the current skecies
sk_data_pval <- as.data.frame(skin_sample_32_corr[["Tursiops aduncus"]][["p_values"]])

# Ensure skin_true_corr rownames are a column for joining
#skin_true_sig <- skin_true_sig %>%
# rownames_to_column(var = "Beta_ID")

# Initialize a result storage
Tursiops_fraction_32 <- data.frame(Replicate = colnames(sk_data_pval), False_positive = NA, percent = NA)

# Loop over each column in sk_data_pval
for (i in seq_along(colnames(sk_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sk_data_pval)[i]
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sk_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) < 0.05) %>%  # Filter p < 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with skin_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, Tursiops_true_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  #  Calculate the fraction of matches
  total_rows <- nrow(sk_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  # Store the result
  Tursiops_fraction_32$False_positive[i] <- fraction_matching
  Tursiops_fraction_32$percent[i] <- percent
}

# Print the final results
print(Tursiops_fraction_32)
view(Tursiops_fraction_32)
Tursiops_fraction_32$skecies <- "Tursiops aduncus"
Tursiops_fraction_32$sample_size <- "32"
Tursiops_fraction_32$tissue <- "skin"


####################### Tursiops aduncus skin samples N = 50 ################################
# Extract correlation and p-value data for the skecies
rownames(skin_sample_50_corr[["Tursiops aduncus"]][["correlations"]]) <- rownames(skin_metadata_betas)
sk_data_corr <- as.data.frame(skin_sample_50_corr[["Tursiops aduncus"]][["correlations"]])
rownames(skin_sample_50_corr[["Tursiops aduncus"]][["p_values"]]) <- rownames(skin_metadata_betas)
sk_data_pval <- as.data.frame(skin_sample_50_corr[["Tursiops aduncus"]][["p_values"]])

library(dplyr)
library(tibble)

# Extract p-values data frame for the current skecies
sk_data_pval <- as.data.frame(skin_sample_50_corr[["Tursiops aduncus"]][["p_values"]])

# Ensure skin_true_corr rownames are a column for joining
#skin_true_sig <- skin_true_sig %>%
# rownames_to_column(var = "Beta_ID")

# Initialize a result storage
Tursiops_fraction_50 <- data.frame(Replicate = colnames(sk_data_pval), False_positive = NA, percent = NA)

# Loop over each column in sk_data_pval
for (i in seq_along(colnames(sk_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sk_data_pval)[i]
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sk_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) < 0.05) %>%  # Filter p < 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with skin_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, Tursiops_true_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  #  Calculate the fraction of matches
  total_rows <- nrow(sk_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  # Store the result
  Tursiops_fraction_50$False_positive[i] <- fraction_matching
  Tursiops_fraction_50$percent[i] <- percent
}

# Print the final results
print(Tursiops_fraction_50)
view(Tursiops_fraction_50)
Tursiops_fraction_50$skecies <- "Tursiops aduncus"
Tursiops_fraction_50$sample_size <- "50"
Tursiops_fraction_50$tissue <- "skin"

######################### Trusiops aduncus N = 70  #############################

# Extract correlation and p-value data for the skecies
rownames(skin_sample_70_corr[["Tursiops aduncus"]][["correlations"]]) <- rownames(skin_metadata_betas)
sk_data_corr <- as.data.frame(skin_sample_70_corr[["Tursiops aduncus"]][["correlations"]])
rownames(skin_sample_70_corr[["Tursiops aduncus"]][["p_values"]]) <- rownames(skin_metadata_betas)
sk_data_pval <- as.data.frame(skin_sample_70_corr[["Tursiops aduncus"]][["p_values"]])

library(dplyr)
library(tibble)

# Extract p-values data frame for the current skecies
sk_data_pval <- as.data.frame(skin_sample_70_corr[["Tursiops aduncus"]][["p_values"]])

# Ensure skin_true_corr rownames are a column for joining
#skin_true_sig <- skin_true_sig %>%
# rownames_to_column(var = "Beta_ID")

# Initialize a result storage
Tursiops_fraction_70 <- data.frame(Replicate = colnames(sk_data_pval), False_positive = NA, percent = NA)

# Loop over each column in sk_data_pval
for (i in seq_along(colnames(sk_data_pval))) {
  
  # Extract the current column
  current_column <- colnames(sk_data_pval)[i]
  
  #  Filter for significant p-values in the current column
  filtered_pval <- sk_data_pval %>%
    select(current_column) %>%  # Select only the current column
    filter(!!sym(current_column) < 0.05) %>%  # Filter p < 0.05
    rownames_to_column(var = "Beta_ID")  # Convert rownames to a column
  
  #  Perform an inner join with skin_true_corr to find matching CpG sites
  matching_rows <- inner_join(filtered_pval, Tursiops_true_sig, by = "Beta_ID")
  
  # Count the number of matching row names
  num_matching <- nrow(matching_rows)
  print(num_matching)
  #  Calculate the fraction of matches
  total_rows <- nrow(sk_data_pval)
  print(total_rows)
  fraction_matching <- (num_matching / total_rows)
  percent <- (num_matching / total_rows) * 100
  # Store the result
  Tursiops_fraction_70$False_positive[i] <- fraction_matching
  Tursiops_fraction_70$percent[i] <- percent
}

# Print the final results
print(Tursiops_fraction_70)
view(Tursiops_fraction_70)
Tursiops_fraction_70$skecies <- "Tursiops aduncus"
Tursiops_fraction_70$sample_size <- "70"
Tursiops_fraction_70$tissue <- "skin"

## Tursiops plot
tursiops_skin_combined_fdr <- rbind(Tursiops_fractions, Tursiops_fraction_32, Tursiops_fraction_50, Tursiops_fraction_70)
view(tursiops_skin_combined_fdr)

ggplot(tursiops_skin_combined_fdr, aes(x=percent,  fill=sample_size)) +
  geom_density(alpha=0.2) + theme_bw() + scale_x_continuous(breaks = seq(0, 20, by = 2)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme(panel.background = element_rect(color = "black", # Color of the border
                                        size = 0.2)) + 
  labs(x = "False positive rate (in %)", 
       title = "The false positive rate using downsampled DNAm skin samples from Tursiops aduncus")


ggplot(tursiops_skin_combined_fdr, aes(x=percent,  fill=sample_size)) +
  geom_density(alpha=0.2) + theme_bw() + scale_x_continuous(breaks = seq(0, 20, by = 2)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme(panel.background = element_rect(color = "black", # Color of the border
                                        size = 0.2)) + facet_wrap(~ sample_size) +
  labs(x = "false positive rate (in %)", 
       title = "The false positive rate using downsampled DNAm skin samples from Tursiops aduncus")















