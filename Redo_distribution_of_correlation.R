######### Author : Adeyemi Akinade
######### Place : Clemson University, SC
######### lab : Gopalan's lab
######### Date : 11 February, 2025
######### Project : Comparative Epigenetics of Ageing


############################ loading data with 36 samples in blood samples ################################### 
blood_samples_36_correlation <- readRDS("~/Phd_data/outputs/blood_samples_36_correlation.rds")


## adding the rownames of the blood_correlation of Homo sapiens
rownames(blood_samples_36_correlation[["Homo sapiens"]][["correlations"]]) <- rownames(blood_merged_betas)

Homo_sapiens_corr36 <- (blood_samples_36_correlation[["Homo sapiens"]][["correlations"]])
View(Homo_sapiens_corr36)

##make the data a dataframe
Homo_sapiens_corr36 <- as.data.frame(Homo_sapiens_corr36)

# Calculate the row means
Homo_sapiens_corr36$Mean <- apply(Homo_sapiens_corr36, 1, mean, na.rm = TRUE)
Homo_sapiens_corr36$Min <- apply(Homo_sapiens_corr36, 1, min, na.rm = TRUE)
Homo_sapiens_corr36$Max <- apply(Homo_sapiens_corr36, 1, max, na.rm = TRUE)


# Calculate the row standard deviations
Homo_sapiens_corr36$std <- apply(Homo_sapiens_corr36, 1, sd, na.rm = TRUE)

# If you want to create columns for (mean - 1 SD) and (mean + 1 SD)
Homo_sapiens_corr36$lower_bound_std1 <- Homo_sapiens_corr36$Mean - Homo_sapiens_corr36$std
Homo_sapiens_corr36$upper_bound_std1 <- Homo_sapiens_corr36$Mean + Homo_sapiens_corr36$std

Homo_sapiens_corr36$lower_bound_std2 <- Homo_sapiens_corr36$Mean - (2*Homo_sapiens_corr36$std)
Homo_sapiens_corr36$upper_bound_std2 <- Homo_sapiens_corr36$Mean + (2*Homo_sapiens_corr36$std)


# Optionally, you can view the first few rows of these new columns
head(Homo_sapiens_corr36[, c("Mean", "std", "lower_bound_std1", "lower_bound_std2","upper_bound_std1", "upper_bound_std2")])

## Loading the estimated correlation data
estimated_data <- readRDS("~/Phd_data/estimated_correlation.RDS")

# checking the estimated correlation data
head(estimated_data)
estimated_data <- as.data.frame(estimated_data)

#changing the column name from Beta_ID to cpg_sites
estimated_data <- estimated_data %>% rename(cpg_site = Beta_ID )
head(estimated_data)

## selecting estimated correlation
sapiens_corr36 <- Homo_sapiens_corr36 %>% select(Mean, Min, Max, std, lower_bound_std1, lower_bound_std2,upper_bound_std1, upper_bound_std2)
est_corr <- estimated_data %>% select(cpg_site, correlation)

all_sapiens_36 <- cbind(sapiens_corr36,est_corr)
View(all_sapiens_36)
all_sapiens_36 <- all_sapiens_36 %>% arrange(desc(Mean)) 

plot(all_sapiens_36$Mean, all_sapiens_36$correlation)

## Plotting the points on estimated correlation and replicated correlations
ggplot(all_sapiens_36, aes(x = correlation, y = Mean)) +
  geom_line(color= "blue") #+ geom_hex(bins = 40)

##################### Ovis aries N = 36 #############################  

blood_samples_36_correlation <- readRDS("~/Phd_data/outputs/blood_samples_36_correlation.rds")


## adding the rownames of the blood_correlation of Homo sapiens
rownames(blood_samples_36_correlation[["Ovis aries"]][["correlations"]]) <- rownames(blood_merged_betas)

ovis_corr36 <- (blood_samples_36_correlation[["Ovis aries"]][["correlations"]])
head(ovis_corr36)

##make the data a dataframe
ovis_corr36 <- as.data.frame(ovis_corr36)

# Calculate the row means
ovis_corr36$Mean <- apply(ovis_corr36, 1, mean, na.rm = TRUE)
ovis_corr36$Min <- apply(ovis_corr36, 1, min, na.rm = TRUE)
ovis_corr36$Max <- apply(ovis_corr36, 1, max, na.rm = TRUE)
ovis_corr36$Median <- apply(ovis_corr36, 1, median, na.rm = TRUE)

# Calculate the row standard deviations
ovis_corr36$std <- apply(ovis_corr36, 1, sd, na.rm = TRUE)

# If you want to create columns for (mean - 1 SD) and (mean + 1 SD)
ovis_corr36$lower_bound_std1 <- ovis_corr36$Mean - ovis_corr36$std
ovis_corr36$upper_bound_std1 <- ovis_corr36$Mean + ovis_corr36$std



# Optionally, you can view the first few rows of these new columns
head(ovis_corr36[, c("Mean","Median", "std", "lower_bound_std1","upper_bound_std1")])

## Loading the estimated correlation data
estimated_data <- readRDS("~/Phd_data/estimated_correlation.RDS")

# checking the estimated correlation data
head(estimated_data)
estimated_data <- as.data.frame(estimated_data)

#changing the column name from Beta_ID to cpg_sites
estimated_data <- estimated_data %>% rename(cpg_site = Beta_ID )
head(estimated_data)

## selecting estimated correlation
ovis_sel_corr36 <- ovis_corr36 %>% select(Mean, Median, Min, Max, std, lower_bound_std1,upper_bound_std1)
est_corr <- estimated_data %>% select(cpg_site, correlation)

all_ovis_36 <- cbind(ovis_sel_corr36,est_corr)
View(all_ovis_36)
all_ovis_36 <- all_sapiens_36 %>% arrange(desc(Mean)) 

plot(all_ovis_36$Mean, all_ovis_36$correlation, type = "l", lwd = 2, col = "black")

## Plotting the points on estimated correlation and replicated correlations
ggplot(all_sapiens_36, aes(x = correlation, y = Mean)) +
  geom_line(color= "blue") #+ geom_hex(bins = 40)


############ using pivot longer to see the disctribution of the replicated correlations ##############
Homo_sapiens_corr <- as.data.frame(Homo_sapiens_corr)
ovis_corr_36 <- rownames_to_column(ovis_corr36, var = "cpg_site")
head(ovis_corr_36)

# Convert the 50 replicate columns into long format
ovis_long <- ovis_corr_36 %>%
  pivot_longer(cols = starts_with("Rep_"), names_to = "replicate", values_to = "corr_value")

head(ovis_long)

# Merge with estimated correlation data
ovis_plot_data <- ovis_long %>%
  left_join(est_corr, by = "cpg_site")
head(ovis_plot_data)















###################### skin only ###################################
blood_samples_36_correlation <- readRDS("~/Phd_data/outputs/blood_samples_36_correlation.rds")


## adding the rownames of the blood_correlation of Homo sapiens
rownames(blood_samples_36_correlation[["Homo sapiens"]][["correlations"]]) <- rownames(blood_merged_betas)

Homo_sapiens_corr36 <- (blood_samples_36_correlation[["Homo sapiens"]][["correlations"]])
View(Homo_sapiens_corr36)

##make the data a dataframe
Homo_sapiens_corr36 <- as.data.frame(Homo_sapiens_corr36)

# Calculate the row means
Homo_sapiens_corr36$Mean <- apply(Homo_sapiens_corr36, 1, mean, na.rm = TRUE)
Homo_sapiens_corr36$Min <- apply(Homo_sapiens_corr36, 1, min, na.rm = TRUE)
Homo_sapiens_corr36$Max <- apply(Homo_sapiens_corr36, 1, max, na.rm = TRUE)


# Calculate the row standard deviations
Homo_sapiens_corr36$std <- apply(Homo_sapiens_corr36, 1, sd, na.rm = TRUE)

# If you want to create columns for (mean - 1 SD) and (mean + 1 SD)
Homo_sapiens_corr36$lower_bound_std1 <- Homo_sapiens_corr36$Mean - Homo_sapiens_corr36$std
Homo_sapiens_corr36$upper_bound_std1 <- Homo_sapiens_corr36$Mean + Homo_sapiens_corr36$std

Homo_sapiens_corr36$lower_bound_std2 <- Homo_sapiens_corr36$Mean - (2*Homo_sapiens_corr36$std)
Homo_sapiens_corr36$upper_bound_std2 <- Homo_sapiens_corr36$Mean + (2*Homo_sapiens_corr36$std)


# Optionally, you can view the first few rows of these new columns
head(Homo_sapiens_corr36[, c("Mean", "std", "lower_bound_std1", "lower_bound_std2","upper_bound_std1", "upper_bound_std2")])

## Loading the estimated correlation data
estimated_data <- readRDS("~/Phd_data/estimated_correlation.RDS")

# checking the estimated correlation data
head(estimated_data)
estimated_data <- as.data.frame(estimated_data)

#changing the column name from Beta_ID to cpg_sites
estimated_data <- estimated_data %>% rename(cpg_site = Beta_ID )
head(estimated_data)

## selecting estimated correlation
sapiens_corr36 <- Homo_sapiens_corr36 %>% select(Mean, Min, Max, std, lower_bound_std1, lower_bound_std2,upper_bound_std1, upper_bound_std2)
est_corr <- estimated_data %>% select(cpg_site, correlation)

all_sapiens_36 <- cbind(sapiens_corr36,est_corr)
View(all_sapiens_36)
all_sapiens_36 <- all_sapiens_36 %>% arrange(desc(Mean)) 

plot(all_sapiens_36$Mean, all_sapiens_36$correlation)

## Plotting the points on estimated correlation and replicated correlations
ggplot(all_sapiens_36, aes(x = correlation, y = Mean)) +
  geom_line(color= "blue") #+ geom_hex(bins = 40)
