######### Author : Adeyemi Akinade
######### Place : Clemson University, SC
######### lab : Gopalan's lab
######### Date : 11 February, 2025
######### Project : Comparative Epigenetics of Ageing


### loading data 
blood_samples_12_correlation <- readRDS("~/Phd_data/outputs/blood_samples_12_correlation.rds")

## adding the rownames of the blood_correlation of Homo sapiens
rownames(blood_samples_12_correlation[["Homo sapiens"]][["correlations"]]) <- rownames(blood_merged_betas)

Homo_sapiens_corr <- (blood_samples_12_correlation[["Homo sapiens"]][["correlations"]])
str(Homo_sapiens_corr)

## Extarcting the pvalue
rownames(blood_samples_12_correlation[["Homo sapiens"]][["p_values"]]) <- rownames(blood_merged_betas)
Homo_sapiens_p_values <- blood_samples_12_correlation[["Homo sapiens"]][["p_values"]]
View(Homo_sapiens_p_values)

## converting the cpgsites data structure to dataframe
homo_com_site <-  data.frame (Homo_sapiens_corr["cg16867657", ])
colnames(homo_com_site) <- "corr"
homo_com_site$corr <- as.numeric(homo_com_site$corr)
View(homo_com_site)

## converting the pvalues data structure to dataframe
homo_com_pval <-  data.frame (Homo_sapiens_p_values["cg16867657", ])
colnames(homo_com_pval) <- "pval"
homo_com_pval$pval <- as.numeric(homo_com_pval$pval)
View(homo_com_pval)
homo_com_site <- cbind(homo_com_site, homo_com_pval)



### histogram
 ggplot(homo_com_site, aes(x = corr )) + 
  geom_histogram(
    colour = 1, fill = "white") +  scale_x_continuous(breaks = seq(0.80, 1.0, by = 0.006)) +
  labs(title = " Age correlation with cg16867657 from 12 blood samples in humans", x = "correlation", y = "Count") +
  theme_minimal()
 
 ggplot(homo_com_pval, aes(x = pval )) + 
   geom_histogram(
     colour = 1, fill = "white") +
   labs(title = "pvalues of age correlation with cg16867657 using 12 blood samples from humans", x = "correlation", y = "Count") +
   theme_minimal()
 
 
 
 ############### Ovies aries   #################
 ## adding the rownames of the blood_correlation of Ovies aries
 rownames(blood_samples_12_correlation[["Ovis aries"]][["correlations"]]) <- rownames(blood_merged_betas)
 
 ovis_corr <- (blood_samples_12_correlation[["Ovis aries"]][["correlations"]])
 View(ovis_corr)
 
 ## Extarcting the pvalue
 rownames(blood_samples_12_correlation[["Ovis aries"]][["p_values"]]) <- rownames(blood_merged_betas)
ovis_p_values <- blood_samples_12_correlation[["Ovis aries"]][["p_values"]]
 View(ovis_p_values)
 
 ## converting the cpgsites data structure to dataframe
ovis_com_site <-  data.frame (ovis_corr["cg16867657", ])
 colnames(ovis_com_site) <- "corr"
 ovis_com_site$corr <- as.numeric(ovis_com_site$corr)
 View(ovis_com_site)
 
 ## converting the pvalues data structure to dataframe
 ovis_com_pval <-  data.frame (ovis_p_values["cg16867657", ])
 colnames(ovis_com_pval) <- "pval"
 ovis_com_pval$pval <- as.numeric(ovis_com_pval$pval)
 View(ovis_com_pval)
 ovis_com_site <- cbind(ovis_com_site, ovis_com_pval)
 
 
 
 ### histogram
 ggplot(ovis_com_site, aes(x = corr )) + 
   geom_histogram(
     colour = 1, fill = "white") +
   labs(title = "correlation of age with cg16867657 from 12 blood samples in Ovis aries", x = "correlation", y = "Count") +
   theme_minimal()
 
 ggplot(ovis_com_pval, aes(x = pval )) + 
   geom_histogram(
     colour = 1, fill = "white") +
   labs(title = "pvalues of correlation of age with cg16867657 using 12 blood samples in Ovis aries", x = "correlation", y = "Count") +
   theme_minimal()

  
 ############################ loading data with 36 samples ################################### 
 blood_samples_36_correlation <- readRDS("~/Phd_data/outputs/blood_samples_36_correlation.rds")
 
 ## adding the rownames of the blood_correlation of Homo sapiens
 rownames(blood_samples_36_correlation[["Homo sapiens"]][["correlations"]]) <- rownames(blood_merged_betas)
 
 Homo_sapiens_corr36 <- (blood_samples_36_correlation[["Homo sapiens"]][["correlations"]])
 View(Homo_sapiens_corr36)
 
 ## Extarcting the pvalue
 rownames(blood_samples_36_correlation[["Homo sapiens"]][["p_values"]]) <- rownames(blood_merged_betas)
 Homo_sapiens_p_values36 <- blood_samples_36_correlation[["Homo sapiens"]][["p_values"]]
 View(Homo_sapiens_p_values36)
 
 ## converting the cpgsites data structure to dataframe
 homo_com_site36 <-  data.frame (Homo_sapiens_corr36["cg16867657", ])
 colnames(homo_com_site36) <- "corr"
 homo_com_site36$corr <- as.numeric(homo_com_site36$corr)
 View(homo_com_site36)
 
 ## converting the pvalues data structure to dataframe
 homo_com_pval36 <-  data.frame (Homo_sapiens_p_values36["cg16867657", ])
 colnames(homo_com_pval36) <- "pval"
 homo_com_pval36$pval <- as.numeric(homo_com_pval36$pval)
 View(homo_com_pval36)
 homo_com_site36 <- cbind(homo_com_site36, homo_com_pval36)
 
 
 
 ### histogram
 ggplot(homo_com_site36, aes(x = corr )) + 
   geom_histogram(
     colour = 1, fill = "white") + scale_x_continuous(breaks = seq(0.80, 1.0, by = 0.004)) +
   labs(title = " Age correlation with cg16867657 from 36 blood samples in humans", x = "correlation", y = "Count") +
   theme_minimal()
 
 ggplot(homo_com_pval36, aes(x = pval )) + 
   geom_histogram(
     colour = 1, fill = "white") +
   labs(title = "pvalues of age correlation with cg16867657 using 36 blood samples from humans", x = "correlation", y = "Count") +
   theme_minimal()
  
 
 

############################ loading data with 60 samples ################################### 
 blood_samples_60_correlation <- readRDS("~/Phd_data/outputs/blood_samples_60_correlation.rds")
 
 ## adding the rownames of the blood_correlation of Homo sapiens
 rownames(blood_samples_60_correlation[["Homo sapiens"]][["correlations"]]) <- rownames(blood_merged_betas)
 
 Homo_sapiens_corr60 <- (blood_samples_60_correlation[["Homo sapiens"]][["correlations"]])
 View(Homo_sapiens_corr60)
 
 ## Extarcting the pvalue
 rownames(blood_samples_60_correlation[["Homo sapiens"]][["p_values"]]) <- rownames(blood_merged_betas)
 Homo_sapiens_p_values60 <- blood_samples_60_correlation[["Homo sapiens"]][["p_values"]]
 View(Homo_sapiens_p_values60)
 
 ## converting the cpgsites data structure to dataframe
 homo_com_site60 <-  data.frame (Homo_sapiens_corr60["cg16867657", ])
 colnames(homo_com_site60) <- "corr"
 homo_com_site60$corr <- as.numeric(homo_com_site60$corr)
 View(homo_com_site60)
 
 ## converting the pvalues data structure to dataframe
 homo_com_pval60 <-  data.frame (Homo_sapiens_p_values60["cg16867657", ])
 colnames(homo_com_pval60) <- "pval"
 homo_com_pval60$pval <- as.numeric(homo_com_pval60$pval)
 View(homo_com_pval60)
 homo_com_site60 <- cbind(homo_com_site60, homo_com_pval60)
 
 
 
 ### histogram
 ggplot(homo_com_site60, aes(x = corr )) + 
   geom_histogram(
     colour = 1, fill = "white") +   scale_x_continuous(breaks = seq(0.90, 1.0, by = 0.004)) +
   labs(title = " Age correlation with cg16867657 from 60 blood samples in humans", x = "correlation", y = "Count") +
   theme_minimal()
 
 ggplot(homo_com_pval60, aes(x = pval )) + 
   geom_histogram(
     colour = 1, fill = "white") +
   labs(title = "pvalues of age correlation with cg16867657 using 60 blood samples from humans", x = "correlation", y = "Count") +
   theme_minimal()
 
 
 
 
 ############################ loading data with 90 samples ################################### 
 blood_samples_90_correlation <- readRDS("~/Phd_data/outputs/blood_samples_90_correlation.rds")
 
 ## adding the rownames of the blood_correlation of Homo sapiens
 rownames(blood_samples_90_correlation[["Homo sapiens"]][["correlations"]]) <- rownames(blood_merged_betas)
 
 Homo_sapiens_corr90 <- (blood_samples_90_correlation[["Homo sapiens"]][["correlations"]])
 View(Homo_sapiens_corr90)
 
 ## Extarcting the pvalue
 rownames(blood_samples_90_correlation[["Homo sapiens"]][["p_values"]]) <- rownames(blood_merged_betas)
 Homo_sapiens_p_values90 <- blood_samples_90_correlation[["Homo sapiens"]][["p_values"]]
 View(Homo_sapiens_p_values90)
 
 ## converting the cpgsites data structure to dataframe
 homo_com_site90 <-  data.frame (Homo_sapiens_corr90["cg16867657", ])
 colnames(homo_com_site90) <- "corr"
 homo_com_site90$corr <- as.numeric(homo_com_site90$corr)
 View(homo_com_site90)
 
 ## converting the pvalues data structure to dataframe
 homo_com_pval90 <-  data.frame (Homo_sapiens_p_values90["cg16867657", ])
 colnames(homo_com_pval90) <- "pval"
 homo_com_pval90$pval <- as.numeric(homo_com_pval90$pval)
 View(homo_com_pval90)
 homo_com_site90 <- cbind(homo_com_site90, homo_com_pval90)
 
 
 
 ### histogram
 ggplot(homo_com_site90, aes(x = corr )) + 
   geom_histogram(
     colour = 1, fill = "white") + scale_x_continuous(breaks = seq(0.90, 1.0, by = 0.004)) +
   labs(title = " Age correlation with cg16867657 from 90 blood samples in humans", x = "correlation", y = "Count") +
   theme_minimal()
 
 ggplot(homo_com_pval90, aes(x = pval )) + 
   geom_histogram(
     colour = 1, fill = "white") +
   labs(title = "pvalues of age correlation with cg16867657 using 90 blood samples from humans", x = "correlation", y = "Count") +
   theme_minimal()
 
 
 
 
 
 
 ############################         OVis aries            ##################################
 ############################ loading data with 36 samples ################################### 
 blood_samples_36_correlation <- readRDS("~/Phd_data/outputs/blood_samples_36_correlation.rds")
 
 ## adding the rownames of the blood_correlation of Homo sapiens
 rownames(blood_samples_36_correlation[["Ovis aries"]][["correlations"]]) <- rownames(blood_merged_betas)
 
 ovis_corr36 <- (blood_samples_36_correlation[["Ovis aries"]][["correlations"]])
 View(ovis_corr36)
 
 ## Extarcting the pvalue
 rownames(blood_samples_36_correlation[["Ovis aries"]][["p_values"]]) <- rownames(blood_merged_betas)
 ovis_p_values36 <- blood_samples_36_correlation[["Ovis aries"]][["p_values"]]
 View(ovis_p_values36)
 
 ## converting the cpgsites data structure to dataframe
ovis_com_site36 <-  data.frame (ovis_corr36["cg16867657", ])
 colnames(ovis_com_site36) <- "corr"
 ovis_com_site36$corr <- as.numeric(ovis_com_site36$corr)
 View(ovis_com_site36)
 
 ## converting the pvalues data structure to dataframe
 ovis_com_pval36 <-  data.frame (ovis_p_values36["cg16867657", ])
 colnames(ovis_com_pval36) <- "pval"
 ovis_com_pval36$pval <- as.numeric(ovis_com_pval36$pval)
 View(ovis_com_pval36)
 ovis_com_site36 <- cbind(ovis_com_site36, ovis_com_pval36)
 
 
 
 ### histogram
 ggplot(ovis_com_site36, aes(x = corr )) + 
   geom_histogram(
     colour = 1, fill = "white") + scale_x_continuous(breaks = seq(0.05, 1.0, by = 0.05)) +
   labs(title = " Age correlation with cg16867657 from 36 blood samples in Ovis aries", x = "correlation", y = "Count") +
   theme_minimal()
 
 ggplot(ovis_com_pval36, aes(x = pval )) + scale_x_continuous(breaks = seq(0.01, 1.0, by = 0.01))
   geom_histogram(
     colour = 1, fill = "white") +
   labs(title = "pvalues of age correlation with cg16867657 using 36 blood samples from Ovis aries", x = "correlation", y = "Count") +
   theme_minimal()

   
   ################################sample 60 ################################################# 
   blood_samples_60_correlation <- readRDS("~/Phd_data/outputs/blood_samples_60_correlation.rds")
   
   ## adding the rownames of the blood_correlation of Homo sapiens
   rownames(blood_samples_60_correlation[["Ovis aries"]][["correlations"]]) <- rownames(blood_merged_betas)
   
   ovis_corr60 <- (blood_samples_60_correlation[["Ovis aries"]][["correlations"]])
   View(ovis_corr60)
   
   ## Extarcting the pvalue
   rownames(blood_samples_60_correlation[["Ovis aries"]][["p_values"]]) <- rownames(blood_merged_betas)
   ovis_p_values60 <- blood_samples_60_correlation[["Ovis aries"]][["p_values"]]
   View(ovis_p_values60)
   
   ## converting the cpgsites data structure to dataframe
   ovis_com_site60 <-  data.frame (ovis_corr60["cg16867657", ])
   colnames(ovis_com_site60) <- "corr"
   ovis_com_site60$corr <- as.numeric(ovis_com_site60$corr)
   View(ovis_com_site60)
   
   ## converting the pvalues data structure to dataframe
   ovis_com_pval60 <-  data.frame (ovis_p_values60["cg16867657", ])
   colnames(ovis_com_pval60) <- "pval"
   ovis_com_pval60$pval <- as.numeric(ovis_com_pval60$pval)
   View(ovis_com_pval60)
   ovis_com_site60 <- cbind(ovis_com_site60, ovis_com_pval60)
   
   
   
   ### histogram
   ggplot(ovis_com_site60, aes(x = corr )) + 
     geom_histogram(
       colour = 1, fill = "white") + scale_x_continuous(breaks = seq(0.05, 1.0, by = 0.05)) +
     labs(title = " Age correlation with cg16867657 from 60 blood samples in Ovis aries", x = "correlation", y = "Count") +
     theme_minimal()
   
   ggplot(ovis_com_pval60, aes(x = pval )) +
   geom_histogram(
     colour = 1, fill = "white") +  scale_x_continuous(breaks = seq(0.05, 1.0, by = 0.05)) +
     labs(title = "pvalues of age correlation with cg16867657 using 60 blood samples from Ovis aries", x = "correlation", y = "Count") +
     theme_minimal()
   
   
  
########################### Plot of correlation distribution ###################################
###########################                                 ####################################   
   
## This script code is used to plot the distribution of the estimated correlation with pvalues using samples greater than 10 from blood 
## And the correlation with pvalues based on sample size n = 12, 36, 60, 90
## This was used to determine which of the sample sizes will give the most robust correlation values
   
   
# Load necessary libraries
library(tidyverse)
library(ggdist)  # For stat_lineribbon()

   
# Load the estimated correlation file (CpG sites as rows, single correlation column)
estimated_data <- readRDS("~/Phd_data/estimated_correlation.RDS")
Homo_sapiens_corr <- readRDS("~/Phd_data/Homo_sapiens_corr_n12.RDS")

# checking the data
head(estimated_data)
estimated_data <- as.data.frame(estimated_data)

#changing the column name from Beta_ID to cpg_sites
estimated_data <- estimated_data %>% rename(cpg_site = cpg_sites )
head(estimated_data)

## Checking the data
Homo_sapiens_corr <- as.data.frame(Homo_sapiens_corr)
Homo_sapiens_corr <- rownames_to_column(Homo_sapiens_corr, var = "cpg_site")
rep_corr_mean <- as.data.frame(apply(Homo_sapiens_corr, 1, median, na.rm = TRUE))
head(Homo_sapiens_corr)
head(rep_corr_mean)

colnames(rep_corr_mean) <- "median_corr"
Homo_sapiens_corr1 <- cbind(estimated_data$cpg_site, rep_corr_mean$median_corr)
Homo_sapiens_corr1 <- as.data.frame(Homo_sapiens_corr1)
## adding colnames
colnames(Homo_sapiens_corr1) <- c("cpg_site", "median_corr")
str(Homo_sapiens_corr1)


# Convert the 50 replicate columns into long format
long_rep_corr <- Homo_sapiens_corr %>%
  pivot_longer(cols = starts_with("Rep_"), names_to = "replicate", values_to = "corr_value")

head(long_rep_corr)



# Merge with estimated correlation data
plot_data <- long_rep_corr %>%
  left_join(estimated_data, by = "cpg_site")
head(plot_data)

## Merge with Homo_sapiens_corr1
plot_data1 <- merge(plot_data, Homo_sapiens_corr1, by = "cpg_site")
head(plot_data1)

## small datasets
spl_dat2 <- plot_data1[1:10000,]
head(spl_dat2)

library(ggdist)
library(tidybayes)

ggplot(spl_dat2, aes(x = correlation, y = median_corr)) +
  stat_lineribbon(aes(y = corr_value), .width =  c(0.5, 0.95)) +  # 50%, 80%, 95% credible intervals
  scale_fill_brewer(palette = "Blues") +
  labs(
    title = "Spread of Replicate Correlations vs Estimated Correlation",
    x = "Estimated Correlation",
    y = "Replicate Correlation"
  ) +
  theme_minimal()


####################### for n = 12
# Load the estimated correlation file (CpG sites as rows, single correlation column)
estimated_data <- readRDS("~/Phd_data/estimated_correlation.RDS")
Homo_sapiens_corr <- readRDS("~/Phd_data/Homo_sapiens_corr_n12.RDS")

# checking the data
head(estimated_data)
estimated_data <- as.data.frame(estimated_data)

#changing the column name from Beta_ID to cpg_sites
estimated_data <- estimated_data %>% rename(cpg_site = Beta_ID )

## Checking the data
Homo_sapiens_corr <- as.data.frame(Homo_sapiens_corr)
Homo_sapiens_corr <- rownames_to_column(Homo_sapiens_corr, var = "cpg_site")
head(Homo_sapiens_corr)

# Convert the 50 replicate columns into lo5ng format
long_rep_corr <- Homo_sapiens_corr %>%
  pivot_longer(cols = starts_with("Rep_"), names_to = "replicate", values_to = "corr_value")

head(long_rep_corr)

# Merge with estimated correlation data
plot_data <- long_rep_corr %>%
  left_join(estimated_data, by = "cpg_site")
head(plot_data_36)

spl_dat1 <- plot_data[1:2000,]
head(spl_dat1)


library(ggdist)
library(tidybayes)

ggplot(plot_data, aes(x = correlation, y = corr_value)) +
  stat_lineribbon(aes(y = corr_value), .width =  c(0.5,0.80,0.95)) +  # 50%, 80%, 95% credible intervals
  scale_fill_brewer(palette = "Blues") +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.2)) + 
  labs(
    title = "Spread of Replicate Correlations vs Estimated Correlation N = 12",
    x = "Estimated Correlation",
    y = "Replicate Correlation"
  ) +
  theme_minimal()


ggplot(spl_dat1, aes(x = correlation, y = corr_value)) +
  stat_lineribbon(aes(y = corr_value), .width =  c(0.5,0.80,0.95)) +  # 50%, 80%, 95% credible intervals
  scale_fill_brewer(palette = "Blues") +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.2)) + 
  labs(
    title = "Spread of Replicate Correlations vs Estimated Correlation N = 12",
    x = "Estimated Correlation",
    y = "Replicate Correlation"
  ) +
  theme_minimal()
   
 
####################### for n =36
# Load the estimated correlation file (CpG sites as rows, single correlation column)
estimated_data <- readRDS("~/Phd_data/estimated_correlation.RDS")
Homo_sapiens_corr <- readRDS("~/Phd_data/Homo_sapiens_corr_n12.RDS")

# checking the data
head(estimated_data)
estimated_data <- as.data.frame(estimated_data)

#changing the column name from Beta_ID to cpg_sites
estimated_data <- estimated_data %>% rename(cpg_site = Beta_ID )

## Checking the data
Homo_sapiens_corr36 <- as.data.frame(Homo_sapiens_corr36)
Homo_sapiens_corr36 <- rownames_to_column(Homo_sapiens_corr36, var = "cpg_site")

# Convert the 50 replicate columns into long format
long_rep_corr36 <- Homo_sapiens_corr36 %>%
  pivot_longer(cols = starts_with("Rep_"), names_to = "replicate", values_to = "corr_value")

View(long_rep_corr36)

# Merge with estimated correlation data
plot_data_36 <- long_rep_corr36 %>%
  left_join(estimated_data, by = "cpg_site")
head(plot_data_36)

spl_dat1 <- plot_data[1:50000,]
head(spl_dat1)


library(ggdist)
library(tidybayes)

ggplot(plot_data_36, aes(x = correlation, y = corr_value)) +
  stat_lineribbon(aes(y = corr_value), .width =  c(0.5,0.80,0.95)) +  # 50%, 80%, 95% credible intervals
  scale_fill_brewer(palette = "Blues") +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.2)) + 
  labs(
    title = "Spread of Replicate Correlations vs Estimated Correlation N = 36",
    x = "Estimated Correlation",
    y = "Replicate Correlation"
  ) +
  theme_minimal()



####################### for n = 90
# Load the estimated correlation file (CpG sites as rows, single correlation column)
estimated_data <- readRDS("~/Phd_data/estimated_correlation.RDS")
Homo_sapiens_corr <- readRDS("~/Phd_data/Homo_sapiens_corr_n12.RDS")

# checking the data
head(estimated_data)
estimated_data <- as.data.frame(estimated_data)

#changing the column name from Beta_ID to cpg_sites
estimated_data <- estimated_data %>% rename(cpg_site = Beta_ID )

## Checking the data
Homo_sapiens_corr90 <- as.data.frame(Homo_sapiens_corr90)
Homo_sapiens_corr90 <- rownames_to_column(Homo_sapiens_corr90, var = "cpg_site")
head(Homo_sapiens_corr90)

# Convert the 50 replicate columns into long format
long_rep_corr90 <- Homo_sapiens_corr90 %>%
  pivot_longer(cols = starts_with("Rep_"), names_to = "replicate", values_to = "corr_value")

View(long_rep_corr90)
head(long_rep_corr90)
# Merge with estimated correlation data
plot_data_90 <- long_rep_corr90 %>%
  left_join(estimated_data, by = "cpg_site")
head(plot_data_90)

spl_dat3 <- plot_data_90[1:1000,]
head(spl_dat3)


library(ggdist)


ggplot(plot_data_90, aes(x = correlation)) +
  stat_lineribbon(aes(y = corr_value), .width = c(0.5,0.80,0.95)) +  # 50%, 80%, 95% credible intervals
  scale_fill_brewer(palette = "Blues") + 
  labs(
    title = "Spread of Replicate Correlations vs Estimated Correlation N = 90",
    x = "Estimated Correlation",
    y = "Replicate Correlation"
  ) + theme_minimal()


 


ggplot(spl_dat3, aes(x = correlation)) +
   stat_lineribbon(aes(y = corr_value)) +  # 80%, 95% credible intervals
  scale_fill_brewer(palette = "Blues") + 
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.2)) + 
  labs(
    title = "Spread of Replicate Correlations vs Estimated Correlation N = 90",
    x = "Estimated Correlation",
    y = "Replicate Correlation"
  ) +
  theme_minimal()

plot(plot_data$correlation, plot_data_90$corr_value)
 ####################### for n = 60
# Load the estimated correlation file (CpG sites as rows, single correlation column)
estimated_data <- readRDS("~/Phd_data/estimated_correlation.RDS")
Homo_sapiens_corr <- readRDS("~/Phd_data/Homo_sapiens_corr_n12.RDS")

# checking the data
head(estimated_data)
estimated_data <- as.data.frame(estimated_data)

#changing the column name from Beta_ID to cpg_sites
estimated_data <- estimated_data %>% rename(cpg_site = Beta_ID )

## Checking the data
Homo_sapiens_corr60 <- as.data.frame(Homo_sapiens_corr60)
Homo_sapiens_corr60 <- rownames_to_column(Homo_sapiens_corr60, var = "cpg_site")

# Convert the 50 replicate columns into long format
long_rep_corr60 <- Homo_sapiens_corr60 %>%
  pivot_longer(cols = starts_with("Rep_"), names_to = "replicate", values_to = "corr_value")

View(long_rep_corr60)

# Merge with estimated correlation data
plot_data <- long_rep_corr60 %>%
  left_join(estimated_data, by = "cpg_site")
head(plot_data)

spl_dat2 <- plot_data[1:50000,]
head(spl_dat2)


library(ggdist)
library(tidybayes)

ggplot(plot_data, aes(x = correlation, y = corelation)) +
  stat_lineribbon(aes(y = corr_value), .width = c(0.5,0.80,0.95)) +  # 50%, 80%, 95% credible intervals
  scale_fill_brewer(palette = "Blues") + 
  labs(
    title = "Spread of Replicate Correlations vs Estimated Correlation N = 60",
    x = "Estimated Correlation",
    y = "Replicate Correlation"
  ) +
  theme_minimal()




# ON SAMPLE DATA
set.seed(12345)
tibble(
  x = rep(1:20, 1000),
  y = rnorm(20000, mean = 0.5)
) %>%
  ggplot(aes(x = x, y = y)) +
  stat_lineribbon() +
  scale_fill_brewer()




































































########### comparing with n = 36 ##############  
  # Load the estimated correlation file (CpG sites as rows, single correlation column)
  estimated_data <- final_results %>% filter(species == "Homo sapiens") %>% filter (tissue=="Blood")
  rownames(estimated_data) <- rownames(Homo_sapiens_corr)
  View(estimated_data)
  
  # Calculate the 95th percentile for each row (CpG site)
  replicate_36 <- as.data.frame(apply(Homo_sapiens_corr36, 1, median, na.rm = TRUE))
  colnames(replicate_36) <- "corr_median"
  View(replicate_36)
  
  # Convert to a tibble for easier plotting
  n36_plot <- cbind(replicate_36$corr_median, estimated_data$correlation) 
  colnames(n36_plot) <- c("n36_corr", "est_corr")
  View(n36_plot)
  # Plot using ggplot and stat_lineribbon
  ggplot(n36_plot, aes(x = est_corr, y = n36_corr)) + geom_point(color = "navy")
  
  ## stat_lineribbon() +
  scale_fill_brewer() +
    labs(title = "Estimated Correlation vs. 95th Percentile of Replicate Correlations",
         x = "Estimated Correlation",
         y = "95th Percentile of Replicate Correlations") +
    theme_minimal()
  
  
  ########### comparing with n = 60 ##############  
  # Load the estimated correlation file (CpG sites as rows, single correlation column)
  estimated_data <- final_results %>% filter(species == "Homo sapiens") %>% filter (tissue=="Blood")
  rownames(estimated_data) <- rownames(Homo_sapiens_corr)
  View(estimated_data)
  
  # Calculate the 95th percentile for each row (CpG site)
  replicate_60 <- as.data.frame(apply(Homo_sapiens_corr60, 1, median, na.rm = TRUE))
  colnames(replicate_60) <- "corr_median"
  View(replicate_60)
  
  # Convert to a tibble for easier plotting
  n60_plot <- cbind(replicate_60$corr_median, estimated_data$correlation) 
  colnames(n60_plot) <- c("n60_corr", "est_corr")
  View(n60_plot)
  # Plot using ggplot and stat_lineribbon
  ggplot(n60_plot, aes(x = est_corr, y = n60_corr)) + geom_point(color = "navy")
  
  ## stat_lineribbon() +
  scale_fill_brewer() +
    labs(title = "Estimated Correlation vs. 95th Percentile of Replicate Correlations",
         x = "Estimated Correlation",
         y = "95th Percentile of Replicate Correlations") +
    theme_minimal()
  
  ########### comparing with n = 90 ##############  
  # Load the estimated correlation file (CpG sites as rows, single correlation column)
  estimated_data <- final_results %>% filter(species == "Homo sapiens") %>% filter (tissue=="Blood")
  rownames(estimated_data) <- rownames(Homo_sapiens_corr)
  View(estimated_data)
  
  # Calculate the 95th percentile for each row (CpG site)
  replicate_90 <- as.data.frame(apply(Homo_sapiens_corr90, 1, median, na.rm = TRUE))
  colnames(replicate_90) <- "corr_median"
  View(replicate_90)
  
  # Convert to a tibble for easier plotting
  n90_plot <- cbind(replicate_90$corr_median, estimated_data$correlation) 
  colnames(n90_plot) <- c("n90_corr", "est_corr")
  View(n90_plot)
  # Plot using ggplot and stat_lineribbon
  ggplot(n90_plot, aes(x = est_corr, y = n90_corr)) + 
  stat_lineribbon() +
  scale_fill_brewer() +
    labs(title = "Estimated Correlation vs. 95th Percentile of Replicate Correlations",
         x = "Estimated Correlation",
         y = "95th Percentile of Replicate Correlations") +
    theme_minimal()
  
