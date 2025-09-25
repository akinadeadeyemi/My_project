######### Author : Adeyemi Akinade
######### Place : Clemson University, SC
######### lab : Gopalan's lab
######### Date : June 16, 2025
######### Project : Comparative Epigenetics of Ageing
######### This script was written to determine the linear regression between age and beta values of DNA methylation sites in tissues


######################## VISUALIZING PHYLOGENETIC REGRESSION MODEL IN BLOOD SAMPLES ONLY ##########################################

### Selecting and Analysing the significant
library(dplyr)

#load blood_phylo_result
blood_phylo_result <- readRDS("/Users/aakinad/Downloads/PHD FOLDER/blood_phylo_result.rds")
blood_phylo_rel_age_result  <- readRDS("/Users/aakinad/Downloads/PHD FOLDER/blood_phylo_age_result.rds")
blood_phylo_sig <- blood_phylo_result %>% filter(p_val < 0.05)

## viewing significant CpG sites.
View(blood_phylo_sig)

## selecting the CpG sites with < 0.05
blood_phylo_lambda_sig <- blood_phylo_sig %>% filter(lambda < 1.0 & lambda > 0.05)

blood_phylo_lambda_sig <- blood_phylo_lambda_sig %>% arrange(desc(lambda))
View(blood_phylo_lambda_sig)


### Testing the relationship between raw DNA methylation levels with age_in_years across species.
## topmost CpG site with high lambda and sigma2 cg27413543

## loading the raw beta values
blood_metadata_wt_life_span_betas <- readRDS("/Users/aakinad/Downloads/PHD FOLDER/blood_metadata_wt_life_span_betas.rds")
blood_metadata_wt_life_span <- readRDS("/Users/aakinad/Downloads/PHD FOLDER/blood_metadata_wt_life_span.rds")

## selecting the cpg sites
cg_27413543_betas <- t(blood_metadata_wt_life_span_betas["cg16867657", ])

## rownames to column names
library(tibble)
cg_27413543_betas <- as.data.frame(cg_27413543_betas)
cg_27413543_betas <- rownames_to_column(cg_27413543_betas, var = "geo_accession")

blood_metadata_age_years <- blood_metadata_wt_life_span %>% select(geo_accession, age_years, species, max_lifespan, order)

### adding longevity criteria
blood_metadata_age_years <- blood_metadata_wt_life_span %>%
  select(geo_accession, age_years, species, max_lifespan, order) %>%
  mutate(
    longevity = case_when(
      max_lifespan > 15 & order %in% c("Rodentia", "Chiroptera") ~ "exceptionally long-lived",
      max_lifespan > 15 ~ "long-lived",
      TRUE ~ "short-lived"
    )
  )

### merging the raw betas with their age data.
cg_27413543_betas_with_age <- merge(cg_27413543_betas, blood_metadata_age_years, by = "geo_accession")

cg_27413543_betas_with_age <- cg_27413543_betas_with_age %>%
  mutate(longevity_class = ifelse(max_lifespan > 15, "Long-lived", "Short-lived"))


library(ggplot2)
library(ggrepel)

ggplot(cg_27413543_betas_with_age, aes(x = age_years, y = cg22704389, color = species)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~species, scales = "free_x") +  # <- one panel per order
  theme_minimal() +
  theme(
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "italic"),
    strip.text = element_text(size = 13, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  labs(
    title = "DNA Methylation at cg22704389 vs Age by Mammalian Order",
    x = "Age (years)",
    y = "DNAm methylation level",
    color = "Species"
  ) +
  theme(legend.position = "top")  # Hide legend to reduce clutter



library(dplyr)
library(ggplot2)

# Choose 4 species
four_species <- c("Homo sapiens", "Mus musculus", "Rattus norvegicus", "Heterocephalus glaber")

cg_subset <- cg_27413543_betas_with_age %>%
  filter(species %in% four_species)

ggplot(cg_subset, aes(x = age_years, y = cg16867657, color = species)) +
  geom_point(size = 1.8, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~species, scales = "free_x") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "italic"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  ) +
  labs(
    title = "DNA Methylation at cg00270194 vs Age (4 Species)",
    x = "Age (years)",
    y = "DNAm methylation level",
    color = "Species"
  ) +
  theme(legend.position = "top")


library(ggplot2)
library(dplyr)
library(tibble)

plot_cpg_by_species_facet <- function(cpg_id, betas_df, metadata_df) {
  # Get betas for the selected CpG and make a data.frame
  betas <- t(betas_df[cpg_id, , drop = FALSE])
  betas <- as.data.frame(betas)
  colnames(betas) <- c(cpg_id)
  betas <- rownames_to_column(betas, var = "geo_accession")
  
  # Merge with metadata
  merged <- merge(betas, metadata_df, by = "geo_accession")
  
  # Create longevity class
  # merged <- merged %>%
   # mutate(longevity_class = ifelse(max_lifespan > 15, "Long-lived", "Short-lived"))
  
  # Plot
  p <- ggplot(merged, aes(x = age_years, y = !!sym(cpg_id), color = species)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~species, scales = "free_x") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(color = "black", size = 10),
      axis.text.y = element_text(color = "black", size = 10),
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10, face = "italic"),
      strip.text = element_text(size = 10, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
    ) +
    labs(
      title = paste("DNAm at", cpg_id, "vs Age by Species"),
      x = "Age (years)",
      y = "DNAm methylation level",
      color = "Species"
    ) +
    theme(legend.position = "none")
  
  return(p)
}

# Read the file (assumes CpGs are in rownames)
top20_cpgs <- rownames(blood_phylo_sig)[1:20]

# Create and store plots
library(patchwork)
plot_list <- lapply(top20_cpgs, plot_cpg_by_species_facet,
                    betas_df = blood_metadata_wt_life_span_betas,
                    metadata_df = blood_metadata_age_years)

# Combine and display or save (e.g., 4 at a time)
wrap_plots(plot_list[5], ncol = 1)





##### DNA METHYLATION RATE(FROM LINEAR MODEL) WITH LIFESPAN (PER SPECIES) ####################
#####                                                                     ####################
#####                                                                     ####################
blood_regression_slope_wide <- readRDS("~/Downloads/blood_regression_slope_wide.rds")

## selecting the cpg sites
blood_regression_cg27413543_betas <- t(blood_regression_slope_wide["cg27413543", ])
blood_regression_cg27413543_betas <- as.data.frame(blood_regression_cg27413543_betas)
blood_regression_cg27413543_betas <- rownames_to_column(blood_regression_cg27413543_betas, var = "species")

### merging the raw betas with their age data.
blood_metadata_max_lifespan <- unique(blood_metadata_age_years[, c("species", "max_lifespan")])
## To ensure the species name match in the two data sets
blood_regression_cg27413543_betas$species <- gsub("_", " ", blood_regression_cg27413543_betas$species)


blood_regression_cg27413543_data <- merge(blood_regression_cg27413543_betas, blood_metadata_max_lifespan, by = "species")



library(ggplot2)
library(ggpubr)
library(ggrepel)

ggplot(blood_regression_cg27413543_data, aes(x = (max_lifespan), y = (cg27413543), color = species)) +
  geom_point(size = 2, alpha = 0.7, color ="black") +
  geom_smooth(method = "lm", se = T, color = "black") +
  geom_text_repel(aes(label = species), size = 2.5, max.overlaps = 100) +
  stat_cor(method = "pearson", size = 2, label.x.npc = "left", label.y.npc = "top") +
  scale_x_log10() +
  labs(
    title = "Correlation Between Log(Max Lifespan) and DNAm at cg27413543",
    x = "Log(Maximum Lifespan, years)",
    y = "DNAm Methylation Rate"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.y = element_text(face = "italic"),
    legend.position = "none"
  ) + 
  theme(axis.text.x = element_text(color = "black",
                                   size = 14),
        axis.text.y = element_text(color = "black",
                                   size = 14,
                                   hjust = 1),
        axis.title.x = element_text(size = 15,
                                    color = "black",
                                    face = "bold"),
        axis.title.y = element_text(size = 15,
                                    color = "black",
                                    face = "italic"))
## using a plotting function
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(tibble)

plot_cpg_vs_lifespan <- function(cpg_id, slope_wide, metadata) {
  # Extract CpG methylation rates
  betas <- t(slope_wide[cpg_id, , drop = FALSE])
  betas_df <- as.data.frame(betas)
  betas_df <- rownames_to_column(betas_df, var = "species")
  
  # Clean species names
  betas_df$species <- gsub("_", " ", betas_df$species)
  
  # Merge with max lifespan
  merged_df <- merge(betas_df, metadata, by = "species")
  colnames(merged_df)[colnames(merged_df) == cpg_id] <- "methylation_rate"
  
  # Plot
  p <- ggplot(merged_df, aes(x = max_lifespan, y = methylation_rate)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    geom_text_repel(aes(label = species), size = 2.5, max.overlaps = 100) +
    stat_cor(method = "pearson", size = 2, label.x.npc = "left", label.y.npc = "top") +
    labs(
      title = paste("Lifespan vs DNAm at", cpg_id),
      x = "(Maximum Lifespan, years)",
      y = "DNAm Methylation Rate"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title.y = element_text(face = "italic"),
      legend.position = "none"
    )
  
  return(p)
}

# Select the first 20 CpG IDs
first_20_cpgs <- rownames(blood_phylo_lambda_sig)[1730:1800]

# Create plots using lapply
plot_list <- lapply(first_20_cpgs, function(cpg) {
  plot_cpg_vs_lifespan(cpg, blood_regression_slope_wide, blood_metadata_max_lifespan)
})

# Save all plots
for (i in seq_along(plot_list)) {
  ggsave(filename = paste0("CpG_plot_", first_20_cpgs[i], ".png"),
         plot = plot_list[[i]], width = 6, height = 5, dpi = 300)
}
#install.packages("patchwork")
library(patchwork)
wrap_plots(plot_list[1:4], ncol = 2, nrow = 2)












### histogram of blood phylosig lamda and sigma2
ggplot(blood_phylo_result, aes(x = sig2)) +
  geom_density(binwidth = 3, fill = "blue", color = "black", alpha = 0.7) +
  labs(
    title = expression("Distribution of Evolutionary rate " * sigma^2, "in blood tissue"),
    x = expression("Evolutionary Rate " * sigma^2),
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14)
  )





###### QQPLOT in ggplot2
library(ggplot2)

# Prepare data for relative age
pval <- blood_regression_result_wt_padjust$Bonf_adj
observed <- -log10(sort(pval))
expected <- -log10(ppoints(length(pval)))
qq_data <- data.frame(Expected = expected, Observed = observed)

# Prepare data for age_in_years
pval <- blood_regression_rel_age_result_wt_padjust$Bonf_adj
observed <- -log10(sort(pval))
expected <- -log10(ppoints(length(pval)))
qq_data <- data.frame(Expected = expected, Observed = observed)


# Plot
ggplot(qq_data, aes(x = Expected, y = Observed)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  labs(
    title = expression("QQ Plot of DNAm rate (relative age) adjusted P-values"),
    x = expression(Expected ~ -log[10](italic(p))),
    y = expression(Observed ~ -log[10](italic(p)))
  ) +
  theme_minimal(base_size = 14)



#### Calculating for p-adjust for linear models
blood_regression_result <- readRDS("~/Downloads/blood_regression_result.rds")
blood_regression_rel_age_result <- readRDS("~/Downloads/blood_regression_test_wt_lifespan_result.rds")
View(blood_regression_result)

### Blood regression with age_in_years with p_adjust
library(dplyr)
blood_regression_result_wt_padjust <- blood_regression_result %>%
  group_by(species) %>%
  mutate(Bonf_adj = p.adjust(p_value, method = "bonferroni")) 

### Blood regression with relative age with p_adjust
blood_regression_rel_age_result_wt_padjust <- blood_regression_rel_age_result %>%
  group_by(species) %>%
  mutate(Bonf_adj = p.adjust(p_value, method = "bonferroni")) 

### Selecting the corrected pval significant sites
blood_regression_sig_cpgsites <- blood_regression_result_wt_padjust %>%
  filter(Bonf_adj < 0.05)

blood_regression_rel_age_sig_cpgsites <- blood_regression_rel_age_result_wt_padjust %>%
  filter(Bonf_adj < 0.05)

length(unique(blood_regression_rel_age_sig_cpgsites$cpg_site))

### Do p.adjust and turn rows to columns
library(tibble)
blood_phylo_result_wt_padjust <- blood_phylo_result %>%
  mutate(Bonf_adj = p.adjust(p_val, method = "bonferroni")) %>% 
  rownames_to_column(var = "cpg_site")

blood_phylo_rel_age_result_wt_padjust <- blood_phylo_rel_age_result %>%
  mutate(Bonf_adj = p.adjust(p_val, method = "bonferroni")) %>% 
  rownames_to_column(var = "cpg_site")

#### Histogram of p.adjust coloured by species
library(ggplot2)
custom_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF00FF", "#838B8B", "#0000FF", "#050505","#00CD00","#FF1493", "#FFC1C1","#B03060", "#2F4F4F","#551A8B",
  "#FF7F00", "#A65628", "#00F5FF", "#F781BF", "#CD8500", "#66C2A5", "#7CFC00", "#0000CD","#8B2323","#6B8E23", "#87CEFF", "#A4D3EE", "#CD2990", 
  "#FC8D62", "#458B74", "#CDBA96", "#8B3626", "#006400", "#FFD700", "#8B6508","#00008B", "#9370DB","#96CDCD", "#8B4500", "#2F4F4F",  "#87CEFA"
)

ggplot(blood_regression_sig_cpgsites, aes(x = Bonf_adj, fill = species)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 50) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Distribution of Bonferroni-adjusted lm p-values by Species",
       x = "Bonferroni-adjusted p-value",
       y = "Count") +
  theme_minimal()

  
 ## creating a table for the blood regression significant
blood_regression_sig_cpgsites_table <- table(blood_regression_sig_cpgsites$species)

blood_regression_sig_cpgsites_species <- split(blood_regression_sig_cpgsites,blood_regression_sig_cpgsites$species)


###### Selecting for cpgsites that are significant in linear model but not phylosig analysis
### Step 1. 
# For age_in_years
blood_phylo_lambda_sig <- blood_phylo_result_wt_padjust %>%
  filter(lambda > 0 & Bonf_adj < 0.05)
# For relative age 
blood_phylo_rel_age_lambda_sig <- blood_phylo_rel_age_result_wt_padjust %>%
  filter(lambda > 0 & Bonf_adj < 0.05) 

library(tibble)  # for rownames_to_column
blood_phylo_lambda_sig <- blood_phylo_lambda_sig %>%
  rownames_to_column(var = "cpg_site")

### step 2  
library(dplyr)
blood_species_specific_sig <- blood_regression_sig_cpgsites %>%
  anti_join(blood_phylo_lambda_sig, by = "cpg_site") %>%
  group_by(species)

blood_species_specific_sig_rel_age <- blood_regression_rel_age_sig_cpgsites %>%
  anti_join(blood_phylo_rel_age_lambda_sig, by = "cpg_site") %>%
  group_by(species)

saveRDS(blood_species_specific_sig, 
               file = "/Users/aakinad/Downloads/PHD FOLDER/blood_age_only_associated_cpg_freq.RDS")
length(unique(blood_species_specific_sig$cpg_site))


### step 3: spliting the age_asociated_only cpg sites across species
# Create a named list: one vector of CpGs per species
blood_age_associated_only_cpg_list <- split(blood_species_specific_sig$cpg_site, blood_species_specific_sig$species)

install.packages("UpSetR")
library(UpSetR)
install.packages("ComplexUpset")
library(ComplexUpset)

# Convert list to a binary presence/absence matrix
upset_input <- fromList(blood_age_associated_only_cpg_list)
upset_input_tbl <- as_tibble(upset_input)


UpSetR::upset(
  upset_input,
  nsets = 32,             # Number of species to include
  nintersects = 120,       # Number of top combinations to show
  order.by = "freq",      # Sort by most frequent intersections
  sets.bar.color = "navyblue",
  mainbar.y.label = "Number of CpG Sites",
  sets.x.label = "CpGs per Species",
  text.scale = c(0.5, 0.5, 0.5, 0.5, 0.7, 0.5)
 )

#### going from the most shared 
UpSetR::upset(
  upset_input,
  nsets = 32,
  nintersects = 120,
  order.by = "freq",   # Orders by number of sets (species) in each intersection
  decreasing = TRUE,    # Still from low to high degree (not always honored)
  sets.bar.color = "navyblue",
  mainbar.y.label = "Number of CpG Sites",
  sets.x.label = "CpGs per Species",
  text.scale = c(0.8, 0.8, 0.8, 0.8, 1.0, 0.8)
)


################ SELECTING FOR SIGNIFICANT Age_associated CpG only SITES THAT ARE SHARED AMONG ALL SPECIES IN LINEAR MODEL  ##################
# Here, we are counting how many age-associated_only cpg sites that are shared among species. we are taking cpgsites that are
# age-associated but do not show phylogenetic signal in each species. Then we are extracting the cpgsites are shared across these species
## This helps us to know sites that show convergent evolution (i.e evolve independently/ no shared ancestry) due to functional conservation 
## without phylogenetic dependence or core aging mechanism as well as tissue specific response aging process due to similar environmental exposures
## or pressures hence making them appearing shared.

# For age_in_years
blood_species_specific_cpg_counts <- blood_species_specific_sig %>%
  group_by(cpg_site) %>%
  summarise(n_species = n_distinct(species)) #%>%
  #filter(n_species > 1)   # Shared CpGs: appear in 2+ species

# For relative age
blood_species_specific_cpg_counts_rel_age <- blood_species_specific_sig_rel_age %>%
  group_by(cpg_site) %>%
  summarise(n_species = n_distinct(species))


## save the count of age-associated_only_cpgsites
saveRDS(blood_species_specific_cpg_counts, 
        file = "/Users/aakinad/Downloads/PHD FOLDER/blood_species_specific_cpg_counts.RDS")


#### VISUALIZING THE PATTERNS OF AGE-ASSOICATED-ONLY AND PHYLOSIG SHARED CPGS ACROSS SPECIES
### joining back the age-associated-only-cpgsites to the original age-associated cpg information
### these are age-associated_only cpgs that are shared by at least two or more species
blood_species_specific_shared_cpgs <- blood_species_specific_sig %>%
semi_join(blood_species_specific_cpg_counts, by = "cpg_site") 

 View(blood_species_specific_cpg_counts)


plot_cpg_by_species_facet <- function(cpg_id, betas_df, metadata_df) {
  library(ggplot2)
  library(dplyr)
  library(tibble)
  
  # Get betas for the selected CpG and make a data.frame
  betas <- t(betas_df[cpg_id, , drop = FALSE])
  betas <- as.data.frame(betas)
  colnames(betas) <- c(cpg_id)
  betas <- rownames_to_column(betas, var = "geo_accession")
  
  # Merge with metadata
  merged <- merge(betas, metadata_df, by = "geo_accession")
  
  # Define color palette
  custom_colors <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF00FF", "#838B8B", "#0000FF", "#050505", "#00CD00", "#FF1493",
    "#FFC1C1", "#B03060", "#2F4F4F", "#551A8B", "#FF7F00", "#A65628", "#00F5FF", "#F781BF", "#CD8500", "#66C2A5",
    "#7CFC00", "#0000CD", "#8B2323", "#6B8E23", "#87CEFF", "#A4D3EE", "#CD2990", "#FC8D62", "#458B74", "#CDBA96",
    "#8B3626", "#006400", "#FFD700", "#8B6508", "#00008B", "#9370DB", "#96CDCD", "#8B4500", "#2F4F4F", "#87CEFA"
  )
  
  # Ensure species is a factor with consistent order
  merged$species <- factor(merged$species)
  
  # Plot
  p <- ggplot(merged, aes(x = age_years, y = !!sym(cpg_id), color = species)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~species, scales = "free_x") +
    scale_color_manual(values = custom_colors[1:length(unique(merged$species))]) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(color = "black", size = 10),
      axis.text.y = element_text(color = "black", size = 10),
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10, face = "italic"),
      strip.text = element_text(size = 10, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
    ) +
    labs(
      title = paste("DNAm vs Age at", cpg_id, "aacross species"),
      x = "Age (years)",
      y = "DNAm methylation level",
      color = "Species"
    ) +
    theme(legend.position = "none")
  
  return(p)
}



#### New function with bonferonni correction

plot_cpg_by_species_facet <- function(cpg_id, betas_df, metadata_df, blood_phylo_result_wt_padjust) {
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(rlang)
  
  # Prepare beta data
  betas <- t(betas_df[cpg_id, , drop = FALSE])
  betas <- as.data.frame(betas)
  colnames(betas) <- c(cpg_id)
  betas <- rownames_to_column(betas, var = "geo_accession")
  
  # Merge with metadata
  merged <- merge(betas, metadata_df, by = "geo_accession")
  
  # Ensure species is factor
  merged$species <- factor(merged$species)
  
  ## Adding the phylogenetic info
  blood_phylo_info <- blood_phylo_result_wt_padjust %>% filter(cpg_site == cpg_id)
  
  phylo_label <- paste0("Phylo signal:\nλ = ", signif(blood_phylo_info$lambda, 3),
                        ", sig2 = ", signif(blood_phylo_info$sig2, 3),
                        ", Bonf_adj = ", signif(blood_phylo_info$Bonf_adj, 3))
  
  # Fit linear models per species to get slope, R², and p-value
  pvals <- merged %>%
    group_by(species) %>%
    summarise(
      slope = tryCatch(
        coef(lm(!!sym(cpg_id) ~ age_years, data = pick(everything())))[2],
        error = function(e) NA_real_
      ),
      r2 = tryCatch(
        summary(lm(!!sym(cpg_id) ~ age_years, data = pick(everything())))$r.squared,
        error = function(e) NA_real_
      ),
      slope_p = tryCatch(
        summary(lm(!!sym(cpg_id) ~ age_years, data = pick(everything())))$coefficients[2, "Pr(>|t|)"],
        error = function(e) NA_real_
      )
    )
  
  # Adjust p-values (FDR)
  pvals <- pvals %>%
    mutate(adj_p = p.adjust(slope_p, method = "bonferroni"))
  
  # Prepare labels for annotation
  pvals <- pvals %>%
    mutate(
      label = paste0("slope = ", signif(slope, 3)," ", "R² = ", signif(r2, 3) ," ", "p_val.adj = ", signif(adj_p, 3)))
  
  
  # Define color palette
  custom_colors <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF00FF", "#838B8B", "#0000FF", "#050505", "#00CD00", "#FF1493",
    "#FFC1C1", "#B03060", "#2F4F4F", "#551A8B", "#FF7F00", "#A65628", "#00F5FF", "#F781BF", "#CD8500", "#66C2A5",
    "#7CFC00", "#0000CD", "#8B2323", "#6B8E23", "#87CEFF", "#A4D3EE", "#CD2990", "#FC8D62", "#458B74", "#CDBA96",
    "#8B3626", "#006400", "#FFD700", "#8B6508", "#00008B", "#9370DB", "#96CDCD", "#8B4500", "#2F4F4F", "#87CEFA"
  )
  
  # Plot
  p <- ggplot(merged, aes(x = age_years, y = !!sym(cpg_id), color = species)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~species, scales = "free_x") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
    scale_color_manual(values = custom_colors[1:length(unique(merged$species))]) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(color = "black", size = 10),
      axis.text.y = element_text(color = "black", size = 10),
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10, face = "italic"),
      strip.text = element_text(size = 10, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
    ) +
    labs(
      title = paste("DNAm vs Age patterns at", cpg_id, "across species"),
      caption = phylo_label,
      x = "Age (years)",
      y = "DNAm methylation level",
      color = "Species"
    ) +
    theme(legend.position = "none") +
    geom_text(
      data = pvals,
      aes(x = -Inf, y = Inf, label = label),
      hjust = -0.1,
      vjust = 1.2,
      inherit.aes = FALSE,
      size = 3,
      color = "black"
    )
  
  return(p)
}


##### using files of blood_phylo_lambda_sig: This code provide the analysis of selecting those sites that are both age-associated and 
#### phylogenetically dependent/shared ancestry across 34 mammalian species
# Read the file (assumes CpGs are in rownames)
blood_phylo_lambda_sig <- blood_phylo_lambda_sig %>% arrange(desc(lambda))

### For CpG sites that is significant in only one species
blood_species_specific_cpg_counts_n1 <- blood_species_specific_sig %>%
  group_by(cpg_site) %>%
  summarise(n_species = n_distinct(species)) %>%
  filter(n_species == 1) 

#top20_cpgs <- (blood_species_specific_cpg_counts_n1$cpg_site)[1:20]
top20_cpgs <- (blood_phylo_lambda_sig$cpg_site)[1:20]


blood_species_specific_cpg_counts <- blood_species_specific_cpg_counts %>% arrange(n_species)


### Longevity species - looking at shared age-associated- only cpg sites
longevity_species <- blood_metadata_age_years %>% filter(longevity %in% c("long-lived", "exceptionally long-lived")) 
targeted_long_lived_species <- unique(longevity_species$species)+

long_lived_shared_age_only_cpgs <- blood_species_specific_shared_cpgs %>% filter(species %in% targeted_long_lived_species)

interest_species<- c("Heterocephalus glaber", "Marmota flaviventris", "Homo sapiens", "Macaca mulatta", "Loxodonta africana") 
few_species <- c("Canis lupus familiaris", "Homo sapiens", "Mus musculus", "Rattus norvegicus")

#get the metadata of the few species
few_sp_metadata <- blood_metadata_wt_life_span_and_LQ %>% filter(species %in% few_species )

long_lived <- blood_species_specific_shared_cpgs %>% filter(species %in% interest_species)
long_lived_metadata_age_years <- blood_metadata_age_years %>% filter(species %in% interest_species )

long_lived_species_specific_cpg_counts <- long_lived %>%
  group_by(cpg_site) %>%
  summarise(n_species = n_distinct(species)) %>%
  filter(n_species > 1)  %>% arrange(desc(n_species)) # Shared CpGs: appear in 2+ species

### Short-lived species
short_lived <- blood_metadata_age_years %>% filter(longevity %in% "short-lived") 
short_lived_species <- unique(short_lived$species)
short_lived_shared <- blood_species_specific_shared_cpgs %>% filter(species %in% short_lived_species)
short_lived_metadata_age_years <- blood_metadata_age_years %>% filter(species %in% short_lived_species )

# shared cpg sites across short-lived species
short_lived_species_specific_cpg_counts <- short_lived_shared %>%
  group_by(cpg_site) %>%
  summarise(n_species = n_distinct(species)) %>%
  filter(n_species > 1)  %>% arrange(desc(n_species)) # Shared CpGs: appear in 2+ species




## selecting top 20 sites shared among 
top20_cpgs <- (long_lived_species_specific_cpg_counts$cpg_site)[1:50]
top30_cpgs <- (short_lived_species_specific_cpg_counts$cpg_site)[1:50]



rodent_metadata_age_years <- blood_metadata_age_years %>% filter(order== "Rodentia")

rodent_species_specific_cpg_counts <- rodent_shared_age_only_cpgs %>%
  group_by(cpg_site) %>%
  summarise(n_species = n_distinct(species)) %>%
  filter(n_species > 1)  %>% arrange((n_species)) # Shared CpGs: appear in 2+ species




### selecting the top 100 cpgs with highest sigma2 values and seeing their patterns in rodents
blood_phylo_result_wt_padjust <- blood_phylo_result_wt_padjust %>% arrange(desc(sig2))
blood_phylo_result_wt_padjust <- blood_phylo_result_wt_padjust %>% rownames_to_column(var = "cpg_site")

rodent_phylo_result_wt_padjust <- blood_phylo_result_wt_padjust %>% filter(cpg_site %in% 
                                                                             rodent_shared_age_only_cpgs$cpg_site

top100_cpgs <- (rodent_phylo_result_wt_padjust$cpg_site)[1:100]



# Create and store plots
library(patchwork)
plot_list <- lapply(top20_cpgs, plot_cpg_by_species_facet,
                    betas_df = blood_metadata_wt_life_span_betas,
                    metadata_df = few_sp_metadata,
                    blood_phylo_result_wt_padjust = blood_phylo_result_wt_padjust)

 # Combine and display or save (e.g., 4 at a time)
wrap_plots(plot_list[1], ncol = 1)
wrap_plots(plot_list[2], ncol = 1)
wrap_plots(plot_list[3], ncol = 1)


####### VENN DIGRAM OF AGE_ASSOCIATED AND PHYLOSIG CPGS IN BLOOD DNAm
### showing Venn diagram of shared cpgsites that are phylogenetic and at least one is age-associated.
install.packages("ggVennDiagram")
library(ggVennDiagram)

# preparing the sets as a named list
venn_data <- list(
  Age_Associated = blood_regression_sig_cpgsites$cpg_site,
  Phylogenetic = blood_phylo_lambda_sig$cpg_site
)

# Plotting the Venn diagram
ggVennDiagram(venn_data, label = "count") +
  scale_fill_gradient(low = "lightgreen", high = "darkgreen") +
  labs(title = "Venn Diagram of Age-Associated vs Phylogenetic CpGs") +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "bottom",  # Hide legend if not needed
    text = element_text(size = 8)
  )


library(VennDiagram)
library(grid)
library(futile.logger)

venn.plot <- draw.pairwise.venn(
  area1 = length(unique(univ_clock_cpgs$age_in_years)),  # Age-Associated
  area2 = length(unique(univ_clock_cpgs$rel_age)),         # Phylogenetic
  cross.area = length(intersect(
    unique(univ_clock_cpgs$age_in_years),
    unique(univ_clock_cpgs$rel_age)
  )),
  category = c("Age (in years)",
               "Relative Age"),
  fill = c("purple", "orange"),
  lty = "blank",
  cex = 2,
  cat.cex = 1.5,
  cat.pos = c(-30, 30),   # For horizontal layout
  cat.dist = 0.1,
  cat.just = list(c(1,0), c(0, 1)),
  euler.d = FALSE,
  scaled = FALSE
)




################### COMMON CPG SITES BETWEEN PHYLOGENETIC SIGNAL CPGS AND UNIVERSAL DNAm EPIGENETIC CLOCK-2 CPGs
# Get common CpG sites between phylogenetic signal cpgs and universal mammalian cpgs
phylo_mammalian_shared_cpgs <- intersect(blood_phylo_lambda_sig$cpg_site, cpg_sites_from_mammalian_universal_dnam$cpg_sites)
length(phylo_mammalian_shared_cpgs)

age_associated_only_mammalian_shared_cpgs <- intersect(blood_species_specific_shared_cpgs$cpg_site, cpg_sites_from_mammalian_universal_dnam$cpg_sites)
length(age_associated_only_mammalian_shared_cpgs)
View(age_associated_only_mammalian_shared_cpgs)


### plotting venn diagram
library(VennDiagram)

library(VennDiagram)

# Ensure uniqueness
age_associated <- unique(blood_regression_sig_cpgsites$cpg_site)
phylogenetic_age_yrs <- unique(blood_phylo_lambda_sig$cpg_site)
mammalian_univ_clock_1 <- unique(univ_clock_cpgs$age_in_years)

### Working on ven diagram between significant phylogenetic signal, DNAm rate for relative age, and universal DNAm clock 2 (using rel_age)
rel_age_associated <- unique(blood_regression_rel_age_sig_cpgsites$cpg_site)
phylogenetic_rel_age <- unique(blood_phylo_rel_age_lambda_sig$cpg_site)
rel_age_clock_2 <- unique(univ_clock_cpgs$rel_age)


# Create the plot
venn.plot <- draw.triple.venn(
  area1 = length(age_associated),
  area2 = length(phylogenetic_age_yrs),
  area3 = length(mammalian_univ_clock_1),
  n12 = length(intersect(age_associated, phylogenetic_age_yrs)),
  n23 = length(intersect(phylogenetic_age_yrs, mammalian_univ_clock_1)),
  n13 = length(intersect(age_associated, mammalian_univ_clock_1)),
  n123 = length(Reduce(intersect, list(age_associated, phylogenetic_age_yrs,mammalian_univ_clock_1))),
  category = c("age-associated cpgs only", "Phylogenetic (age_in_years) cpgs", "Mammalian universal Clock 1 cpgs"),
  fill = c("pink", "orange", "lightblue"),
  lty = "blank",
  cex = 2,
  cat.cex = 1.5,
  cat.pos = c(-20, 20, 180),
  cat.dist = c(0.05, 0.05, 0.01),
  euler.d = FALSE,
  scaled = FALSE
)

grid.text(
  label = paste0(
    "Total CpG Sites:\n",
    "Age-Associated CpGs: ", length(age_associated), "\n",
    "Phylogenetic CpGs: ", length(phylogenetic), "\n",
    "Mammalian universal CpGs: ", length(mammalian_universal)
  ),
  x = 0.88, y = 0.1, gp = gpar(fontsize = 10)
)


## change colname
univ_clock_cpgs <- univ_clock_cpgs %>%
  rename(age_in_years = clock_1)



### distribution of lambda values from DNAm rate relative_age and age in years 
df1_labeled <- data.frame(lambda = blood_phylo_age_results, source = "phylo_age_years")
df2_labeled <- data.frame(lambda = blood_phylo_result, source = "phylo_rel_age")


df_combined <- rbind(df1_labeled, df2_labeled)

ggplot(df_combined, aes(x = lambda.lambda, fill = source)) +
  geom_histogram(position = "dodge", bins = 30) +
  labs(title = "Blood DNAm Lambda Distribution by datasets source", x = "Lambda", y = "Count") +
  theme_minimal()

