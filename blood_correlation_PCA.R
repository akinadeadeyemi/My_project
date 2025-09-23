######### Author : Adeyemi Akinade
######### Place : Clemson University, SC
######### lab : Gopalan's lab
######### Date : March 24, 2025
######### Project : Comparative Epigenetics of Ageing
######### This script was written to determine the PCA of age-correlated cpg sites in blood and skin tissues


######################## PCA OF AGE-CORRELATED CPG SITES FROM BLOOD SAMPLES ONLY ##########################################
## loading the data...
blood_correlation_output <- readRDS("~/Phd_data/outputs/blood_correlation_full_sample.rds")

library(dplyr)
library(tidyr)
#### shape the data into a wide format 
blood_correlation_wide_data <- blood_correlation_output %>%
  dplyr::select(species, cpg_site, correlation) %>%
  pivot_wider(names_from = cpg_site, values_from = correlation, values_fill = 0)

# View the reshaped data
head(blood_correlation_wide_data)

# Set species as row names
blood_correlation_wide_data <- blood_correlation_wide_data %>%
  column_to_rownames(var = "species")

View(blood_correlation_wide_data)

############ PCA analysis for the raw beta values in blood samples ##############
# Running PCA 
blood_correlation_pca <- prcomp(blood_correlation_wide_data, scale. = TRUE)


# View the principal component scores (sample-level PCA scores)
head(blood_correlation_pca$x)


# View the proportion of variance explained by each principal component
blood_correlation_pve <- (blood_correlation_pca$sdev^2) / sum(blood_correlation_pca$sdev^2)

blood_correlation_proportion <- blood_correlation_pve * 100

######## Generating a Scree plot ########
# Create a data frame for plotting
blood_correlation_scree_data <- data.frame(
  PC = 1:length(blood_correlation_pve),
  Variance = blood_correlation_pve,
  proportion = blood_correlation_proportion
)

# selecting the first 20 PCs
blood_correlation_PC20 <- blood_correlation_scree_data[1:20, ]




## Visualizing the percentage variance explained by each PCs in blood samples only

ggplot(blood_correlation_PC20, aes(x = PC, y = proportion)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(1, 20, by = 1)) + 
  scale_y_continuous(breaks = seq(2, 20, by = 2), limits = c(0, 20)) +  
  labs(
    title = "Percentage variance explained in the first 20 PCs of correlation between age and cpg_sites from blood samples",
    x = "Principal Component (PCs)",
    y = "Proportion of Variance Explained (%)"
  ) + 
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = "transparent", color = 1, size = 1)
  ) +
  theme(plot.title = element_text(family = "serif",              # Font family
                                  face = "bold",                 # Font face
                                  color = 1,                     # Font color
                                  size = 10,                     # Font size
                                  hjust = 0.5,                     # Horizontal adjustment
                                  vjust = 2.5,                     # Vertical adjustment
                                  angle = 0,                   # Font angle
                                  lineheight = 1,                # Line spacing
                                  margin = margin(20, 0, 0, 0)))

                                                          


# Extract the scores for plotting
blood_correlation_pca_scores <- as.data.frame(blood_correlation_pca$x)
head(blood_correlation_pca_scores)


# Add sample IDs as a new column
blood_correlation_pca_scores$species <- rownames(blood_correlation_pca_scores)

### Specify the species to be labelled
blood_correlation_pca_scores <- blood_correlation_pca_scores %>%
  mutate(label = ifelse(species %in% c("Acinonyx jubatus", "Aonyx cinereus", "Apodemus sylvaticus", "Bos taurus",            
                                        "Callithrix jacchus", "Canis lupus familiaris", "Capreolus capreolus", "Ceratotherium simum simum",  
                                        "Cervus canadensis" , "Cervus elaphus", "Chlorocebus sabaeus", "Delphinapterus leucas"      ,
                                        "Elephas maximus", "Equus caballus", "Equus quagga", "Felis catus"  ,
                                        "Heterocephalus glaber", "Homo sapiens", "Lagenorhynchus obliquidens", "Loxodonta africana",        
                                        "Macaca mulatta", "Macropus giganteus", "Marmota flaviventris", "Mus musculus", 
                                      "Odobenus rosmarus divergens", "Orcinus orca", "Osphranter rufus", "Ovis aries",               
                                       "Phoca vitulina", "Rattus norvegicus", "Sus scrofa", "Sus scrofa domesticus",   
                                       "Tursiops truncatus", "Zalophus californianus"
                                       ), species, ""))

library(ggrepel)
# Plot PCA with grouping by order in blood samples only

custom_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF00FF", "#838B8B", "#0000FF", "#050505","#00CD00","#FF1493", "#FFC1C1","#B03060",
  "#FF7F00", "#A65628", "#00F5FF", "#F781BF", "#CD8500", "#66C2A5", "#7CFC00", "#0000CD","#8B2323","#6B8E23", "#87CEFF", "#A4D3EE",
  "#FC8D62", "#458B74", "#CDBA96", "#8B3626", "#006400", "#FFD700", "#8B6508","#00008B", "#9370DB","#96CDCD"  
)

ggplot(blood_correlation_pca_scores, aes(x = PC3, y = PC4, color = species)) +
  geom_point(size = 2.0, alpha = 1) +
  scale_color_manual(values = custom_colors) + # Use your custom colors
  labs(
    title = "PCA of correlation of age with cpg_sites of different species from blood samples",
    x = paste0("Principal Component 3 (", round(blood_correlation_proportion[3], 2), "%)"),
    y = paste0("Principal Component 4 (", round(blood_correlation_proportion[4], 2), "%)")
  ) + geom_text_repel(aes(label = label), size = 2, max.overlaps = 50) +
  #geom_text(aes(label = label), vjust = -1, hjust = 2, size = 3, check_overlap = TRUE) +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(panel.border = element_rect(fill = "transparent", 
                                    color = 1,            
                                    size = 0.5)) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(family = "serif",              # Font family
                                  face = "bold",                 # Font face
                                  color = 1,                     # Font color
                                  size = 10,                     # Font size
                                  hjust = 0.5,                     # Horizontal adjustment
                                  vjust = 2.5,                     # Vertical adjustment
                                  angle = 0,                   # Font angle
                                  lineheight = 1,                # Line spacing
                                  margin = margin(20, 0, 0, 0)))





