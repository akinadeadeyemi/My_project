######################## PCA OF AGE-CORRELATED CPG SITES FROM BLOOD SAMPLES ONLY ##########################################
## loading the data...
skin_correlation_output <- readRDS("~/Phd_data/outputs/skin_correlation_full_sample_output.rds")

library(dplyr)
#### shape the data into a wide format 
skin_correlation_wide_data <- skin_correlation_output %>%
  select(species, cpg_site, correlation) %>%
  pivot_wider(names_from = cpg_site, values_from = correlation, values_fill = 0)

# View the reshaped data
View(skin_correlation_wide_data)

# Set species as row names
skin_correlation_wide_data <- skin_correlation_wide_data %>%
  column_to_rownames(var = "species")

View(skin_correlation_wide_data)

############ PCA analysis for the raw beta values in blood samples ##############
# Running PCA 
skin_correlation_pca <- prcomp(skin_correlation_wide_data, scale. = TRUE)


# View the principal component scores (sample-level PCA scores)
head(skin_correlation_pca$x)


# View the proportion of variance explained by each principal component
skin_correlation_pve <- (skin_correlation_pca$sdev^2) / sum(skin_correlation_pca$sdev^2)

skin_correlation_proportion <- skin_correlation_pve * 100

######## Generating a Scree plot ########
# Create a data frame for plotting
skin_correlation_scree_data <- data.frame(
  PC = 1:length(skin_correlation_pve),
  Variance = skin_correlation_pve,
  proportion = skin_correlation_proportion
)
View(skin_correlation_scree_data)

# selecting the first 20 PCs
skin_correlation_PC20 <- skin_correlation_scree_data[1:20, ]




## Visualizing the percentage variance explained by each PCs in blood samples only

ggplot(skin_correlation_PC20, aes(x = PC, y = proportion)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(1, 20, by = 1)) + 
  scale_y_continuous(breaks = seq(2, 20, by = 2), limits = c(0, 20)) +  
  labs(
    title = "Percentage variance explained in the first 20 PCs of correlation between age and cpg_sites from skin samples",
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
skin_correlation_pca_scores <- as.data.frame(skin_correlation_pca$x)
head(skin_correlation_pca_scores)


# Add sample IDs as a new column
skin_correlation_pca_scores$species <- rownames(skin_correlation_pca_scores)

### Specify the species to be labelled
skin_correlation_pca_scores <- skin_correlation_pca_scores %>%
  mutate(label = ifelse(species %in% c("Antrozous pallidus", "Bos taurus",  "Carollia perspicillata",        
                                        "Cephalorhynchus hectori hectori", "Cynopterus brachyotis", "Delphinapterus leucas",          
                                        "Desmodus rotundus", "Eidolon helvum", "Eptesicus fuscus",               
                                       "Equus quagga", "Heterocephalus glaber", "Homo sapiens",                  
                                        "Lagenorhynchus obliquidens", "Leptonycteris yerbabuenae", "Macaca mulatta",                 
                                        "Megaptera novaeangliae", "Molossus molossus", "Mus musculus",                 
                                        "Myotis lucifugus", "Myotis myotis" , "Myotis vivesi",               
                                        "Orcinus orca", "Oreamnos americanus" , "Phoca vitulina",               
                                        "Phyllostomus discolor", "Phyllostomus hastatus",  "Pteropus hypomelanus",           
                                        "Pteropus poliocephalus", "Pteropus pumilus", "Pteropus rodricensis",           
                                        "Pteropus vampyrus", "Rattus norvegicus", "Rhinolophus ferrumequinum",      
                                       "Rhynchonycteris naso", "Rousettus aegyptiacus" , "Saccopteryx bilineata",          
                                        "Tadarida brasiliensis" , "Tursiops aduncus", "Tursiops truncatus",           
                                        "Zalophus californianus"    
  ), species, ""))

library(ggrepel)
# Plot PCA with grouping by order in blood samples only

custom_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF00FF", "#838B8B", "#0000FF", "#050505","#00CD00","#FF1493", "#FFC1C1","#B03060", "#2F4F4F","#551A8B",
  "#FF7F00", "#A65628", "#00F5FF", "#F781BF", "#CD8500", "#66C2A5", "#7CFC00", "#0000CD","#8B2323","#6B8E23", "#87CEFF", "#A4D3EE", "#CD2990", 
  "#FC8D62", "#458B74", "#CDBA96", "#8B3626", "#006400", "#FFD700", "#8B6508","#00008B", "#9370DB","#96CDCD", "#8B4500", "#2F4F4F",  "#87CEFA"
)

ggplot(skin_correlation_pca_scores, aes(x = PC3, y = PC4, color = species)) +
  geom_point(size = 2.0, alpha = 1) +
  scale_color_manual(values = custom_colors) + # Use your custom colors
  labs(
    title = "PCA of correlation of age with cpg_sites of different species from blood samples",
    x = paste0("Principal Component 3 (", round(skin_correlation_proportion[3], 2), "%)"),
    y = paste0("Principal Component 4 (", round(skin_correlation_proportion[4], 2), "%)")
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


