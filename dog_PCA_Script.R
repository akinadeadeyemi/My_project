#################### PCA for canis lupus familiaris blood samples ######################
####################                                              ######################

#### extract the sample IDs
dog_sample_id <- blood_samples %>% filter(species == "Canis lupus familiaris") %>% select (geo_accession, sample_barcode)
View(dog_sample_id)

### Selcting dog samples
dog_sample <- blood_samples %>% filter(species == "Canis lupus familiaris")

# Create a named vector for mapping
dog_id_mapping <- setNames(dog_sample_id$geo_accession, dog_sample_id$sample_barcode)
View(dog_id_mapping)

# Replace column names
new_column_names <- ifelse(colnames(merged_betas) %in% names(dog_id_mapping),
                           dog_id_mapping[colnames(merged_betas)],
                           colnames(merged_betas))

# Assign new column names
colnames(merged_betas) <- new_column_names
View(merged_betas) 

####### merging sample IDs in both beta and blood_skin_samples to generate  ########
dog_betas <- merged_betas[, colnames(merged_betas) %in% dog_sample_id$geo_accession]
View(dog_betas)

#### starting PCA  ######
############ PCA analysis for the raw beta values ##############
dog_beta_transpose <- t(dog_betas)


# Running PCA 
dog_pca <- prcomp(dog_beta_transpose, scale. = TRUE)

###### Check the structure of the PCA result
str(dog_pca)

# View the principal component scores (sample-level PCA scores)
head(dog_pca)


# View the proportion of variance explained by each principal component
dog_variance <- (dog_pca$sdev^2) / sum(dog_pca$sdev^2)
dog_explained<- dog_variance * 100

######## Generating a Scree plot ########
library(ggplot2)

# Create a data frame for plotting
dog_scree <- data.frame(
  PC = 1:length(dog_variance),
  Variance = dog_variance,
  proportion = dog_explained
)
# selecting the first 20 PCs
dog_20_PC <- dog_scree[1:20, ]


## Visualizing the percentage variance explained by each PCs


ggplot(dog_20_PC, aes(x = PC, y = proportion)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(1, 20, by = 1)) + 
  scale_y_continuous(breaks = seq(2, 24, by = 2), limits = c(0, 25)) +  
  labs(title = " Percentage variance explained in the First 20 PCs of dog blood samples",
       x = "Principal Component (PCs)",
       y = "Proportion of Variance Explained") +
  theme_minimal()




# Extract the scores for plotting
dog_pca_scores <- as.data.frame(dog_pca$x)
str(dog_pca_scores)
View(dog_pca_scores)

# Add sample IDs as a new column
dog_pca_scores$geo_accession <- rownames(dog_beta_transpose)

#### Visualizing PCA 1 and 2 by merging PCA scores with sample metadata
dog_pca_scores <- dog_pca_scores %>%
  merge(dog_sample, by='geo_accession', all.x=TRUE) %>% drop_na()

View(dog_pca_scores)


# Plot PCA with grouping by species
ggplot(dog_pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(size = 1.0, alpha = 1, col = "navy") +
  labs(
    title = "PCA of Beta values by tissues",
    x = paste0("Principal Component 1 (", round(proportion_explained[1], 2), "%)"),
    y = paste0("Principal Component 2 (", round(proportion_explained[2], 2), "%)"),
  ) +
  theme(panel.grid = element_blank()) + theme(legend.title = element_text(family = "Roboto",
                                                                          color = "black", 
                                                                          size = 12,
                                                                          face = 2)) + 
  theme(panel.border = element_rect(fill = "transparent", 
                                    color = 1,            
                                    size = 0.5)) + theme_minimal()

