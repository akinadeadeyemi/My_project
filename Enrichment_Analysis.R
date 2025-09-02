######### Author : Adeyemi Akinade
######### Place : Clemson University, SC
######### lab : Gopalan's lab
######### Date : September 2, 2025
######### Project : Comparative Epigenetics of Ageing
######### This script was written to determine the gene enrichment analysis of the 

############### GENE ENRICHMENT ANALYSIS ######################
## loading the data
blood_result_LQ <- readRDS("/Users/aakinad/Downloads/PHD FOLDER/Evolution of Epigenetic Aging across mammalian species /Blood_tissues/blood_results_DNAm_rate_LQ_mvBM.rds")

### Visualizing the distribution of the data 
library(ggplot2)

mu <- mean(result_LQ$cor)
sigma <- sd(result_LQ$cor)

ggplot(data.frame(x = result_LQ$cor), aes(x)) +
  geom_histogram(aes(y = ..density..), bins = 25, fill = "lightgray", color = "black") +
  stat_function(fun = dnorm, args = list(mean = mu, sd = sigma), 
                color = "red", size = 1.0)

## selecting top cpg sites for enrichment analysis
library(dplyr)
Top_cpg_LQ <- blood_result_LQ %>%
  filter(abs(cor) > 0.5)


### loading and unzipping the mammalian DNAm manifest file
manifest <- read.csv("/Users/aakinad/Downloads/PHD FOLDER/Manifest_chromosome.csv")
View(manifest)

## extracting top cpg site in LQ/DNAm rate covariance analysis
Top_cpg_sites <- manifest[manifest$IlmnID %in% Top_cpg_LQ$CpG_site, ]

top_site_bedfile <- Top_cpg_sites %>% select ("IlmnID", "Human.Hg38_seqnames", "Human.Hg38_end", "Human.Hg38_start")

### converting the extracted table into bedfile/
write.table(
  top_site_bedfile, 
  file = "top_cpgsite.bed", 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE, 
  col.names = FALSE
)





















