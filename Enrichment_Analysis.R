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
## extract  cpgsites with evolutionary positive correlation between DNAm rate vs LQ --> corr = or > 0.5
positive_cpg_LQ <- blood_result_LQ %>%
  filter((cor) >= 0.5)
## Evolutionary negative correlation
negative_cpg_LQ <- blood_result_LQ %>%
  filter((cor) <= -0.5)

### loading and unzipping the mammalian DNAm manifest file
manifest <- read.csv("/Users/aakinad/Downloads/PHD FOLDER/Manifest_chromosome.csv")
View(manifest)

## extracting top cpg site in LQ/DNAm rate covariance analysis
Top_cpg_sites <- manifest[manifest$IlmnID %in% Top_cpg_LQ$CpG_site, ]
positive_cpg_bed <- manifest[manifest$IlmnID %in% positive_cpg_LQ$CpG_site, ]
negative_cpg_bed <- manifest[manifest$IlmnID %in% negative_cpg_LQ$CpG_site, ]

## selecting the columns with important columns
top_site_bedfile <- Top_cpg_sites %>% select ("Human.Hg38_seqnames", "Human.Hg38_start", "Human.Hg38_end", "IlmnID")
positive_bed <- positive_cpg_bed %>% select ("Human.Hg38_seqnames", "Human.Hg38_start", "Human.Hg38_end", "IlmnID")
negative_bed <- negative_cpg_bed %>% select ("Human.Hg38_seqnames", "Human.Hg38_start", "Human.Hg38_end", "IlmnID")

## Adding header 
colnames(positive_bed) <- c("chr", "start", "end", "name")
colnames(negative_bed) <- c("chr", "start", "end", "name")


### converting the extracted table into bedfile/
write.table(
  negative_bed, 
  file = "negative_cpg.bed", 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE, 
  col.names = FALSE
)

## creating a background bedfile for GREAT Enrichment Analysis
background_cpgsites <- manifest %>% select ("Human.Hg38_seqnames", "Human.Hg38_start", "Human.Hg38_end", "IlmnID") %>% slice(1:37554)

background_cpgs <- background_cpgsites %>% 
  filter(!is.na(Human.Hg38_seqnames) & !is.na(Human.Hg38_start) & !is.na(Human.Hg38_end) & !is.na(IlmnID))

write.table(
  background_cpgs, 
  file = "background_cpgs.bed", 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE, 
  col.names = FALSE
)

############# VISUALIZING THE LOCATION OF BOTH POSITIVE AND NEGATIVE CPGS #########################
### labelling cpgs that are evolutionary positive and negative correlated beetween DNAm rate with LQ
manifest <- manifest %>%
  mutate(Group = case_when(
    IlmnID %in% positive_cpgs ~ "Positive",
    IlmnID %in% negative_cpgs ~ "Negative",
    TRUE ~ "Other"
  ))














######### VISUALIZING THE OUTPUT OF GREAT #############
############## GO: BIOLOGICAL PROCESS #################
library(readr)
install.packages("data.table")
library(data.table)

GO_BP <- fread("/Users/aakinad/Downloads/PHD FOLDER/Evolution of Epigenetic Aging across mammalian species /Blood_tissues/greatExportAll.tsv")

### extracting top hits using HyperFDRQ < 0.05 for stringency.

# Positive CpGs
library(dplyr)
pos_bp_sig <- GO_BP %>%
  filter(HyperFdrQ < 0.05) %>%
  arrange(HyperFdrQ)

# For plotting: top 20 by fold enrichment
top_pos_bp <- pos_sig %>%
  arrange(desc(RegionFoldEnrich)) %>%
  head(32) %>%
  mutate(logP = -log10(HyperP))

# Positive CpGs
library(ggplot2)
ggplot(top_pos_bp, aes(x = reorder(Desc, RegionFoldEnrich), y = RegionFoldEnrich, fill = logP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "GO Biological Process", y = "Fold Enrichment", fill = "-log10(HyperP)",
       title = "Top Enriched GO:BP Terms - Positive CpGs") +
  theme_minimal(base_size = 10)


############# ENRRICHED GENES #####################
GO_GENE <- fread("/Users/aakinad/Downloads/PHD FOLDER/Evolution of Epigenetic Aging across mammalian species /Blood_tissues/greatExportAll (1).tsv")

## selecting top hits
pos_gene_sig <- GO_GENE %>%
  filter(HyperFdrQ < 0.05) %>%
  arrange(HyperFdrQ)

# For plotting: top 20 by fold enrichment
top_pos_gene <- pos_gene_sig %>%
  arrange(desc(RegionFoldEnrich)) %>%
  mutate(logP = -log10(HyperP))

# Positive CpGs
library(ggplot2)
ggplot(top_pos_gene, aes(x = reorder(Desc, RegionFoldEnrich), y = RegionFoldEnrich, fill = logP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "GO Ensembl Genes", y = "Fold Enrichment", fill = "logP",
       title = "Top Enriched GO:GENE Terms - Positive CpGs") +
  theme_minimal(base_size = 10)


############# ENRRICHED CELLULAR PATHWAYS #####################
GO_CC <- fread("/Users/aakinad/Downloads/PHD FOLDER/Evolution of Epigenetic Aging across mammalian species /Blood_tissues/greatExportAll (3).tsv")

## selecting top hits
pos_cc_sig <- GO_CC %>%
  filter(HyperFdrQ < 0.05) %>%
  arrange(HyperFdrQ)

# For plotting: top 20 by fold enrichment
top_pos_cc <- pos_cc_sig %>%
  arrange(desc(RegionFoldEnrich)) %>%
  mutate(logP = -log10(HyperP))

# Positive CpGs
library(ggplot2)
ggplot(top_pos_cc, aes(x = reorder(Desc, RegionFoldEnrich), y = RegionFoldEnrich, fill = RegionFoldEnrich)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "GO Cellular Component", y = "Fold Enrichment", fill = "-log10(HyperP)",
       title = "Top Enriched GO:CC Terms - Positive CpGs") +
  theme_minimal(base_size = 12)

############# ENRRICHED MOLECULAR FUNCTIONS #####################
GO_MF <- fread("/Users/aakinad/Downloads/PHD FOLDER/Evolution of Epigenetic Aging across mammalian species /Blood_tissues/greatExportAll (4).tsv")

## selecting top hits
pos_mf_sig <- GO_MF %>%
  filter(HyperFdrQ < 0.05) %>%
  arrange(HyperFdrQ)

# For plotting: top 20 by fold enrichment
top_pos_mf <- pos_mf_sig %>%
  arrange(desc(RegionFoldEnrich)) %>%
  mutate(logP = -log10(HyperP))

# Positive CpGs
library(ggplot2)
ggplot(top_pos_mf, aes(x = reorder(Desc, RegionFoldEnrich), y = RegionFoldEnrich, fill = RegionFoldEnrich)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "GO Molecular Function", y = "Fold Enrichment", fill = "-log10(HyperP)",
       title = "Top Enriched GO:MF Terms - Positive CpGs") +
  theme_minimal(base_size = 10)

######### HUMAN PHENOTYPES ######################
GO_HP <- fread("/Users/aakinad/Downloads/PHD FOLDER/Evolution of Epigenetic Aging across mammalian species /Blood_tissues/greatExportAll (5).tsv")

## selecting top hits
pos_hp_sig <- GO_HP %>%
  filter(HyperFdrQ < 0.05) %>%
  arrange(HyperFdrQ)

# For plotting: top 20 by fold enrichment
top_pos_hp <- pos_hp_sig %>%
  arrange(desc(RegionFoldEnrich)) %>%
  mutate(logP = -log10(HyperP))

# Positive CpGs
library(ggplot2)
ggplot(top_pos_hp, aes(x = reorder(Desc, RegionFoldEnrich), y = RegionFoldEnrich, fill = RegionFoldEnrich)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "GO Human Phenotype", y = "Fold Enrichment", fill = "-log10(HyperP)",
       title = "Top Enriched GO:HP Terms - Positive CpGs") +
  theme_minimal(base_size = 10)





