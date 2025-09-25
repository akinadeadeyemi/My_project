######### Author : Adeyemi Akinade
######### Place : Clemson University, SC
######### lab : Gopalan's lab
######### Date : September 2, 2025
######### Project : Comparative Epigenetics of Ageing
######### This script was written to determine the gene enrichment analysis of the cpg sites that show strong 
######### evolutionary correlation between DNAm rate vs LQ

############### GENE ENRICHMENT ANALYSIS ######################
## loading the data
setwd("/Users/aakinad/Downloads/PHD FOLDER/Evolution of Epigenetic Aging across mammalian species /Blood_tissues")
blood_result_LQ <- readRDS("/Users/aakinad/Downloads/PHD FOLDER/Evolution of Epigenetic Aging across mammalian species /Blood_tissues/blood_results_DNAm_rate_LQ_mvBM.rds")

### Visualizing the distribution of the data 
library(ggplot2)

mu <- mean(blood_result_LQ$cor)
sigma <- sd(blood_result_LQ$cor)

ggplot(data.frame(x = blood_result_LQ$cor), aes(x)) +
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
Final_top_cpg_manifest <- Top_cpg_sites %>%
  mutate(Group = case_when(
    IlmnID %in% positive_cpg_LQ$CpG_site ~ "Positive",
    IlmnID %in% negative_cpg_LQ$CpG_site ~ "Negative",
    TRUE ~ "Other"
  ))

region_counts_group <- Final_top_cpg_manifest %>%
  group_by(Group, Human.Hg19_main_Categories) %>%
  summarise(Count = n()) %>%
  ungroup()


### ggplot2 figure
ggplot(region_counts_group, aes(x = Group, y = Count, fill = Human.Hg19_main_Categories)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "CpG Sites", 
       y = "Number of CpGs", 
       fill = "Genomic Region",
       title = expression(paste("Distribution of CpGs (correlation > 0.5 or < -0.5)"))) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c(
    "Exon" = "#1b9e77",
    "Intron" = "#CD1076",
    "Intergenic_downstream" = "#7570b3",
    "Intergenic_upstream" = "#EE3B3B",
    "Promoter" = "#0000CD",
    "fiveUTR" = "#551A8B",
    "threeUTR" = "#66a61e"
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
top_pos_bp <- pos_bp_sig %>%
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
GO_GENE_1000 <- fread ("/Users/aakinad/Downloads/PHD FOLDER/Evolution of Epigenetic Aging across mammalian species /Blood_tissues/Enrichment_Top_1000:least_1000_CpGsites/greatExportAll.tsv")
## selecting top hits
pos_gene_sig <- GO_GENE_1000 %>%
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
       title = "Top enriched GO:GENE terms in top 1000 positive CpGs") +
  theme_minimal(base_size = 16)


############# ENRRICHED CELLULAR PATHWAYS #####################
GO_CC <- fread("/Users/aakinad/Downloads/PHD FOLDER/Evolution of Epigenetic Aging across mammalian species /Blood_tissues/greatExportAll (3).tsv")
GO_BP_1000 <- fread ("/Users/aakinad/Downloads/PHD FOLDER/Evolution of Epigenetic Aging across mammalian species /Blood_tissues/Enrichment_Top_1000:least_1000_CpGsites/greatExportAll (1).tsv")

## selecting top hits
pos_cc_sig <- GO_BP_1000 %>%
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
  labs(x = "GO Biological Pathway", y = "Fold Enrichment", fill = "-log10(HyperP)",
       title = "Top Enriched GO:BP Terms - Positive CpGs") +
  theme_minimal(base_size = 9)

############# ENRRICHED MOLECULAR FUNCTIONS #####################
GO_MF <- fread("/Users/aakinad/Downloads/PHD FOLDER/Evolution of Epigenetic Aging across mammalian species /Blood_tissues/greatExportAll (4).tsv")
GO_CC_1000 <- fread ("/Users/aakinad/Downloads/PHD FOLDER/Evolution of Epigenetic Aging across mammalian species /Blood_tissues/Enrichment_Top_1000:least_1000_CpGsites/greatExportAll (2).tsv")

## selecting top hits
pos_mf_sig <- GO_CC_1000 %>%
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
  labs(x = "GO Cellular Function", y = "Fold Enrichment", fill = "-log10(HyperP)",
       title = "Top Enriched GO:CC Terms - Positive CpGs") +
  theme_minimal(base_size = 12)

######### Molecular Functions ######################
GO_HP <- fread("/Users/aakinad/Downloads/PHD FOLDER/Evolution of Epigenetic Aging across mammalian species /Blood_tissues/greatExportAll (5).tsv")
GO_MF_1000 <- fread ("/Users/aakinad/Downloads/PHD FOLDER/Evolution of Epigenetic Aging across mammalian species /Blood_tissues/Enrichment_Top_1000:least_1000_CpGsites/greatExportAll (3).tsv")

## selecting top hits
pos_hp_sig <- GO_MF_1000 %>%
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
  labs(x = "GO Molecular Function", y = "Fold Enrichment", fill = "-log10(HyperP)",
       title = "Top Enriched GO:MF Terms - CpGs (correlation > 0.5)") +
  theme_minimal(base_size = 10)

######### Molecular Functions ######################
GO_HP_1000 <- fread ("/Users/aakinad/Downloads/PHD FOLDER/Evolution of Epigenetic Aging across mammalian species /Blood_tissues/Enrichment_Top_1000:least_1000_CpGsites/greatExportAll (4).tsv")
## selecting top hits
pos_hph_sig <- GO_HP_1000 %>%
  filter(HyperFdrQ < 0.05) %>%
  arrange(HyperFdrQ)

# For plotting: top 20 by fold enrichment
top_pos_hph <- pos_hph_sig %>%
  arrange(desc(RegionFoldEnrich)) %>%
  mutate(logP = -log10(HyperP))

# Positive CpGs
library(ggplot2)
ggplot(top_pos_hph, aes(x = reorder(Desc, RegionFoldEnrich), y = RegionFoldEnrich, fill = RegionFoldEnrich)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "GO Human Phenotype", y = "Fold Enrichment", fill = "-log10(HyperP)",
       title = "Top Enriched GO:HP Terms - CpGs (correlation > 0.5)") +
  theme_minimal(base_size = 12)


######################################################################################
################# NEGATIVELY CORRELATED CPGS #########################################
GO_GENE_NEG <- fread("/Users/aakinad/Downloads/PHD FOLDER/Evolution of Epigenetic Aging across mammalian species /Blood_tissues/greatExportAll (6).tsv")
GO_GENE_NEG_1000 <- fread ("/Users/aakinad/Downloads/PHD FOLDER/Evolution of Epigenetic Aging across mammalian species /Blood_tissues/Enrichment_Top_1000:least_1000_CpGsites/bottom_1000_cpgsites/greatExportAll.tsv")

## selecting top hits
neg_gene_sig <- GO_GENE_NEG_1000 %>%
  filter(HyperFdrQ < 0.05) %>%
  arrange(HyperFdrQ)

# For plotting: top 20 by fold enrichment
top_neg_gene <- neg_gene_sig %>%
  arrange(desc(RegionFoldEnrich)) %>%
  mutate(logP = -log10(HyperP))

# Positive CpGs
library(ggplot2)
ggplot(top_neg_gene, aes(x = reorder(Desc, RegionFoldEnrich), y = RegionFoldEnrich, fill = logP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "GO Ensembl Genes", y = "Fold Enrichment", fill = "logP",
       title = "Top Enriched GO:GENE Terms - CpGs (correlation < -0.5)") +
  theme_minimal(base_size = 14)

############ Cellular Component ##########################
GO_CC_NEG_1000 <- fread ("/Users/aakinad/Downloads/PHD FOLDER/Evolution of Epigenetic Aging across mammalian species /Blood_tissues/Enrichment_Top_1000:least_1000_CpGsites/bottom_1000_cpgsites/greatExportAll (1).tsv")
neg_cc_sig <- GO_CC_NEG_1000 %>%
  filter(HyperFdrQ < 0.05) %>%
  arrange(HyperFdrQ)

# For plotting: top 20 by fold enrichment
top_neg_cc <- neg_cc_sig %>%
  arrange(desc(RegionFoldEnrich)) %>%
  mutate(logP = -log10(HyperP))

# Positive CpGs
library(ggplot2)
ggplot(top_neg_cc, aes(x = reorder(Desc, RegionFoldEnrich), y = RegionFoldEnrich, fill = logP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "GO Cellular Component", y = "Fold Enrichment", fill = "logP",
       title = "Top Enriched GO:CC Terms - CpGs (correlation < -0.5)") +
  theme_minimal(base_size = 14)


############ Cellular Component ##########################
GO_MF_NEG_1000 <- fread ("/Users/aakinad/Downloads/PHD FOLDER/Evolution of Epigenetic Aging across mammalian species /Blood_tissues/Enrichment_Top_1000:least_1000_CpGsites/bottom_1000_cpgsites/greatExportAll (2).tsv")
neg_cc_sig <- GO_CC_NEG_1000 %>%
  filter(HyperFdrQ < 0.05) %>%
  arrange(HyperFdrQ)

# For plotting: top 20 by fold enrichment
top_neg_cc <- neg_cc_sig %>%
  arrange(desc(RegionFoldEnrich)) %>%
  mutate(logP = -log10(HyperP))

# Positive CpGs
library(ggplot2)
ggplot(top_neg_cc, aes(x = reorder(Desc, RegionFoldEnrich), y = RegionFoldEnrich, fill = logP)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(x = "GO Cellular Component", y = "Fold Enrichment", fill = "logP",
       title = "Top Enriched GO:CC Terms - CpGs (correlation < -0.5)") +
  theme_minimal(base_size = 14)

#######################Selecting for just top 100 CPGs ##########################
top_2000 <- blood_result_LQ %>% 
  arrange(desc(cor)) %>% 
  slice(1:2000)

### bottom 100
bottom_2000 <- blood_result_LQ %>% 
  arrange(cor) %>% 
  slice(1:2000)

## extracting top cpg site in LQ/DNAm rate covariance analysis
Top_2000_sites <- manifest[manifest$IlmnID %in% top_2000$CpG_site, ]  %>% 
  select ("Human.Hg38_seqnames", "Human.Hg38_start", "Human.Hg38_end", "IlmnID")

bottom_2000_sites <- manifest[manifest$IlmnID %in% bottom_2000$CpG_site, ] %>% 
  select ("Human.Hg38_seqnames", "Human.Hg38_start", "Human.Hg38_end", "IlmnID")
## Adding header 
colnames(Top_2000_sites) <- c("chr", "start", "end", "name")
colnames(bottom_2000_sites) <- c("chr", "start", "end", "name")


### converting the extracted table into bedfile/
write.table(
Top_2000_sites, 
  file = "Top_2000_sites.bed", 
  sep = "\t", 
  quote = FALSE, 
  row.names = FALSE, 
  col.names = FALSE
)

