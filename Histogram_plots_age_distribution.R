######################## Histogram of age distribution from blood samples of mammals #######################
####. Author:Adeyemi Akinade
####  Place : Clemson University, SC
####  lab   : Gopalan's lab
####  Project:  Comparative Epigenetics of Aging
####  Date   : 29 January, 2025



### Histogram of age distribution in blood samples
### Retrieving samples from blood samples extracted from the sample metadata ###
blood_samples <- blood_samples %>% mutate(species= ifelse(species == "Sus scrofa domesticus", "Sus scrofa", species)) # changing sus scrofa domesticus

### installing (ggpubr)
install.packages("ggpubr")
library(ggpubr)

## Plotting histogram
canis <- blood_samples %>% filter(species == "Canis lupus familiaris") %>% select(age_years) %>% drop_na()
hist_canis <- ggplot(canis, aes(x = age_years)) + 
  geom_histogram(
                 colour = 1, fill = "white") +
  geom_vline(xintercept = 1.83, color = "navy", linetype = "dashed", size = 0.5) +
  labs(title = "Histogram of canis lupus age in years (blood samples)", x = "Age in years", y = "Count") +
  geom_text(aes(x = 1.82 + 3.8, y = 60, label = "Age at sexual maturity(1.83)"), color = "black", angle = 0, vjust = 1) +
  theme_minimal()
hist_canis

## Histogram of Mus musculus 
musculus <- blood_samples %>% filter(species == "Mus musculus") %>% select(age_years) %>% drop_na()

musculus_hist <- ggplot(musculus, aes(x = age_years)) + 
    geom_histogram(
      colour = 1, fill = "white") +
    geom_vline(xintercept = 0.12, color = "navy", linetype = "dashed", size = 0.5) +
    labs(title = "Histogram of Mus musculus age in years (blood samples)", x = "Age in years", y = "Count") +
    geom_text(aes(x = 0.12 + 0.65, y = 45, label = "Age at sexual maturity (0.12)"), color = "black", angle = 0, vjust = 1) +
    theme_minimal()

musculus_hist
  
## Histogram of Homo sapiens 
homo <- blood_samples %>% filter(species == "Homo sapiens") %>% select(age_years) %>% drop_na()
  
homo_hist <- ggplot(homo, aes(x = age_years)) + 
    geom_histogram(
      colour = 1, fill = "white") +
    geom_vline(xintercept = 13.5, color = "navy", linetype = "dashed", size = 0.5) +
    labs(title = "Histogram of Homo sapiens Age in years (blood samples)", x = "Age in years", y = "Count") +
  scale_x_continuous(breaks = seq(2, 120, by = 4)) + scale_y_continuous(breaks = seq(0, 90, by = 10)) +
    geom_text(aes(x = 13.5 + 15, y = 90, label = "Age at sexual maturity (13.5)"), color = "black", angle = 0, vjust = 1) +
  theme_minimal()  
  
homo_hist  


## Histogram of Bos taurus 
bos <- blood_samples %>% filter(species == "Bos taurus") %>% select(age_years) %>% drop_na()

bos_hist <- ggplot(bos, aes(x = age_years)) + 
  geom_histogram(
    colour = 1, fill = "white") +
  geom_vline(xintercept = 1.5, color = "navy", linetype = "dashed", size = 0.5) +
  labs(title = "Histogram of Bos taurus age in years (blood samples)", x = "Age in years", y = "Count") +
  scale_x_continuous(breaks = seq(2, 20, by = 2)) + scale_x_continuous(breaks = seq(0, 50, by = 10)) +
  geom_text(aes(x = 1.5 + 2.4, y = 40, label = "Age at sexual maturity (1.5)"), color = "black", angle = 0, vjust = 1) +
  theme_minimal()  

bos_hist


## Histogram of Equus caballus 
caballus <- blood_samples %>% filter(species == "Equus caballus") %>% select(age_years) %>% drop_na()

caballus_hist <- ggplot(caballus, aes(x = age_years)) + 
  geom_histogram(
    colour = 1, fill = "white") +
  geom_vline(xintercept = 2.58, color = "navy", linetype = "dashed", size = 0.5) +
  labs(title = "Histogram of Equus caballus age in years (blood samples)", x = "Age in years", y = "Count") +
  scale_x_continuous(breaks = seq(2, 120, by = 4)) +
  geom_text(aes(x = 2.58 + 5.5, y = 25, label = "Age at sexual maturity (2.58)"), color = "black", angle = 0, vjust = 1) +
  theme_minimal()  

caballus_hist


## Histogram of Macaca mulatta 
macaca <- blood_samples %>% filter(species == "Macaca mulatta") %>% select(age_years) %>% drop_na()

macaca_hist <- ggplot(macaca, aes(x = age_years)) + 
  geom_histogram(
    colour = 1, fill = "white") +
  geom_vline(xintercept = 4.44, color = "navy", linetype = "dashed", size = 0.5) +
  labs(title = "Histogram of Macaca mulatta age in years (blood samples)", x = "Age in years", y = "Count") +
  scale_x_continuous(breaks = seq(5, 30, by = 5)) +
  geom_text(aes(x = 4.44 + 8.0, y = 30, label = "Age at sexual maturity (4.44)"), color = "black", angle = 0, vjust = 1) +
  theme_minimal()  

macaca_hist

############## combined plots. #################################
library(ggpubr)
combined_plot_1 <- ggarrange(hist_canis, bos_hist, musculus_hist, 
                             homo_hist, caballus_hist, macaca_hist, 
                             nrow = 3,
                             ncol = 2)
combined_plot_1
###############################################################



## Histogram of Capreolus capreolus
capreolus <- blood_samples %>% filter(species == "Capreolus capreolus") %>% select(age_years) %>% drop_na()

capreolus_hist <- ggplot(capreolus, aes(x = age_years)) + 
  geom_histogram(
    colour = 1, fill = "white") +
  geom_vline(xintercept = 1.46, color = "navy", linetype = "dashed", size = 0.5) +
  labs(title = "Histogram of Capreolus capreolus age in years (blood samples)", x = "Age in years", y = "Count") +
  scale_x_continuous(breaks = seq(1, 15, by = 1)) + scale_y_continuous(breaks = seq(0, 60, by = 5)) +
  geom_text(aes(x = 1.46 + 2.8, y = 55, label = "Age at sexual maturity (1.46)"), color = "black", angle = 0, vjust = 1) +
  theme_minimal()  

capreolus_hist



## Histogram of Rattus norvegicus
rattus <- blood_samples %>% filter(species == "Rattus norvegicus") %>% select(age_years) %>% drop_na()

rattus_hist <- ggplot(rattus, aes(x = age_years)) + 
  geom_histogram(
    colour = 1, fill = "white") +
  geom_vline(xintercept = 0.22, color = "navy", linetype = "dashed", size = 0.5) +
  labs(title = "Histogram of Rattus norvegicus age in years (blood samples)", x = "Age in years", y = "Count") +
  scale_x_continuous(breaks = seq(0.1, 5, by = 0.2)) + scale_y_continuous(breaks = seq(0, 35, by = 5)) +
  geom_text(aes(x = 0.2 + 0.45, y = 35, label = "Age at maturity (0.22)"), color = "black", angle = 0, vjust = 1) +
  theme_minimal()  

rattus_hist



## Histogram of Ovis aries
ovis <- blood_samples %>% filter(species == "Ovis aries") %>% select(age_years) %>% drop_na()

ovis_hist <- ggplot(ovis, aes(x = age_years)) + 
  geom_histogram(
    colour = 1, fill = "white") +
  geom_vline(xintercept = 2, color = "navy", linetype = "dashed", size = 0.5) +
  labs(title = "Histogram of Ovis aries age in years (blood samples)", x = "Age in years", y = "Count") +
  scale_x_continuous(breaks = seq(1, 15, by = 1)) +
  scale_y_continuous(breaks = seq(0, 80, by = 10)) +
  geom_text(aes(x = 2 + 1.2, y = 60, label = "Age at sexual maturity (2)"), color = "black", angle = 0, vjust = 1) +
  theme_minimal()  

ovis_hist



## Histogram of Marmota flaviventris
marmota <- blood_samples %>% filter(species == "Marmota flaviventris") %>% select(age_years) %>% drop_na()

marmota_hist <- ggplot(marmota, aes(x = age_years)) + 
  geom_histogram(
    colour = 1, fill = "white") +
  geom_vline(xintercept = 2, color = "navy", linetype = "dashed", size = 0.5) +
  labs(title = "Histogram of Marmota flaviventris age in years (blood samples)", x = "Age in years", y = "Count") +
  scale_x_continuous(breaks = seq(0, 15, by = 2)) + scale_y_continuous(breaks = seq(0, 18, by = 2)) +
  geom_text(aes(x = 2 + 2.6, y = 17, label = "Age at sexual maturity (2)"), color = "black", angle = 0, vjust = 1) +
  theme_minimal()  

marmota_hist




## Histogram of Tursiops truncatus
tursiops <- blood_samples %>% filter(species == "Tursiops truncatus") %>% select(age_years) %>% drop_na()

tursiops_hist <- ggplot(tursiops, aes(x = age_years)) + 
  geom_histogram(
    colour = 1, fill = "white") +
  geom_vline(xintercept = 8.93, color = "navy", linetype = "dashed", size = 0.5) +
  labs(title = "Histogram of Tursiops truncatus age in years (blood samples)", x = "Age in years", y = "Count") +
  scale_x_continuous(breaks = seq(0, 60, by = 4)) +
  geom_text(aes(x = 8.93 + 14, y = 16, label = "Age at sexual maturity (8.93)"), color = "black", angle = 0, vjust = 1) +
  theme_minimal()  

tursiops_hist



## Histogram of Chlorocebus sabaeus
chlorocebus <- blood_samples %>% filter(species == "Chlorocebus sabaeus") %>% select(age_years) %>% drop_na()

chlorocebus_hist <- ggplot(chlorocebus, aes(x = age_years)) + 
  geom_histogram(
    colour = 1, fill = "white") +
  geom_vline(xintercept = 3.92, color = "navy", linetype = "dashed", size = 0.5) +
  labs(title = "Histogram of Chlorocebus sabaeus age in years (blood samples)", x = "Age in years", y = "Count") +
  scale_x_continuous(breaks = seq(0, 28, by = 4)) +
  geom_text(aes(x = 3.92 + 6, y = 20, label = "Age at sexual maturity (3.92)"), color = "black", angle = 0, vjust = 1) +
  theme_minimal()  

chlorocebus_hist

############################# combined plot 2 ####################################
combined_plot_2 <- ggarrange(capreolus_hist, rattus_hist, ovis_hist, 
                             marmota_hist,tursiops_hist, chlorocebus_hist,
                             nrow = 3,
                             ncol = 2)
combined_plot_2
##################################################################################




## Histogram of Felis catus
catus <- blood_samples %>% filter(species == "Felis catus") %>% select(age_years) %>% drop_na()

catus_hist <- ggplot(catus, aes(x = age_years)) + 
  geom_histogram(
    colour = 1, fill = "white") +
  geom_vline(xintercept = 0.79, color = "navy", linetype = "dashed", size = 0.5) +
  labs(title = "Histogram of Chlorocebus sabaeus age in years (blood samples)", x = "Age in years", y = "Count") +
  scale_x_continuous(breaks = seq(0, 28, by = 2)) + scale_y_continuous(breaks = seq(0, 28, by = 4)) +
  geom_text(aes(x = 0.79 + 6, y = 30, label = "Age at sexual maturity (0.79)"), color = "black", angle = 0, vjust = 1) +
  theme_minimal()  

catus_hist




## Histogram of Sus scrofa 
sus <- blood_samples %>% filter(species == "Sus scrofa") %>% select(age_years) %>% drop_na()

sus_hist <- ggplot(sus, aes(x = age_years)) + 
  geom_histogram(
    colour = 1, fill = "white") +
  geom_vline(xintercept = 0.46, color = "navy", linetype = "dashed", size = 0.5) +
  labs(title = "Histogram of Sus scrofa age in years (blood samples)", x = "Age in years", y = "Count") +
  scale_x_continuous(breaks = seq(0, 7, by = 0.5)) + scale_y_continuous(breaks = seq(0, 30, by = 5)) +
  geom_text(aes(x = 0.46 + 2.0, y = 28, label = "Age at sexual maturity (0.46)"), color = "black", angle = 0, vjust = 1) +
  theme_minimal()  

sus_hist


## Histogram of Callithrix jacchus
jacchus <- blood_samples %>% filter(species == "Callithrix jacchus") %>% select(age_years) %>% drop_na()

jacchus_hist <- ggplot(jacchus, aes(x = age_years)) + 
  geom_histogram(
    colour = 1, fill = "white") +
  geom_vline(xintercept = 1.18, color = "navy", linetype = "dashed", size = 0.5) +
  labs(title = "Histogram of Callithrix jacchus age in years (blood samples)", x = "Age in years", y = "Count") +
  scale_x_continuous(breaks = seq(0, 20, by = 1)) + scale_y_continuous(breaks = seq(0, 16, by = 2)) +
  geom_text(aes(x = 1.18 + 3.5, y = 16, label = "Age at sexual maturity (1.18)"), color = "black", angle = 0, vjust = 1) +
  theme_minimal()  

jacchus_hist




## Histogram of Equus quagga
quagga <- blood_samples %>% filter(species == "Equus quagga") %>% select(age_years) %>% drop_na()

quagga_hist <- ggplot(quagga, aes(x = age_years)) + 
  geom_histogram(
    colour = 1, fill = "white") +
  geom_vline(xintercept = 2.47, color = "navy", linetype = "dashed", size = 0.5) +
  labs(title = "Histogram of Equus quagga age in years (blood samples)", x = "Age in years", y = "Count") +
  scale_x_continuous(breaks = seq(0, 24, by = 1)) + scale_y_continuous(breaks = seq(0, 20, by = 4)) +
  geom_text(aes(x = 2.47 + 5.5, y = 20, label = "Age at sexual maturity (2.47)"), color = "black", angle = 0, vjust = 1) +
  theme_minimal()  

quagga_hist



## Histogram of Heterocephalus glaber
heterocephalus <- blood_samples %>% filter(species == "Heterocephalus glaber") %>% select(age_years) %>% drop_na()

heterocephalus_hist <- ggplot(heterocephalus, aes(x = age_years)) + 
  geom_histogram(
    colour = 1, fill = "white") +
  geom_vline(xintercept = 5, color = "navy", linetype = "dashed", size = 0.5) +
  labs(title = "Histogram of Heterocephalus glaber age in years (blood samples)", x = "Age in years", y = "Count") +
  scale_x_continuous(breaks = seq(0, 35, by = 5)) + scale_y_continuous(breaks = seq(0, 24, by = 4)) +
  geom_text(aes(x = 5 + 5, y = 24, label = "Age at sexual maturity (5)"), color = "black", angle = 0, vjust = 1) +
  theme_minimal()  

heterocephalus_hist


## Histogram of Delphinapterus leucas
leucas <- blood_samples %>% filter(species == "Delphinapterus leucas") %>% select(age_years) %>% drop_na()

leucas_hist <- ggplot(leucas, aes(x = age_years)) + 
  geom_histogram(
    colour = 1, fill = "white") +
  geom_vline(xintercept = 11.01, color = "navy", linetype = "dashed", size = 0.5) +
  labs(title = "Histogram of Delphinapterus leucas age in years (blood samples)", x = "Age in years", y = "Count") +
  scale_x_continuous(breaks = seq(0, 60, by = 5)) + scale_y_continuous(breaks = seq(0, 12, by = 4)) +
  geom_text(aes(x = 11.01 + 15, y = 10, label = "Age at sexual maturity (11.01)"), color = "black", angle = 0, vjust = 1) +
  theme_minimal()  

leucas_hist



################################ combined plot 3 #########################################
combined_plot_3 <- ggarrange(catus_hist, sus_hist, jacchus_hist, 
                             quagga_hist, heterocephalus_hist, leucas_hist,
                             nrow = 3,
                             ncol = 2)
combined_plot_3
##################################################################################


## Histogram of Elephas maximus
elephas <- blood_samples %>% filter(species == "Elephas maximus") %>% select(age_years) %>% drop_na()

elephas_hist <- ggplot(elephas, aes(x = age_years)) + 
  geom_histogram(
    colour = 1, fill = "white") +
  geom_vline(xintercept = 9.01, color = "navy", linetype = "dashed", size = 0.5) +
  labs(title = "Histogram of Elephas maximus age in years (blood samples)", x = "Age in years", y = "Count") +
  scale_x_continuous(breaks = seq(0, 100, by = 5)) + scale_y_continuous(breaks = seq(0, 12, by = 4)) +
  geom_text(aes(x = 9.01 + 12, y = 10, label = "Age at sexual maturity 9.01"), color = "black", angle = 0, vjust = 1) +
  theme_minimal()  

elephas_hist

