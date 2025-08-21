


setwd("/Users/aakinad/Downloads")

install.packages("rotl")
library(rotl)

# Original species vector
blood_species <- c("Acinonyx jubatus", "Aonyx cinereus", "Apodemus sylvaticus", "Bos taurus",
                   "Callithrix jacchus", "Canis lupus familiaris", "Capreolus capreolus", "Ceratotherium simum simum",
                   "Cervus canadensis", "Cervus elaphus", "Chlorocebus sabaeus", "Delphinapterus leucas",
                   "Elephas maximus", "Equus caballus", "Equus quagga", "Felis catus",
                   "Heterocephalus glaber", "Homo sapiens", "Lagenorhynchus obliquidens", "Loxodonta africana",
                   "Macaca mulatta", "Macropus giganteus", "Marmota flaviventris", "Mus musculus",
                   "Odobenus rosmarus divergens", "Orcinus orca", "Osphranter rufus", "Ovis aries",
                   "Phoca vitulina", "Rattus norvegicus", "Sus scrofa", "Sus scrofa domesticus",
                   "Tursiops truncatus", "Zalophus californianus")


resolved_blood_species <- tnrs_match_names(blood_species)
resolved_blood_species
blood_species_tree <- tol_induced_subtree(ott_ids = ott_id(resolved_blood_species))

plot(blood_species_tree, cex = 0.7, col = "navy")

### Getting the tree format
library(ape)
write.tree(blood_species_tree, file = "blood_species_tree.newick")
cat(write.tree(blood_species_tree))

class(blood_species_tree)

### Read the tree
tree <- read.tree("/Users/aakinad/Downloads/blood_species_tree.newick")
plot (tree)
tree$tip.label



###### SKIN samples species
skin_species <- c("Antrozous pallidus", "Bos taurus", "Carollia perspicillata", "Cephalorhynchus hectori hectori",
                  "Cynopterus brachyotis", "Delphinapterus leucas", "Desmodus rotundus", "Eidolon helvum", "Eptesicus fuscus", 
                  "Equus quagga", "Heterocephalus glaber", "Homo sapiens","Lagenorhynchus obliquidens", "Leptonycteris yerbabuenae",
                  "Macaca mulatta", "Megaptera novaeangliae", "Molossus molossus", "Mus musculus", "Myotis lucifugus", "Myotis myotis",     
                  "Myotis vivesi", "Orcinus orca", "Oreamnos americanus", "Phoca vitulina", "Phyllostomus discolor", "Phyllostomus hastatus", 
                  "Pteropus hypomelanus", "Pteropus poliocephalus", "Pteropus pumilus", "Pteropus rodricensis", "Pteropus vampyrus", "Rattus norvegicus",
                  "Rhinolophus ferrumequinum", "Rhynchonycteris naso", "Rousettus aegyptiacus", "Saccopteryx bilineata", "Tadarida brasiliensis",
                  "Tursiops aduncus", "Tursiops truncatus", "Zalophus californianus")

library(rotl)
#### building the skin species tree
resolved_skin_species <- tnrs_match_names(skin_species)
resolved_skin_species
skin_species_tree <- tol_induced_subtree(ott_ids = ott_id(resolved_skin_species))

plot(skin_species_tree, cex = 0.7, col = "navy")

### Getting the tree format
library(ape)
write.tree(skin_species_tree, file = "skin_species_tree.newick")

View(skin_species_tree)
skin_species_tree$tip.label

