#Project: HISTORICAL BIOGEOGRAPHY AND NICHE EVOLUTION IN EUGENIINAE (MYRTACEAE)
#Author: Paulo Henrique Gaem
#University of Michigan, Ann Arbor

# 1. MANIPULATING PHYLOGENY AND GETTING GBIF DATA

# Setting working directory
setwd("~/Documents/Projetos/3_PhD/1_Project/Eugenia Biogeography")

#loading packages
library(ape)
library(rgbif)

#loading functions
source('00_functions_eug.R')

# 1. Reading Neotropical Myrteae tree of NMWG (2024) and keeping only Eugeniinae 

#reading and plotting myrteae phylogenetic tree with one terminal per species (NMWG 2024 Appendix S7)
myrtre <- read.tree('data/1_mmc_target_common_Oct17_pruned.tre')

#retrieving Eugeniinae tip labels
eug.labels <- myrtre$tip.label[172:440]
eug.labels <- as.vector(eug.labels$x)

#pruning phylogeny in order to only keep Eugeniinae
eugtre <- keep.tip(myrtre, eug.labels)

#removing duplicate species
tip.drop <- c('Eugenia.monticola_Flickinger.2015DR27', 'Eugenia.stirpiflora_VIIS_1658',
              'Eugenia.uniflora_Genome.KR867678_ICN193277')

eugtre <- drop.tip(eugtre, tip.drop)

#Saving filtered tree as tree file
write.tree(eugtre, file='data/2_eugtree.tre')

#Extracting and saving tip labels
tiplabb <- eugtre$tip.label

#Extracting taxon names from tip labels
spp_namesNMWG <- sapply(strsplit(tiplabb, "_"), function(x) {
  species <- unlist(strsplit(x[1], "\\."))
  species <- gsub("\\.var\\.", " var. ", species)
  paste(species, collapse = " ")
})

# Print result
print(spp_namesNMWG)

#Changing a wrong name
spp_namesNMWG[spp_namesNMWG == "Myrcianthes cruciatus"] <- "Myrcianthes cruciata"

#Applying the function to all species names
spp_namesGBIF <- resolveGBIF_synonyms(spp_namesNMWG)

save(tiplabb, spp_namesNMWG, spp_namesGBIF, file='data/3_spp_names_across.RData')

#Preparing the list of names to be searched on GBIF
resolved_names <- unlist(spp_namesGBIF)

#Searching and requesting occurrence data from GBIF using both accepted names and synonyms
email <- "email"
user <- "user"
pwd <- "pwd"

occ_result <- rgbif::occ_download(rgbif::pred_in("scientificName", resolved_names),
                                  pred_in("basisOfRecord", c('PRESERVED_SPECIMEN','MATERIAL_CITATION')),
                                  pred("hasCoordinate", TRUE),
                                  format = "SIMPLE_CSV", 
                                  user=user,pwd=pwd,email=email) # Sending request to GBIF

c <- occ_download_get('0005129-250415084134356') %>%
  occ_download_import()
