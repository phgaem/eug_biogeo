#Project: HISTORICAL BIOGEOGRAPHY AND NICHE EVOLUTION IN EUGENIINAE (MYRTACEAE)
#Author: Paulo Henrique Gaem
#University of Michigan, Ann Arbor

# 3. PREPARING GEOGRAPHIC DATASET FOR INFOMAP BIOREGIONS

# Setting working directory
setwd("~/Documents/Projetos/3_PhD/1_Project/Eugenia Biogeography")

#loading packages
library(dplyr)

#loading functions
source('00_functions_eug.R')

# standardising species names in the phylogeny and the geographical occurrence matrix
# to input in Infomap Bioregions

occs <- read.csv(file="data/12_gbif_data_NativeOccs.csv")
occs$PhylAcceptedName <- NA
occs <- occs %>%
  select(gbifID, scientificName, acceptedName, PhylAcceptedName, everything())

reference_table <- read.csv('data/9_reference_table.csv')

occs <- occs %>%
  left_join(reference_table, by = c("acceptedName" = "GBIF.names")) %>%
  mutate(PhylAcceptedName = tip.labels) %>%
  select(-tip.labels, -NMWG.names, -WCVP.manual)

which(!is.na(occs$PhylAcceptedName))
occs <- occs[which(!is.na(occs$PhylAcceptedName)),]

occs$scientificName <- occs$PhylAcceptedName

write.csv(occs, 'data/13_infomap_occs2.csv')

# Run Infomap Bioregions on the website.

# Now, loading the presence-absence matrix given by the program

infomap <- read.delim('results/1_infomap bioregions/13_infomap_occs_presence-absence.txt',
                      header=F, sep="\t", colClasses = "character") %>% as.data.frame

# Replace spaces with underscores in infomap$V1
infomap$V1 <- gsub(" ", "_", infomap$V1)

reference_table <- read.csv('data/9_reference_table.csv')
ref_tips <- reference_table$tip.labels
ref_tips <- data.frame(tip.labels=ref_tips, stringsAsFactors=F)

# Use the match function to find indices where ref_tips$tip.labels matches infomap$V1
match_indices <- match(ref_tips$tip.labels, infomap$V1)

# Assign matched values again in the same way
ref_tips$pres.abs <- infomap$V2[match_indices]

# getting presence-absence table with all terminals saving it to edit manually
write.csv(ref_tips, file='data/14_pre.biogeobears_presabs.csv', row.names=F)

# Now, editing the file to add the presence-absences manually