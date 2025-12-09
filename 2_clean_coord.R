#Project: HISTORICAL BIOGEOGRAPHY AND NICHE EVOLUTION IN EUGENIINAE (MYRTACEAE)
#Author: Paulo Henrique Gaem
#University of Michigan, Ann Arbor

# 2. CLEANING GEOGRAPHIC COORDINATES

# Setting working directory
setwd("~/Documents/Projetos/3_PhD/1_Project/Eugenia Biogeography")

#loading packages
library(dplyr)
library(readr)
library(ggplot2)
library(countrycode)
library(CoordinateCleaner)
library(sf)
library(rWCVP)

#loading functions
source('00_functions_eug.R')

#loading raw GBIF data
dat <- read_tsv('data/4_gbif_data.csv')

#Selecting columns of interest
dat <- dat %>%
  dplyr::select(gbifID, scientificName, decimalLongitude, decimalLatitude,
                countryCode, taxonRank, basisOfRecord)

#plot data raw to get an overview
wm <- borders("world", colour = "gray50", fill = "gray50")
ggplot() +
  coord_fixed() +
  wm +
  geom_point(data = dat,
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "darkred",
             size = 0.5) +
  theme_bw()

#Converting country code from ISO2c to ISO3c
dat$countryCode <-  countrycode(dat$countryCode, 
                                origin =  'iso2c',
                                destination = 'iso3c')

#Detecting "red flags" (possible problems)
dat <- data.frame(dat)
flags <- clean_coordinates(x = dat, 
                           lon = "decimalLongitude", 
                           lat = "decimalLatitude",
                           countries = "countryCode",
                           species = "scientificName",
                           tests = c("capitals", "centroids", "equal",
                                     "gbif", "institutions","zeros"))

summary(flags)

#Plotting flags to visualise
plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")

#Excluding some problematic records
dat_cl <- subset(flags, flags$.equ=='TRUE'); dat_cl$.equ=NULL #Identical coordinates
dat_cl <- subset(dat_cl, dat_cl$.zer=='TRUE'); dat_cl$.zer=NULL #Coordinates equal to zero
dat_cl$.val=NULL

sum(dat_cl$.gbf=="FALSE") #no occurrences at the headquarters of GBIF
dat_cl$.gbf=NULL

#Removing coordinates of taxa on genus level (we only want to use species and below)
dat_cl <- subset(dat_cl, !dat_cl$taxonRank=='GENUS'); dat_cl$taxonRank=NULL

#Exporting all occurrence points with CoordCl flags tests results
write.csv(dat_cl, file="data/5_geo_flags_CoordCl.csv", row.names=F)

#Excluding coordinates in capitals of countries and territories, except those of small insular ones
for (i in 1:length(excl_cap)){
  dat_cl <- subset(dat_cl, !.cap=='FALSE' | !countryCode==excl_cap[i])
}
dat_cl$.cap=NULL

#Excluding coordinates in centroids of countries, except those of small insular ones
for (i in 1:length(excl_cen)){
  dat_cl <- subset(dat_cl, !.cen=='FALSE' | !countryCode==excl_cen[i])
}
dat_cl$.cen=NULL

#Excluding coordinates in botanical institutions
dat_cl <- subset(dat_cl, !.inst=='FALSE')
dat_cl$.inst=NULL

dat_cl$.summary=NULL #excluding summary flags column

#Exporting dataset cleaned with CoordinateCleaner
write.csv(dat_cl, file="data/6_geodata_cleaned_CoordCl.csv")

load('data/3_spp_names_across.RData')

#Getting accepted scientific names for all occurrences

dat_cl$acceptedName <- NA #creating the accepted name column
dat_cl <- dat_cl %>% select(1:2, acceptedName, everything())  #Rearranging columns

for(i in 1:length(spp_namesGBIF)) {
  one_species <- spp_namesGBIF[[i]]
  dat_cl$acceptedName[which(dat_cl$scientificName %in% one_species)] <- one_species[1]
  cat(i, "\r")
}

#Some species are considered synonyms on GBIF, but they are treated separately on
#the NMWG paper. We need to change the accepted names in dat_cl for those separately

#A vector with names to be written on the column acceptedNames
accNames <- c("Eugenia pisonis O.Berg", "Eugenia fajardensis (Krug & Urb.) Urb.")

#A vector with names accepted by GBIF, but to be changed to names of 'accNames'
rejNames <- c("Eugenia moschata (Aubl.) Nied. ex T.Durand & B.D.Jacks.",
              "Myrcianthes fragrans (Sw) Mc Vaugh")

for (aa in 1:length(accNames)){
  indexNames <- as.vector(which(dat_cl$scientificName==accNames[aa]))
  for (i in 1:length(indexNames)){
    dat_cl[indexNames[i], "acceptedName"] <- accNames[aa]
  }
}

#All occurrences of Myrcianthes cruciata are named as Myrcianthes pseudomato on GBIF. Lets change this.

#Getting the native range of Myrcianthes cruciata from POWO
native_range <- wcvp_distribution("Myrcianthes cruciata", taxon_rank="species",
                                  introduced=FALSE, extinct=FALSE, 
                                  location_doubtful=FALSE)

#Visualising the native distribution of Myrcianthes cruciata
(p <- wcvp_distribution_map(native_range, crop_map=TRUE) + 
    theme(legend.position="none"))

# Creating a subset of my spatial data containing points named "Myrcianthes pseudomato" by GBIF
occs_myr <- dat_cl %>% 
  select(acceptedName, decimalLatitude, decimalLongitude) %>%
  st_as_sf(coords=c("decimalLongitude", "decimalLatitude"), crs = st_crs(4326)) %>%
  subset(acceptedName == "Myrcianthes pseudomato (D.Legrand) Mc Vaugh")

#creating a column on the subset spatial data frame containing the information of whether
#the points fall within the native range of Myrcianthes cruciata or not

occs_myr$native <- st_intersects(occs_myr,
                                 st_union(native_range), #merging all native polygons
                                 sparse=FALSE)[,1]

#retrieving row names of all occurrences "TRUE" for native range of Myrcianthes cruciata
Mcru <- subset(occs_myr, occs_myr$native=="TRUE") %>% row.names()

#now changing the acceptedName of these occurrences on the dat_cl, the geographical dataset with accepted names of species
dat_cl[Mcru, "acceptedName"] <- "Myrcianthes cruciata M.Ibrahim & Proença"

#Ok. Now all Brazilian specimens previously known as Myrcianthes pseudomato have the correct name Myrcianthes cruciata.

#Saving the geographical coordinate table with correct acceptedNames
write.csv(dat_cl, file="data/7_gbif_data_accName.csv")

#Putting GBIF names in a vector
acc_namesGBIF <- rep(NA, 266)

for(i in 1:length(acc_namesGBIF)){
  one_species <- spp_namesGBIF[[i]]
  acc_namesGBIF[i] <- one_species[1]
}
# Creating table containing tip labels, NMWG names and GBIF names reference names
reference_table <- data.frame(tip.labels=tiplabb, NMWG.names=spp_namesNMWG,
                              GBIF.names=acc_namesGBIF)

# Adding POWO (WCVP) names to reference_table
reference_table$WCVP.manual <- NA
WCVP_manual <- read.csv('data/8_WCVP_manual.csv')
reference_table$WCVP.manual <- WCVP_manual
colnames(reference_table[,4]) <- "WCVP.manual"

write.csv(reference_table, "data/9_reference_table.csv", row.names=F)

#Cleaning coordinates to species' native range with rWCVP (adapted from https://matildabrown.github.io/rWCVP/articles/coordinate-cleaning.html)

#Selecting species name, latitude and longitude columns and transforming into a spatial data frame
occs <- 
  dat_cl %>% 
  select(acceptedName, decimalLatitude, decimalLongitude) %>%
  st_as_sf(coords=c("decimalLongitude", "decimalLatitude"), crs = st_crs(4326))

#Creating a subset of the reference table in which only rows with taxa present at POWO and GBIF are present
ref_table2 <- subset(reference_table, !reference_table$WCVP.manual=='NONE')

for (i in 1:nrow(ref_table2)){
  # Getting the native range of species
  native_range <- wcvp_distribution(ref_table2[i,'WCVP.manual'], taxon_rank="species",
                                    introduced=FALSE, extinct=FALSE, 
                                    location_doubtful=FALSE)
  
  # Plotting the native range of the species
  p <- wcvp_distribution_map(native_range, crop_map=TRUE) + 
    theme(legend.position="none")
  
  # Creating a subset of my spatial data containing only data from the species in analysis
  occs2 <- subset(occs, occs$acceptedName == ref_table2[i,'GBIF.names'])
  
  # Creating a buffer around the native range of the species
  buffered_dist <- native_range %>%
    st_union() %>%
    st_buffer(0.54) #60 km buffer (near the equator)
  
  # Creating a column on the spatial data frame containing the information of whether
  # the points fall within the buffer around the native range of the species or not
  occs2$native_buffer <- st_intersects(occs2, buffered_dist, sparse=FALSE)[,1]
  
  # Determine if there is only one category in native_buffer
  n_distinct_native_buffer <- n_distinct(factor(occs2$native_buffer))
  
  # Create the base plot with the buffer
  p <- p + 
    geom_sf(data=buffered_dist, fill="transparent", col="gold") +
    labs(title = paste(i, '. ', ref_table2[i, 'NMWG.names'], sep=''))
  
  # Conditionally add the points layer
  if (n_distinct_native_buffer == 1) {
    p <- p + geom_sf(data=occs2, 
                     fill=c("gold")[factor(occs2$native_buffer)],
                     col="black", 
                     shape=21, size=1.3)
  } else {
    p <- p + geom_sf(data=occs2, 
                     fill=c("red","gold")[factor(occs2$native_buffer)],
                     col="black", 
                     shape=21, size=1.3)
  }
  
  # Save to PDF
  filename <- paste("figures/distr_plots/", i, "_", ref_table2[i,'NMWG.names'], ".pdf", sep="")
  pdf(filename, 14, 9)
  print(p)
  dev.off()
  cat(i, "\r")
}

#Expect plots visually

#All plots have been visually expected and it seems that everything is OK. Now, we need to delete the points that are
#too far from the native range of each species (i.e., outside the native range and its buffer)

dat_cl$Native <- NA #creating the native/exotic column
dat_cl <- dat_cl %>% select(1:6, Native, everything())  #Rearranging columns

#Getting each species' occurrence points statuses as either native or exotic
for (i in 1:nrow(ref_table2)){
  
  # Getting the native range of species
  native_range <- wcvp_distribution(ref_table2[i,'WCVP.manual'], taxon_rank="species",
                                    introduced=FALSE, extinct=FALSE, 
                                    location_doubtful=FALSE)
  
  # Creating a subset of my spatial data containing only data from the species in analysis
  occs2 <- subset(occs, occs$acceptedName == ref_table2[i,'GBIF.names'])
  
  if(nrow(occs2)==0){ cat(i, "\r") } else {
    
    # Creating a buffer around the native range of the species
    buffered_dist <- native_range %>%
      st_union() %>%
      st_buffer(0.54) #60 km buffer (near the equator)
    
    # Creating a column on the spatial data frame containing the information of whether
    # the points fall within the buffer around the native range of the species or not
    occs2$native_buffer <- st_intersects(occs2, buffered_dist, sparse=FALSE)[,1]
    
    # Checking which row names of occs2 are present in dat_cl
    common_rows <- rownames(occs2)[rownames(occs2) %in% rownames(dat_cl)]
    
    # Using these common row names to update the Native column in dat_cl
    dat_cl[common_rows, "Native"] <- occs2[common_rows, "native_buffer"]
    cat(i, "\r")
  }
}

dat_cl3 <- dat_cl[is.na(dat_cl$Native),] #Retrieving rows for which neither True nor False was placed in dat_cl$Native

#These species do not belong to our study. They have been retrieved with our GBIF dataset because a species affinis
# (e.g., Myrcia aff. epithet) was present in the list of names to be harmonised with the GBIF database.

#Now, filtering from dat_cl only the native occurrences
dat_clNative <- subset(dat_cl, dat_cl$Native=="TRUE")

occsNative <- dat_clNative %>% #transforming in a spatial data frame
  select(acceptedName, decimalLatitude, decimalLongitude) %>%
  st_as_sf(coords=c("decimalLongitude", "decimalLatitude"), crs = st_crs(4326))

#Inspecting native occurrences visually

for (i in 1:nrow(ref_table2)){
  # Getting the native range of species
  native_range <- wcvp_distribution(ref_table2[i,'WCVP.manual'], taxon_rank="species",
                                    introduced=FALSE, extinct=FALSE, 
                                    location_doubtful=FALSE)
  
  # Plotting the native range of the species
  p <- wcvp_distribution_map(native_range, crop_map=TRUE) + 
    theme(legend.position="none")
  
  # Creating a subset of my spatial data containing only data from the species in analysis
  occs2 <- subset(occsNative, occsNative$acceptedName == ref_table2[i,'GBIF.names'])
  
  # Creating a buffer around the native range of the species
  buffered_dist <- native_range %>%
    st_union() %>%
    st_buffer(0.54) #60 km buffer (near the equator)
  
  # Create the base plot with the buffer
  p <- p + 
    geom_sf(data=buffered_dist, fill="transparent", col="gold") +
    labs(title = paste(i, '. ', ref_table2[i, 'NMWG.names'], sep=''))
  
  p <- p + geom_sf(data=occs2, 
                   fill=c("gold"),
                   col="black", 
                   shape=21, size=1.3)
  
  # Save to PDF
  filename <- paste("figures/distr_native/", i, "_", ref_table2[i,'NMWG.names'], ".pdf", sep="")
  pdf(filename, 14, 9)
  print(p)
  dev.off()
  cat(i, "\r")
}

#Now, coordinates of species not present on GBIF have to be added manually to dat_clNative.
#These refer to taxa with uncertain names (typically with one geographical point per
#species, e.g. Eugenia sp3).

#reading manual occurrences and adding them to the dat_clNative object
manual_occur <- read.csv('data/10_manual_occurrences.csv')
dat_clNativeComplete <- rbind(dat_clNative, manual_occur)

#All good. Now, saving geographic datasets again:
write.csv(dat_cl, file="data/11_gbif_data_NativeStatusSpp.csv", row.names=F) #all occurrences with Native/Exotic statuses
write.csv(dat_clNativeComplete, file="data/12_gbif_data_NativeOccs.csv", row.names=F) #only native occurrences