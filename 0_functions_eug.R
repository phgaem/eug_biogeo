#Project: HISTORICAL BIOGEOGRAPHY AND NICHE EVOLUTION IN EUGENIINAE (MYRTACEAE)
#Author: Paulo Henrique Gaem
#University of Michigan, Ann Arbor

# 0. FUNCTIONS AND VECTORS USED IN THE PIPELINE

#Retrieving all accepted and synonym names of each species the phylogeny
resolveGBIF_synonyms <- function(name) {
  gnr_resolve_x <- function(x) {
    taxon_info <- name_backbone(name = x)
    if("scientificName" %in% colnames(taxon_info)) {
      all_names <- taxon_info$scientificName
      synonyms <- name_usage(key = taxon_info$usageKey, data = "synonyms")$data
      if(nrow(synonyms)>0) {
        synonym_names <- synonyms$scientificName
        all_names <- c(all_names, synonym_names)
      }
      if(is.null(all_names)) {
        all_names <- paste0("UNMATCHED_",x)
      }  
    } else {
      all_names <- paste0("UNMATCHED_",x)
    }
    return(all_names)
  }
  new.names <- list(unname(pbapply::pbsapply(name, gnr_resolve_x)))
  return(new.names[[1]])
}

#Countries and territories that should be included in step
#'excluding coordinates in capital cities'

#####
excl_cap <- c("FRA", #France
              "BLZ", #Belize
              "GBR", #United Kingdom
              "DZA", #Algeria
              "BRA", #Brazil
              "ZWE", #Zimbabwe
              "GUY", #Guyana
              "CRI", #Costa Rica
              "JPN", #Japan
              "SLV", #El Salvador
              "VEN", #Venezuela
              "NIC", #Nicaragua
              "PAN", #Panama
              "HTI", #Haiti
              "SUR", #Suriname
              "HND", #Honduras
              "GTM", #Guatemala
              "CUB", #Cuba
              "ECU", #Ecuador
              "DOM", #Dominican Republic
              "COL", #Colombia
              "CMR", #Cameroon
              "AUS", #Australia
              "THA", #Thailand
              "GHA", #Ghana
              "BEL", #Belgium
              "BEN", #Benin
              "BOL", #Bolivia
              "PER", #Peru
              "PRY", #Paraguay
              "BDI", #Burundi
              "MDG", #Madagascar
              "COD", #Congo DR
              "ESP", #Spain
              "SWZ", #Eswatini
              "LUX", #Luxembourg
              "URY", #Uruguay
              "SLE", #Sierra Leone
              "TGO", #Togo
              "MOZ", #Mozambique
              "USA", #United States
              "GIN", #Guinea
              "MLI", #Mali
              "TZA", #Tanzania
              "TCD", #Chad
              "LBN", #Lebanon
              "TWN", #Taiwan
              "NZL", #New Zealand
              "CHN", #China
              "PAK", #Pakistan
              "KEN", #Kenya
              "NPL", #Nepal
              "GEO") #Georgia
#####

#Countries and territories that should be included in step
#'excluding coordinates in centroids'

#####
excl_cen <- c("PRY", #Paraguay
              "BRA", #Brazil
              "VEN", #Venezuela
              "PER", #Peru
              "BOL", #Bolivia
              "BLZ", #Belize
              "GUY", #Guyana
              "MEX", #Mexico
              "DOM", #Dominican Republic
              "IND", #India
              "TUN", #Tunisia
              "PAN", #Panama
              "GUF", #French Guiana
              "SUR", #Suriname
              "COL", #Colombia
              "URY", #Uruguay
              "CHL", #Chile
              "ESP", #Spain
              "CRI", #Costa Rica
              "SLV", #El Salvador
              "BDI", #Burundi
              "GIN", #Guinea
              "ARG", #Argentina
              "FRA", #France
              "CHN", #China
              "NZL", #New Zealand
              "CUB") #Cuba