#Project: HISTORICAL BIOGEOGRAPHY AND NICHE EVOLUTION IN EUGENIINAE (MYRTACEAE)
#Author: Paulo Henrique Gaem
#University of Michigan, Ann Arbor

# 5. EXTRACTING NODES AGES FROM PHYLOGENY

# Setting working directory
setwd("~/Documents/Projetos/3_PhD/1_Project/Eugenia Biogeography")

#loading packages
library(ape)
library(treeio)
library(tidytree)
library(dplyr)
library(tidyr)
library(writexl)

# reading Eugeniinae tree used for biogeographical reconstruction
eugtre <- read.tree('data/2_eugtree.tre')

eug.labels <- eugtre$tip.label # getting tip.labels

# some tips are a little bit different between the eugeniinae tree made with NMWG's Appendix 7 and 
#their tree with nodes ages and confidence intervals (Appendix 5). We need to fix that before matching tips

to_replace <- c( # eugeniinae tree tips
  "Eugenia.talbotii_Parmentier.5037", "Eugenia.cocosensis_VS.7857",
  "Eugenia.retinadenia_Oviedo.sn", "Eugenia.verticillata_Duarte.sn_ESA85678",
  "Eugenia.calcadensis_Erkings.sn", "Eugenia.capensis.subsp.zeyheri_Maurin.1800",
  "Eugenia.varia_Oviedo.sn", "Eugenia.membranifolia_Duarte.sn_ESA85677",
  "Eugenia.pacifica_VS.7887", "Eugenia.hypargyrea_CarillaReyes.5231",
  "Eugenia.walkerae_Walker.95.016", "Eugenia.pseudopsidium_Santiago.2017.7",
  "Eugenia.capensis.subsp.natalitia_Maurin.1796", "Eugenia.sp_Faria.3140",
  "Eugenia.macrobracteolata_Faria.3051"
)

replacements <- c( # wording used in NMWG's Appendix 5
  "Eugenia.talbotii_PM.5037", "Eugenia.cocosensis_IgeaVS.7857",
  "Eugenia.retinadenia_OviedoChaves.sn", "Eugenia.verticillata_Duarte.ESA85678",
  "Eugenia.calcadensis_DharmarErkings.sn", "Eugenia.zeyheri_Maurin.1800",
  "Eugenia.varia_OviedoChaves.sn", "Eugenia.membranifolia_Duarte.ESA85677",
  "Eugenia.pacifica_IgeaVS.7887", "Eugenia.hypargyrea_CarilloReyes.5231",
  "Eugenia.walkerae_Walker.95016", "Eugenia.pseudopsidium_Santiago.20177",
  "Eugenia.natalitia_Maurin.1796", "Eugenia.luschnathiana_Faria.3140",
  "Eugenia.macrobracteolata_Faria.3050"
)

# fixing tips
eug.labels.fixed <- eug.labels

for (i in seq_along(to_replace)) {
  idx <- which(eug.labels.fixed == to_replace[i])
  if (length(idx) > 0) eug.labels.fixed[idx] <- replacements[i]
}

# reading NMWG (2024) phylogenetic 
myrtre2 <- read.beast('data/20_mmc_target_common_HPD.tre')

# defining function to clean myrtre2 tip.labels (i.e. drops all characters before the capital letter)
clean_label <- function(x) sub("^[^A-Z]*", "", x)

myrtre2@phylo$tip.label <- clean_label(myrtre2@phylo$tip.label) # cleaning myrtre2's tip labels
tips_to_drop <- setdiff(myrtre2@phylo$tip.label, eug.labels.fixed) # finding tips to drop
pruned_myrtre2 <- drop.tip(myrtre2, tips_to_drop) #dropping tips in myrtre2

length(which(eug.labels.fixed %in% pruned_myrtre2@phylo$tip.label)) # checking. all set.

# now, eugtre used in biogeobears corresponds to the pruned Appendix 5 tree.
# let's plot this tree in order to identify nodes from which we want to extract information.

pdf(file = 'figures/pruned_myrtre2_with_nodes.pdf', width = 8, height = 35) #create a pdf

# plot tree with node numbers
plot(
  pruned_myrtre2@phylo,
  main = "Eugeniinae tree with nodes labels",
  cex = 0.5,
  type = "phylogram"
)

# Calculate node numbers for labelling
num_tips <- Ntip(pruned_myrtre2@phylo)
num_nodes <- Nnode(pruned_myrtre2@phylo)

# Check the number of tips just to be safe (it should be 266)
message("The tree has ", num_tips, " tips.")

# Node labels start from N_tips + 1
node_labels <- (num_tips + 1):(num_tips + num_nodes)

# Add node labels (identifiers)
nodelabels(
  text = node_labels,
  cex = 0.5, 
  bg = "yellow",
  frame = "rect"
)

dev.off()

system("open figures/pruned_myrtre2_with_nodes.pdf")

# Extracting informations from nodes

nodes <- c(
  267, 523, 531, 268, 520, 507, 500, 479, 465, 458, 430, 424,
  415, 382, 274, 501, 378, 416, 439, 316, 324, 294, 312, 382, 383,
  412, 401, 384, 396, 397, 300
)

nodes <- unique(nodes)

tree_tbl <- as_tibble(pruned_myrtre2)

node_info <- tree_tbl %>%
  filter(node %in% nodes)

node_info <- node_info %>%
  select(node, CAheight_0.95_HPD, CAheight_mean)

node_info <- node_info %>%
  unnest_wider(CAheight_0.95_HPD, names_sep = "_") %>%
  rename(
    CAheight_0.95_HPD_lower = CAheight_0.95_HPD_1,
    CAheight_0.95_HPD_upper = CAheight_0.95_HPD_2
  )

node_info <- node_info %>%
  select(node, CAheight_mean, CAheight_0.95_HPD_lower, CAheight_0.95_HPD_upper)

write_xlsx(node_info, "results/2_biogeobears/node_info.xlsx")