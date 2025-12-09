#Project: HISTORICAL BIOGEOGRAPHY AND NICHE EVOLUTION IN EUGENIINAE (MYRTACEAE)
#Author: Paulo Henrique Gaem
#University of Michigan, Ann Arbor

# Setting working directory
setwd("~/Documents/Projetos/3_PhD/1_Project/Eugenia Biogeography")

# PLOTTING FIGURE 3

#loading packages
library(dplyr)
library(ape)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(colorspace)

#loading native occurrences together with environmental data
load('data/21_occs_env.RData')

# calculating the mean of each soil variable across depths

occs <- occs %>%
  select(
    gbifID, acceptedName, lon, lat, bio1, bio4, bio12, bio15
  ) %>%
  mutate(
    nitrogen = rowMeans(select(occs, starts_with("nitrogen_")), na.rm = TRUE),
    ocd = rowMeans(select(occs, starts_with("ocd_")), na.rm = TRUE),
    phh2o = rowMeans(select(occs, starts_with("phh2o_")), na.rm = TRUE),
    sand = rowMeans(select(occs, starts_with("sand_")), na.rm = TRUE),
    silt = rowMeans(select(occs, starts_with("silt_")), na.rm = TRUE),
    clay = rowMeans(select(occs, starts_with("clay_")), na.rm = TRUE),
    cfvo = rowMeans(select(occs, starts_with("cfvo_")), na.rm = TRUE),
    bdod = rowMeans(select(occs, starts_with("bdod_")), na.rm = TRUE)
  )

#sand (8) is highly correlated with silt (9) and clay (10). let's remove sand from the database
occs <- occs %>% select(-sand)

# adding tip.labels to the table
ref.table <- read.csv('data/9_reference_table.csv')

occs <- occs %>%
  left_join(ref.table %>% select(tip.labels, GBIF.names), 
            by = c("acceptedName" = "GBIF.names")) %>%
  select(acceptedName, tip.labels, everything())

# creating sectional information

eugtre <- read.tree('data/2_eugtree.tre')

tiplab <- as.data.frame(eugtre$tip.label)
tiplab$sections <- NA

tiplab[1:108,2] <- 'Umbellatae'
tiplab[109:142,2] <- 'Jossinia'
tiplab[143:152,2] <- 'Excelsae'
tiplab[153:159,2] <- 'Speciosae'
tiplab[160:187,2] <- 'Racemosae'
tiplab[188:195,2] <- 'Schizocalomyrtus'
tiplab[196:209,2] <- 'Phyllocalyx'
tiplab[210:231,2] <- 'Pilothecium'
tiplab[232:238,2] <- 'Eugenia'
tiplab[239:252,2] <- 'Pseudeugenia'
tiplab[253:256,2] <- 'Hexachlamys'
tiplab[257:266,2] <- 'Myrcianthes'

tiplab <- tiplab %>%
  rename(tips = `eugtre$tip.label`)

occs <- occs %>%
  left_join(tiplab,
            by=c('tip.labels'='tips')) %>%
  select(acceptedName, tip.labels, sections, everything())


# PCA of environmental variables
envirn <- occs %>% select(-acceptedName, -tip.labels, -sections,  -gbifID,
                          -lon, -lat)

envirn <- envirn %>% #cleaning
  mutate(across(everything(), ~ ifelse(is.na(.) | is.infinite(.), 
                                       mean(., na.rm = TRUE), .)))

pca.res <- prcomp(envirn, center=T, scale.=T)

# variance explained by eigenvectors
sdev <- pca.res$sdev
variance_explained <- sdev^2
percentage_explained <- variance_explained / sum(variance_explained) * 100
explained_df <- data.frame(
  PC = seq_along(percentage_explained),
  Percentage = percentage_explained
)

#loadings: how much of each variable contributes to the principal components

loadings_matrix <- pca.res$rotation
print(loadings_matrix)

squared_loadings <- loadings_matrix^2
print(squared_loadings)

pca_loadings <- list(
  Loadings = as.data.frame(loadings_matrix),
  Squared_Loadings = as.data.frame(squared_loadings)
)

#pca data
pca.dat <- as.data.frame(pca.res$x)

#merging occs.avg with pca data
occs.pca <- cbind(occs[,1:6], pca.dat)

#pca loadings for arrows
loadings <- as.data.frame(pca.res$rotation)
loadings$variable <- rownames(loadings) 

scale_factor <- max(abs(pca.dat$PC1), abs(pca.dat$PC2), abs(pca.dat$PC3)) * 0.5
#scale_factor <- 6

# defining colours for PCA plots

pca.dat$sections <- occs$sections
pca.dat$sections <- factor(pca.dat$sections, levels = c("Myrcianthes", "Pseudeugenia", "Hexachlamys", 'Pilothecium',
                                                        'Eugenia', 'Jossinia', 'Schizocalomyrtus', 'Phyllocalyx',
                                                        'Racemosae', 'Speciosae', 'Excelsae', 'Umbellatae'))

custom_colours <- c(
  "Myrcianthes" = "#4e79a7",
  "Pseudeugenia" = "#f28e2b",
  "Hexachlamys" = "#e15759",
  "Pilothecium" = "#76b7b2",
  "Eugenia" = "#59a14f",
  "Jossinia" = "#edc948",
  "Schizocalomyrtus" = "#AF7AC5",
  "Phyllocalyx" = "#ff9da7",
  "Racemosae" = "#9c755f",
  "Speciosae" = "#bab0ac",
  "Excelsae" = "#d4a6c8",
  "Umbellatae" = "#a8d9a2"
)

# PLOTTING PC1 vs PC2

x_pc <- 'PC1'
y_pc <- 'PC2'

# Create the PCA plot

pca_plot <- ggplot(pca.dat, aes_string(x = x_pc, y = y_pc, colour = "sections")) +
  geom_point() +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.text.y = element_text(size = 10, angle = 90),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 10, face = "italic"),
    legend.title = element_text(size = 10),
    legend.position = "none",
    axis.line = element_line(color = "grey30", size = 0.5),
    axis.line.x.bottom = element_line(color = "grey30", linewidth = 0.2),
    axis.line.y.left = element_line(color = "grey30", linewidth = 0.2)
  ) +
  scale_color_manual(
    values = custom_colours,
    name = "Clades"
  ) +
  scale_y_continuous(
    breaks = c(-2.5, 2.5)
  ) +
  geom_segment(
    data = loadings, aes_string(
      x = "0", y = "0",
      xend = paste0(x_pc, " * scale_factor"),
      yend = paste0(y_pc, " * scale_factor")
    ),
    arrow = arrow(length = unit(0.1, "cm")),
    color = "black"
  ) +
  geom_text_repel(
    data = loadings, aes_string(
      x = paste0(x_pc, " * scale_factor * 1.1"),
      y = paste0(y_pc, " * scale_factor * 1.1"),
      label = "variable"
    ),
    size = 9 / .pt,
    color = "black",
    max.overlaps = Inf  # Optional: show all labels even if crowded
  )

pca_plot

# Save the plot to a file
ggsave(filename = 'figures/pc1-pc2.pdf', plot = pca_plot, device = "pdf", width = 3.39, height = 3,
       units = 'in', dpi = 600, bg = "white")

# PLOTTING PC2 vs PC3

x_pc <- 'PC2'
y_pc <- 'PC3'

# Create the PCA plot

pca_plot2 <- ggplot(pca.dat, aes_string(x = x_pc, y = y_pc, colour = "sections")) +
  geom_point() +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.text.y = element_text(size = 10, angle = 90),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 10, face = "italic"),
    legend.title = element_text(size = 10),
    legend.position = "none",
    axis.line = element_line(color = "grey30", size = 0.5),
    axis.line.x.bottom = element_line(color = "grey30", linewidth = 0.2),
    axis.line.y.left = element_line(color = "grey30", linewidth = 0.2)
  ) +
  scale_color_manual(
    values = custom_colours,
    name = "Clades"
  ) +
  scale_y_continuous(
    breaks = c(-2.5, 2.5)
  ) +
  scale_x_continuous(
    breaks = c(-2.5, 2.5)
  ) +
  geom_segment(
    data = loadings, aes_string(
      x = "0", y = "0",
      xend = paste0(x_pc, " * scale_factor"),
      yend = paste0(y_pc, " * scale_factor")
    ),
    arrow = arrow(length = unit(0.1, "cm")),
    color = "black"
  ) +
  geom_text_repel(
    data = loadings, aes_string(
      x = paste0(x_pc, " * scale_factor * 1.1"),
      y = paste0(y_pc, " * scale_factor * 1.1"),
      label = "variable"
    ),
    size = 9 / .pt,
    color = "black",
    max.overlaps = Inf  # Optional: show all labels even if crowded
  )

pca_plot2

# Save the plot to a file
ggsave(filename = 'figures/pc2-pc3.pdf', plot = pca_plot2, device = "pdf", width = 3.39, height = 3,
  units = 'in', dpi = 600, bg = "white")

# PLOTTING LEGEND

# Plotting full plot with legend

pca_plot2_with_legend <- ggplot(pca.dat, aes_string(x = x_pc, y = y_pc, colour = "sections")) +
  geom_point() +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.text.y = element_text(size = 10, angle = 90),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 10, face = "italic"),
    legend.title = element_text(size = 10),
    legend.position = "right"
  ) +
  scale_color_manual(
    values = custom_colours,
    name = "Clades"
  ) +
  scale_y_continuous(breaks = c(-2.5, 2.5)) +
  scale_x_continuous(breaks = c(-2.5, 2.5)) +
  geom_segment(
    data = loadings, aes_string(
      x = "0", y = "0",
      xend = paste0(x_pc, " * scale_factor"),
      yend = paste0(y_pc, " * scale_factor")
    ),
    arrow = arrow(length = unit(0.1, "cm")),
    color = "black"
  ) +
  geom_text_repel(
    data = loadings, aes_string(
      x = paste0(x_pc, " * scale_factor * 1.1"),
      y = paste0(y_pc, " * scale_factor * 1.1"),
      label = "variable"
    ),
    size = 9 / .pt,
    color = "black",
    max.overlaps = Inf
  )

# extracting legend

legend_only <- cowplot::get_legend(pca_plot2_with_legend)

library(grid)
grid.newpage()
grid.draw(legend_only)

ggsave("figures/legend.pdf",
       cowplot::ggdraw(legend_only),
       width = 1.5, height = 3.5, units = "in", bg = "white")

library(colorspace)

custom_colours2 <- darken(custom_colours, amount = 0.5) # colours of points and contours

# --- VIOLIN PLOT FOR PC1 ---

pc_indices <- 1:3

for (i in pc_indices) {
  
  # Create the ggplot object
  pc_violin_plot <- ggplot(pca.dat, aes(x = sections, y = !!sym(paste0("PC", i)))) + 
    geom_jitter(aes(color = sections), width = 0.2, size = 0.5, alpha = 0.9) +
    geom_violin(aes(fill = sections, color = sections), trim = FALSE, alpha = 0.8, size = 0.3) +
    scale_color_manual(values = custom_colours2, guide = "none") +
    scale_fill_manual(values = custom_colours, name = "Clades") +
    labs(
      x = "Clade (Section)",
      y = paste0("PC", i, " (", round(explained_df$Percentage[i], 1), "%)")
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 10),
      legend.position = "none",
      panel.grid.major.x = element_blank()
    )
  
  # Print the plot to the console
  print(pc_violin_plot)
  
  # Save the plot with a dynamically generated filename
  ggsave(filename = paste0('figures/violin_pc', i, '.pdf'),
         plot = pc_violin_plot, device = "pdf",
         width = 3.22, height = 2.32,
         units = 'in', dpi = 600, bg = "white")
}
