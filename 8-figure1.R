#Project: HISTORICAL BIOGEOGRAPHY AND NICHE EVOLUTION IN EUGENIINAE (MYRTACEAE)
#Author: Paulo Henrique Gaem
#University of Michigan, Ann Arbor

# Setting working directory
setwd("~/Documents/Projetos/3_PhD/1_Project/Eugenia Biogeography")

# PLOTTING FIGURE 1

#loading packages
library(ape)
library(aplot)
library(ggtree)
library(ggtreeExtra)
library(ggforce)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(BioGeoBEARS)

# FIGURE 1. Phylogenetic hypothesis, biogeographical reconstruction and variation in environmental niches across species

eugtre <- read.tree('data/2_eugtree.tre') # loading phylogenetic tree

#loading continuous trait table (from PCA summarisation)
envirn <- read.csv('data/24_occs_envPCA.csv') %>%
  group_by(tip.labels) %>%
  summarise(
    PC1.mean = mean(PC1, na.rm = TRUE),
    PC1.sd = sd(PC1, na.rm = TRUE),
    PC2.mean = mean(PC2, na.rm = TRUE),
    PC2.sd = sd(PC2, na.rm = TRUE),
    PC3.mean = mean(PC3, na.rm = TRUE),
    PC3.sd = sd(PC3, na.rm = TRUE)
  )

trait_aligned <- data.frame(tip.label = eugtre$tip.label, stringsAsFactors = FALSE) %>%
  left_join(envirn, by = c("tip.label" = "tip.labels"))

#loading biogeobears object corresponding to the the lowest AICc model
load('results/2_biogeobears/5_Eugenia_DECJ-East.Rdata')

# Getting the biogeographical areas
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=np(res$inputs$geogfn))
areas = getareas_from_tipranges_object(tipranges)
state_indices_0based = areas_list_to_states_list_new(
  areas=areas,
  maxareas=res$inputs$max_range_size,
  include_null_range=TRUE
)
state_indices_0based

#extracting probabilities on each node from 'res'
probs <- as.data.frame(res$ML_marginal_prob_each_state_at_branch_top_AT_node)

#assigning range names to probs' column names
names(probs) <- sapply(state_indices_0based, paste, collapse = "")

probs <- as.data.frame(probs[, -1]) # drop null range ("_")
probs$node <- 1:nrow(probs) # Add a 'node' column with the correct node numbers

probs_long <- probs %>%
  pivot_longer(-node, names_to = "state", values_to = "prob")
probs_long$node <- as.integer(probs_long$node)

# Base colours for single areas and null range
base_colours <- c(
  "_" = "#D3D3D3", "A" = "#e31a1c", "B" = "#E3EF3E", "C" = "gold",
  "D" = "#4daf4a", "E" = "darkolivegreen3", "F" = "#90FCE6", "G" = "#1f78b4",
  "H" = "#C766FF", "I" = "#FF66D1", "J" = "#8466FF"
)

# Blending helper (for ranges consisting in more than one area)
blend_colours <- function(hex_colors) {
  rgb_matrix <- t(col2rgb(hex_colors))
  blended_rgb <- apply(rgb_matrix, 2, mean)
  rgb(blended_rgb[1], blended_rgb[2], blended_rgb[3], maxColorValue = 255)
}

all_ranges <- colnames(probs)
my_colours <- base_colours
for (range in all_ranges) {
  if (nchar(range) <= 1) next
  individual_areas <- strsplit(range, "")[[1]]
  colours_to_blend <- base_colours[individual_areas]
  my_colours[range] <- blend_colours(colours_to_blend)
}

scale_fill_manual(values = my_colours)

# Preparing biogeographical pies
p_ggtree <- ggtree(eugtre)

tree_dat <- p_ggtree$data
visual_tip_order <- tree_dat$label[tree_dat$isTip][order(tree_dat$y[tree_dat$isTip], decreasing = F)]

pies_plotdata <- probs_long %>%
  left_join(tree_dat[,c("node","x","y")], by = "node") %>%
  group_by(node) %>%
  arrange(state, .by_group = TRUE) %>%
  mutate(
    prob_sum = sum(prob),
    prob_frac = prob / prob_sum,
    start = 2 * pi * cumsum(lag(prob_frac, default = 0)),
    end = 2 * pi * cumsum(prob_frac)
  ) %>%
  ungroup()

pies_plotdata <- pies_plotdata %>%
  mutate(
    pie_radius = ifelse(node <= Ntip(eugtre), 0.6, 1)
  )

# getting ancestral areas present in any pies in the ancestral reconstruction occupying at least 30% of the chart
# (for the legend of the biogeographical analysis)
state_maxprob <- pies_plotdata %>%
  group_by(state) %>%
  summarise(max_prob = max(prob, na.rm = TRUE))

signif_states <- state_maxprob$state[state_maxprob$max_prob >= 0.3]

my_colours_signif <- my_colours[signif_states]
my_colours_signif

# building base pie tree 

pie_plot <- ggtree(eugtre, size = 0.04, colour = 'grey48') +
#  geom_tiplab(size=0.5) + hexpand(.6) +
  geom_arc_bar(
    data = pies_plotdata,
    aes(x0 = x, y0 = y, r0 = 0, r = pie_radius, start = start, end = end, fill = state),
    inherit.aes = FALSE,
    color = NA,
    linewidth = 0.06,
    alpha = 1
  ) +
  scale_fill_manual(values = my_colours) +
  coord_fixed() +
  theme_tree2() +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8, color = "grey48"),
        axis.line.x = element_line(color = "grey48", linewidth = 0.2),
        axis.ticks.x = element_line(color = "grey48", linewidth = 0.2),
        plot.margin = margin(5, 50, 5, 15)
        )

# Labelling clades with sections names

section_ranges <- list( #tip ranges according to visual_tip_order
  Myrc = c(1, 10),
  H = c(11, 14),
  Pseud = c(15, 28),
  Eu = c(29, 35),
  Pilothecium = c(36, 57),
  Sch = c(58, 65),
  Phylloc = c(66, 79),
  Racemosae = c(80, 107),
  Sp = c(108, 114),
  Exc = c(115, 124),
  Jossinia = c(125, 158),
  Umbellatae = c(159, 266)
)

sections_df <- do.call(rbind, lapply(names(section_ranges), function(sec) {
  idx <- section_ranges[[sec]]
  data.frame(
    sections = sec,
    tiplab1 = visual_tip_order[idx[1]],
    tiplab2 = visual_tip_order[idx[2]],
    stringsAsFactors = FALSE
  )
}))

# Convert the tree to a data frame to get tip coordinates
tree_data <- ggtree::fortify(eugtre)
max_x <- max(tree_data$x)

# Join your sections_df with tree data to get the 'y' coordinates
sections_df_coords <- sections_df %>%
  dplyr::left_join(tree_data, by = c("tiplab1" = "label")) %>%
  dplyr::rename(y_start = y) %>%
  dplyr::left_join(tree_data, by = c("tiplab2" = "label")) %>%
  dplyr::rename(y_end = y) %>%
  dplyr::select(sections, y_start, y_end)

# Calculate the midpoint of each clade for the label position
sections_df_coords <- sections_df_coords %>%
  dplyr::rowwise() %>%
  dplyr::mutate(y_mid = mean(c(y_start, y_end)))

# adding bars and labels to pie_plot
pie_plot_sect <- pie_plot +
  geom_rect(
    data = sections_df_coords,
    aes(
      xmin = max_x + 0.9,
      xmax = max_x + 5.1,
      ymin = y_start,
      ymax = y_end
    ),
    fill = "#3A56A4",
    color = NA,
    linewidth = 0.2,
    inherit.aes = F
  ) +
  geom_text(
    data = sections_df_coords,
    aes(
      x = max_x + 3,
      y = y_mid,
      label = sections
    ),
    hjust = 0.5,
    angle = -90,
    size = 8 / ggplot2::.pt,
    color = "white",
    fontface = 'italic',
    inherit.aes = F
  ) +
  coord_fixed(ratio = 1) +
  hexpand(0.3)


# plotting environmental niche values for all species (PC1-3)
fruit_long <- trait_aligned %>%
  pivot_longer(cols = c(PC1.mean, PC2.mean, PC3.mean, PC1.sd, PC2.sd, PC3.sd),
               names_to = c("PC", ".value"),
               names_pattern = "(PC\\d)\\.(mean|sd)")

fruit_long$PC <- factor(fruit_long$PC, levels = c("PC1", "PC2", "PC3"))
fruit_long$tip.label <- factor(fruit_long$tip.label, levels = visual_tip_order)

n_tips <- length(visual_tip_order)

rect_df <- fruit_long %>%
  group_by(PC) %>%
  summarise(
    min_x = min(mean - sd, na.rm = TRUE) - 0.2,
    max_x = max(mean + sd, na.rm = TRUE) + 0.2,
    .groups = "drop"
  ) %>%
  mutate(
    ymin = 0,
    ymax = n_tips + 1
  )

fruit_panel <- ggplot(fruit_long) +
  geom_rect(
    data = rect_df,
    aes(xmin = min_x, xmax = max_x, ymin = ymin, ymax = ymax, group = PC),
    fill = NA,
    color = "grey88",
    linewidth = 0.2,
    inherit.aes = FALSE
  ) +
  geom_linerange(aes(y = tip.label, xmin = mean - sd, xmax = mean + sd), linewidth = 0.2, color = 'grey48', show.legend = FALSE) +
  geom_point(aes(y = tip.label, x = mean), fill = 'grey48', shape = 21, color = 'grey48', stroke = 0.12, size = 0.8, show.legend = FALSE) +
  scale_color_manual(values = c(PC1 = 'grey48', PC2 = 'grey48', PC3 = 'grey48')) +
  scale_fill_manual(values = c(PC1 = 'grey48', PC2 = 'grey48', PC3 = 'grey48')) +
  facet_wrap(vars(PC), ncol = 3, strip.position = "bottom", scales = 'free_x') +
  scale_x_continuous(
    breaks = sort(c(seq(0, 6, by = 2.5), seq(-2.5, -6, by = -2.5))),
    labels = function(x) format(x, nsmall = 1)
  ) +
  scale_y_discrete(limits = visual_tip_order) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 8, color = "grey48", angle = 90, vjust = 0.5, hjust = 1),
    axis.line.x = element_line(color = "grey88", linewidth = 0.2),
    axis.ticks.x = element_line(color = "grey88", linewidth = 0.2),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.placement = "outside",
    strip.text = element_blank()
  ) +
  labs(x = NULL, y = NULL)

#creating the final plot
final_plot <- fruit_panel %>% insert_left(pie_plot_sect, width = 0.999999999)

fig_width  <- 168 / 25.4
fig_height <- 230 / 25.4

ggsave("figures/biogeo_traitDisp_photos.pdf", final_plot,
       width = fig_width, height = fig_height, units = "in", dpi = 600)

system("figures/biogeo_traitDisp_photos.pdf")

# making the legend
legend_df <- data.frame(state = names(my_colours_signif))

# Create the standalone legend plot
legend_plot <- ggplot(legend_df, aes(x = state, y = 1, fill = state)) +
  # Use geom_tile to create a visible, but tiny, plot element
  geom_tile(width = 0.01, height = 0.01, color = NA) +
  
  # Set the colors and name for the legend
  scale_fill_manual(
    values = my_colours_signif,
    name = NULL
  ) +
  
  # Use guides() to customize the legend keys
  guides(fill = guide_legend(
    override.aes = list(
      shape = 15, # Use a solid square for the legend key
      size = 3 # Make the legend key larger
    )
  )) +
  
  # Theme to hide all plot elements and leave only the legend
  theme_void() +
  theme(
    legend.position = "left",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8, colour = "grey48"),
    legend.key.size = unit(0.25, "cm")
  )

legend_plot

ggsave("figures/legend_biogeo.pdf", legend_plot,
       width = fig_width, height = fig_height, units = "in", dpi = 600)