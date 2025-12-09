#Project: HISTORICAL BIOGEOGRAPHY AND NICHE EVOLUTION IN EUGENIINAE (MYRTACEAE)
#Author: Paulo Henrique Gaem
#University of Michigan, Ann Arbor

# Setting working directory
setwd("~/Documents/Projetos/3_PhD/1_Project/Eugenia Biogeography")

# PLOTTING FIGURE 5.

#loading packages
library(ggplot2)
library(patchwork)
library(ggrepel)
library(openxlsx)
library(dplyr)
library(tidyr)
library(colorspace)

#loading data

df <- read.xlsx('results/4_nichevol_ouwieB.xlsx', sheet=1) %>%
  select(-lnL, -AIC, -AICc) %>% separate(clade, into=c('env_var', 'clade'), sep='_')

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

# PC1
p1 <- ggplot(subset(df, env_var == "PC1"),
             aes(x = sig.sq, y = n.shifts, label = clade, color = clade)) +
  geom_errorbar(
    aes(xmin = sig.sq.LB, xmax = sig.sq.UB, y = n.shifts),
    orientation = "y",
    width = 0.2
  ) +
  geom_point() +
  geom_text_repel(color = "black", size = 11/.pt, fontface = "italic", max.overlaps = Inf) +
  scale_color_manual(values = custom_colours) +
  scale_x_continuous(
    breaks = c(-0.0, 0.2, 0.4, 0.6, 0.8),
    limits = c(-0.015, NA)
  ) +
  labs(
    x = expression(sigma^2~"(PC1)"),
    y = "Total shifts/expansions"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    legend.position = "none",
    plot.title = element_blank(),
    axis.title = element_text(size = 11),
    axis.text  = element_text(size = 11)
  ) +
  coord_cartesian(xlim = c(-0.1, NA), ylim = c(-2, NA))

# PC2
p2 <- ggplot(subset(df, env_var == "PC2"),
             aes(x = sig.sq, y = n.shifts, label = clade, color = clade)) +
  geom_errorbar(
    aes(xmin = sig.sq.LB, xmax = sig.sq.UB, y = n.shifts),
    orientation = "y",
    width = 0.2,
  ) +
  geom_point() +
  geom_text_repel(color = "black", size = 11/.pt, fontface = "italic", max.overlaps = Inf) +
  scale_color_manual(values = custom_colours) +
  scale_x_continuous(
    breaks = c(0, 0.03, 0.06, 0.09),
    limits = c(-0.025, NA)
  ) +
  labs(
    x = expression(sigma^2~"(PC2)"),
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    legend.position = "none",
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 11),
    axis.text  = element_text(size = 11)
  ) +
  coord_cartesian(xlim = c(-0.025, NA), ylim = c(-2, NA))

# PC3
p3 <- ggplot(subset(df, env_var == "PC3"),
             aes(x = sig.sq, y = n.shifts, label = clade, color = clade)) +
  geom_errorbar(
    aes(xmin = sig.sq.LB, xmax = sig.sq.UB, y = n.shifts),
    orientation = "y",
    width = 0.2,
  ) +
  geom_point() +
  geom_text_repel(color = "black", size = 11/.pt, fontface = "italic", max.overlaps = Inf) +
  scale_color_manual(values = custom_colours) +
  labs(
    x = expression(sigma^2~"(PC3)"),
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    legend.position = "none",
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 11),
    axis.text  = element_text(size = 11)
  ) +
  coord_cartesian(xlim = c(-0.005, NA), ylim = c(-2, NA))

#PC1 part 2

p4 <- ggplot(subset(df, env_var == "PC1"),
            aes(x = sig.sq, y = n.shf.new, label = clade, color = clade)) +
  geom_errorbar(
    aes(xmin = sig.sq.LB, xmax = sig.sq.UB, y = n.shf.new),
    orientation = "y",
    width = 0.2,
  ) +
  geom_point() +
  geom_text_repel(color = "black", size = 11/.pt, fontface = "italic", max.overlaps = Inf) +
  scale_color_manual(values = custom_colours) +
  scale_x_continuous(
    breaks = c(-0.0, 0.2, 0.4, 0.6, 0.8),
    limits = c(-0.015, NA)
  ) +
  labs(
    x = expression(sigma^2~"(PC1)"),
    y = "Shifts/expansions to new areas"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    legend.position = "none",
    plot.title = element_blank(),
    axis.title = element_text(size = 11),
    axis.text  = element_text(size = 11)
  ) +
  coord_cartesian(xlim = c(-0.1, NA), ylim = c(-2, NA))

# PC2
p5 <- ggplot(subset(df, env_var == "PC2"),
             aes(x = sig.sq, y = n.shf.new, label = clade, color = clade)) +
  geom_errorbar(
    aes(xmin = sig.sq.LB, xmax = sig.sq.UB, y = n.shf.new),
    orientation = "y",
    width = 0.2,
  ) +
  geom_point() +
  geom_text_repel(color = "black", size = 11/.pt, fontface = "italic", max.overlaps = Inf) +
  scale_color_manual(values = custom_colours) +
  scale_x_continuous(
    breaks = c(0, 0.03, 0.06, 0.09),
    limits = c(-0.025, NA)
  ) +
  labs(
    x = expression(sigma^2~"(PC2)"),
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    legend.position = "none",
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 11),
    axis.text  = element_text(size = 11)
  ) +
  coord_cartesian(xlim = c(-0.025, NA), ylim = c(-2, NA))

# PC3
p6 <- ggplot(subset(df, env_var == "PC3"),
             aes(x = sig.sq, y = n.shf.new, label = clade, color = clade)) +
  geom_errorbar(
    aes(xmin = sig.sq.LB, xmax = sig.sq.UB, y = n.shf.new),
    orientation = "y",
    width = 0.2,
  ) +
  geom_point() +
  geom_text_repel(color = "black", size = 11/.pt, fontface = "italic", max.overlaps = Inf) +
  scale_color_manual(values = custom_colours) +
  labs(
    x = expression(sigma^2~"(PC3)"),
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    legend.position = "none",
    plot.title = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 11),
    axis.text  = element_text(size = 11)
  ) +
  coord_cartesian(xlim = c(-0.005, NA), ylim = c(-2, NA))

combined_plot <- (p1 | p2 | p3) / (p4 | p5 | p6) 

ggsave("figures/nichevol_b.pdf", combined_plot, width = 210,   # 85 mm to inches
       height = 170, # 174 mm to inches
       units = "mm", dpi = 600)
