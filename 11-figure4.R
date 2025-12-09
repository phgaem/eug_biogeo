#Project: HISTORICAL BIOGEOGRAPHY AND NICHE EVOLUTION IN EUGENIINAE (MYRTACEAE)
#Author: Paulo Henrique Gaem
#University of Michigan, Ann Arbor

# Setting working directory
setwd("~/Documents/Projetos/3_PhD/1_Project/Eugenia Biogeography")

# PLOTTING FIGURE 4

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(grid)
library(tidyverse)
library(colorspace)
library(cowplot)

# creating colours

palaeo_colour <- "#FF8000"
neotr_colour <- "#8000FF"
placeholder_colour <- NA #"#999999"

palaeo_rgb <- hex2RGB(palaeo_colour)
neotr_rgb <- hex2RGB(neotr_colour)

palaeo_luv <- as(palaeo_rgb, "LUV")
neotr_luv <- as(neotr_rgb, "LUV")

mixed_color_final <- hex(mixcolor(0.5, palaeo_luv, neotr_luv))
mixed_color_final

pantr_colour <- "#C55BA3"

# creating each model data frame

#alpha

pc1_alpha_data <- tibble(
  PC = "PC1 (OUVM)", Parameter = "Alpha", Region = c("Palaeo.", "Neotr."),
  Value = 0.063440385, LB = 0.041561739, UB = 0.109374032
)

pc2_alpha_data <- tibble(
  PC = "PC2 (BMS)", Parameter = "Alpha", Region = c("Palaeo.", "Neotr."),
  Value = 0, LB = 0, UB = 0
)

pc3_alpha_data <- tibble(
  PC = "PC3 (BM1)", Parameter = "Alpha", Region = c("Palaeo.", "Neotr."),
  Value = 0, LB = 0, UB = 0
)

all_alpha_data <- bind_rows(pc1_alpha_data, pc2_alpha_data, pc3_alpha_data) %>%
  mutate(PC = factor(PC, levels = c("PC3 (BM1)", "PC2 (BMS)", "PC1 (OUVM)"))) %>%
  mutate(ColorGroup = case_when(
    Value != 0 ~ "Pantr.",
    TRUE ~ "Placeholder"
  ))

# sigma-square

pc1_sigsq_data <- tibble(
  PC = "PC1 (OUVM)", Parameter = "Sig-sq", Region = c("Palaeo.", "Neotr."),
  Value = c(0.053622612, 0.201685229),
  LB = c(0.03016095, 0.1444505),
  UB = c(0.099348292, 0.302407772)
)

pc2_sigsq_data <- tibble(
  PC = "PC2 (BMS)", Parameter = "Sig-sq", Region = c("Palaeo.", "Neotr."),
  Value = c(0.356494658, 0.042330084),
  LB = c(0.194773853, 0.035533382),
  UB = c(0.530882549, 0.049285966)
)

pc3_sigsq_data <- tibble(
  PC = "PC3 (BM1)", Parameter = "Sig-sq", Region = c("Palaeo.", "Neotr."),
  Value = 0.03380513, LB = 0.028285126, UB = 0.039523609
)

all_sigsq_data <- bind_rows(pc1_sigsq_data, pc2_sigsq_data, pc3_sigsq_data) %>%
  mutate(PC = factor(PC, levels = c("PC3 (BM1)", "PC2 (BMS)", "PC1 (OUVM)"))) %>%
  # NEW: Create ColorGroup based on PC
  mutate(ColorGroup = case_when(
    PC == "PC3 (BM1)" ~ "Pantr.",
    TRUE ~ Region # For PC1 and PC2, use the Region name
  ))

# theta

pc1_optim_data <- tibble(
  PC = "PC1 (OUVM)", Parameter = "Theta", Region = c("Palaeo.", "Neotr."),
  Value = c(-0.638704377, -0.299194771),
  LB = c(-1.852869981, -0.757970517),
  UB = c(0.113035168, 0.053771977)
)

pc2_optim_data <- tibble(
  PC = "PC2 (BMS)", Parameter = "Theta", Region = c("Palaeo.", "Neotr."),
  Value = 0, LB = 0, UB = 0
)

pc3_optim_data <- tibble(
  PC = "PC3 (BM1)", Parameter = "Theta", Region = c("Palaeo.", "Neotr."),
  Value = 0, LB = 0, UB = 0
)

all_optim_data <- bind_rows(pc1_optim_data, pc2_optim_data, pc3_optim_data) %>%
  mutate(PC = factor(PC, levels = c("PC3 (BM1)", "PC2 (BMS)", "PC1 (OUVM)"))) %>%
  mutate(ColorGroup = case_when(
    PC == "PC1 (OUVM)" ~ Region,
    TRUE ~ "Placeholder"
  ))

# creating alpha plot

p_alpha <- ggplot(all_alpha_data, aes(y = PC, x = Value, color = ColorGroup)) +
  
  geom_errorbarh(aes(xmin = LB, xmax = UB), 
                 width = 0.3, size = 0.5, 
                 position = position_identity()) +
  
  geom_point(size = 3, shape = 16, position = position_identity()) +
  
  scale_x_continuous(
    name = expression(alpha),
    breaks = c(0.0, 0.033, 0.067, 0.1)
  ) +
  
  scale_color_manual(
    values = c("Pantr." = pantr_colour, "Placeholder" = placeholder_colour),
    breaks = c("Palaeo.", "Neotr."), 
    labels = c("Palaeo.", "Neotr.")
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    axis.ticks.y = element_line(),
    axis.text.x = element_text(size = 11),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.title.y = element_blank()
  )

p_alpha

# creating sigma plot

p_sigsq <- ggplot(all_sigsq_data, aes(y = PC, x = Value, color = ColorGroup)) +
  
  geom_errorbarh(aes(xmin = LB, xmax = UB), 
                 height = 0.3, size = 0.5, 
                 position = position_identity()) +
  
  geom_point(size = 3, shape = 16, position = position_identity()) +
  
  scale_x_continuous(
    name = expression(sigma^2),
    breaks = scales::pretty_breaks(n = 5)
  ) +
  
  scale_color_manual(
    # Sigma-squared uses all three distinct color groups
    values = c("Palaeo." = palaeo_colour, 
               "Neotr." = neotr_colour, 
               "Pantr." = pantr_colour)
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 11),
    plot.margin = margin(t = 10, r = 5, b = 10, l = 5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.title.y = element_blank()
  )

p_sigsq

# creating theta plot 

p_theta <- ggplot(all_optim_data, aes(y = PC, x = Value, color = ColorGroup)) +
  
  geom_errorbarh(aes(xmin = LB, xmax = UB), 
                 height = 0.3, size = 0.5, 
                 position = position_identity()) +
  
  geom_point(size = 3, shape = 16, position = position_identity()) +
  
  scale_x_continuous(
    name = expression(theta),
    breaks = c(-5, -2.5, 0, 2.5),
    limits = c(-5, 2.5)
  ) +
  
  scale_color_manual(
    values = c("Palaeo." = palaeo_colour, 
               "Neotr." = neotr_colour, 
               "Placeholder" = placeholder_colour)
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 14),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 11),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.title.y = element_blank()
  )

p_theta

# combining plots

combined_plots <- p_alpha | p_sigsq | p_theta

pdf("figures/nichevol_a.pdf", width = 17/2.54, height = 6/2.54)
print(combined_plots)
dev.off()

system("open figures/nichevol_a.pdf")

# Creating legend plot

legend_labels <- c("Palaeotropical", "Neotropical", "Pantropical", "Not calculated")

legend_data <- tibble(
  Category = factor(legend_labels, levels = legend_labels),
  X = 1, Y = 1
)

p_legend <- ggplot(legend_data, aes(x = X, y = Y, color = Category)) +
  geom_point(size = 3, shape = 16) +
  scale_color_manual(
    name = NULL,
    # Map the full labels (Category) to their corresponding colors
    values = c("Palaeotropical" = palaeo_colour,
               "Neotropical" = neotr_colour,
               "Pantropical" = pantr_colour,
               "Not calculated" = placeholder_colour)
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    # Remove plot elements
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p_legend

tiff("figures/nichevol_a_leg.tiff", width = 17, height = 7, units = "cm", res = 600)
print(p_legend)
dev.off()

system("open figures/nichevol_a_leg.tiff")


##############################################################################
##############################################################################

# Add the extracted legend manually below the plots
final_plot <- combined_plots / 
  wrap_elements(legend_guide) +
  plot_layout(heights = c(1, 0.1)) 

# Display the final combined plot
print(final_plot)


combined_plots <- p_alpha | p_sigsq | p_theta

final_plot <- combined_plots / guide_area() +
  plot_layout(guides = "collect", 
              heights = c(1, 0.1), 
              tag_level = 'new')

print(final_plot)

# creating function to make the nine plots

plot_param_panel <- function(df, xlab_expr = NULL, show_yaxis = TRUE, show_xlabel = TRUE,
                             num_xticks = 4, dot_color = "#FF8000") {
  ggplot(df, aes(y = Region, x = Value)) +
    geom_point(size = 3, color = dot_color) +
    geom_errorbarh(aes(xmin = LB, xmax = UB), width = 0.4, size = 0.5, color = dot_color) +
    scale_x_continuous(
      name = if (show_xlabel) xlab_expr else NULL,
      breaks = scales::pretty_breaks(n = num_xticks)
    ) +
    labs(
      title = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.title.x = element_text(size = 12),
      axis.text.y = if (show_yaxis) element_text(size = 12) else element_blank(),
      axis.ticks.y = if (show_yaxis) element_line() else element_blank(),
      axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# creating the nine plots 

p1 <- plot_param_panel(pc1_alpha, xlab_expr = expression(alpha), show_yaxis = T, show_xlabel = F)
p2 <- plot_param_panel(pc1_sigsq, xlab_expr = expression(sigma^2), show_yaxis = F, show_xlabel = F)
p3 <- plot_param_panel(pc1_optim, xlab_expr = expression(theta), show_yaxis = F, show_xlabel = F)

p4 <- plot_param_panel(zero_param, xlab_expr = expression(alpha), show_yaxis = T, show_xlabel = F,
                       dot_color = "#999999")
p5 <- plot_param_panel(pc2_sigsq, xlab_expr = expression(sigma^2), show_yaxis = F, show_xlabel = F)
p6 <- plot_param_panel(zero_param, xlab_expr = expression(theta), show_yaxis = F, show_xlabel = F,
                       dot_color = "#999999")

p7 <- plot_param_panel(zero_param, xlab_expr = expression(alpha), show_yaxis = T, show_xlabel = T,
                       dot_color = "#999999")
p8 <- plot_param_panel(pc3_sigsq, xlab_expr = expression(sigma^2), show_yaxis = F, show_xlabel = T)
p9 <- plot_param_panel(zero_param, xlab_expr = expression(theta), show_yaxis = F, show_xlabel = T,
                       dot_color = "#999999")

pc1_label <- wrap_elements(full = grid::textGrob("    PC1 (OUVM)", x = 0, hjust = 0,
                                                 gp = gpar(fontsize = 12, fontface = "bold")))
pc2_label <- wrap_elements(full = grid::textGrob("   PC2 (BMS)", x = 0, hjust = 0,
                                                 gp = gpar(fontsize = 12, fontface = "bold")))
pc3_label <- wrap_elements(full = grid::textGrob("    PC3 (BM1)", x = 0, hjust = 0,
                                                 gp = gpar(fontsize = 12, fontface = "bold")))

final_plot <- pc1_label / (p1 + p2 + p3) /
  pc2_label / (p4 + p5 + p6) /
  pc3_label / (p7 + p8 + p9) +
  plot_layout(heights = c(0.06, 1, 0.06, 1, 0.06, 1))

tiff("figures/nichevol_a.tiff", width = 17, height = 14, units = "cm", res = 600)
print(final_plot)
dev.off()

system("open figures/nichevol_a.tiff")
