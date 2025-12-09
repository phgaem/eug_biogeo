#Project: HISTORICAL BIOGEOGRAPHY AND NICHE EVOLUTION IN EUGENIINAE (MYRTACEAE)
#Author: Paulo Henrique Gaem
#University of Michigan, Ann Arbor

# 6. GETTING ENVIRONMENTAL DATA BASED ON OCCURRENCE POINTS

# Setting working directory
setwd("~/Documents/Projetos/3_PhD/1_Project/Eugenia Biogeography")

#loading packages
library(ape)
library(terra)
library(dplyr)
library(ggplot2)
library(plotly)
library(openxlsx)

#loading native occurrences
occs <- read.csv("data/12_gbif_data_NativeOccs.csv")
occs <- occs %>% select(gbifID, acceptedName, lon=decimalLongitude, lat=decimalLatitude)
occs <- subset(occs, !is.na(occs$lat))
length(unique(occs$acceptedName)) #246 (out of 266) Eugeniinae species with occurrence points

# extracting CHELSA BioClim values from data points

chelsa.vars <- c("bio1", "bio4", "bio12", "bio15")

for (i in 1:length(chelsa.vars)){
  
  r <- rast(paste0('data/0_EnvRasters/CHELSA/CHELSA_', chelsa.vars[i], '_1981-2010_V.2.1.tif')) #loading raster
  
  extr <- terra::extract(r, occs[, 3:4], search_radius=10000) #extracting values
  
  occs[[paste(chelsa.vars[i])]] <- NA #creating new column
  
  occs[,ncol(occs)] <- extr[,2] #paste values in new column
  
  cat("Completed processing:", chelsa.vars[i], "\n") #progress
}

#defining reference raster (chelsa bioclim) whose CRS will be used to transform the soilgrids rasters 
b <- rast('data/0_EnvRasters/CHELSA/CHELSA_bio1_1981-2010_V.2.1.tif')
  
# extracting SoilGrid values from data points

attributes <- c("bdod", "cfvo", "clay", "sand", "silt", "nitrogen", "ocd", "phh2o")
depths <- c("0-5", "5-15", "15-30", "30-60", "60-100", "100-200")

for (attribute in attributes) {
  for (depth in depths) {
    # Construct file name based on attribute and depth
    file_name <- paste0('data/0_EnvRasters/SoilGrids/', attribute, '_', depth, 'cm_mean_5000.tif')
    
    # Load raster
    r <- rast(file_name)
    
    #transform raster
    r_transformed <- project(r, crs(b), method = "bilinear")
    
    # Extract values
    extr <- terra::extract(r_transformed, occs[, 3:4], search_radius = 50000)
    
    # Create new column name
    column_name <- paste(attribute, depth, sep = "_")
    
    # Add new column to 'occs' and paste extracted values
    occs[[column_name]] <- NA
    occs[, ncol(occs)] <- extr[, 2]
    
    cat("Completed processing:", attribute, depth, "\n") # Progress
  }
}

save(occs, file='data/21_occs_env.RData')

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

# test for correlation between variables using Pearson's index and removing one
# in each highly correlated (> 0.7) pair

num_vars <- select(occs, -gbifID, -acceptedName, -lon, -lat)

# Compute the correlation matrix
cor_matrix <- cor(num_vars, use = "pairwise.complete.obs", method = "pearson")

# Identify highly correlated pairs (absolute correlation > 0.7), same threshold for negative correlations
corr.pairs <- which(abs(cor_matrix) > 0.7, arr.ind = T)
corr.pairs <- corr.pairs[which(!corr.pairs[,1]==corr.pairs[,2]),]

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

write.csv(explained_df, 'data/22_PCA_eigenExplain2.csv', row.names=F)

#loadings: how much of each variable contributes to the principal components

loadings_matrix <- pca.res$rotation
print(loadings_matrix)

squared_loadings <- loadings_matrix^2
print(squared_loadings)

pca_loadings <- list(
  Loadings = as.data.frame(loadings_matrix),
  Squared_Loadings = as.data.frame(squared_loadings)
)

write.xlsx(pca_loadings, "data/23_PCA_loadings.xlsx", rowNames=T)

#pca data
pca.dat <- as.data.frame(pca.res$x)

#merging occs.avg with pca data
occs.pca <- cbind(occs[,1:6], pca.dat)

#saving
write.csv(occs, file='data/24_occs_envRaw.csv', row.names=F)
write.csv(occs.pca, file='data/25_occs_envPCA.csv', row.names=F)

#pca loadings for arrows
loadings <- as.data.frame(pca.res$rotation)
loadings$variable <- rownames(loadings) 

scale_factor <- max(abs(pca.dat$PC1), abs(pca.dat$PC2), abs(pca.dat$PC3)) * 0.5
#scale_factor <- 6

# defining colours for PCA plots

pca.dat$sections <- occs$sections
  
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

# Define the PC pairs you want to plot
pc_pairs <- list(c("PC1", "PC2"), c("PC1", "PC3"), c("PC2", "PC3"))

# 1. plotting PCA in two dimensions for individuals

# Iterate over each pair
for (pair in pc_pairs) {
  # Extract the x and y PCs
  x_pc <- pair[1]
  y_pc <- pair[2]
  
  # Create a dynamic plot title based on the PCs
  plot_title <- paste("PCA of Chelsa/Bioclim + SoilGrids Variables:", x_pc, "vs", y_pc)
  
  # Create the PCA plot
  pca_plot <- ggplot(pca.dat, aes_string(x = x_pc, y = y_pc, colour = "sections")) +
    geom_point() +
    theme_minimal() +
    labs(title = plot_title, x = x_pc, y = y_pc) +
    scale_color_manual(values = custom_colours) +
    # Plot arrows for loadings
    geom_segment(
      data = loadings, aes_string(
        x = "0", y = "0",
        xend = paste0(x_pc, " * scale_factor"),
        yend = paste0(y_pc, " * scale_factor")
      ),
      arrow = arrow(length = unit(0.2, "cm")),
      color = "black"
    ) +
    # Annotate with variable names
    geom_text(
      data = loadings, aes_string(
        x = paste0(x_pc, " * scale_factor * 1.1"),
        y = paste0(y_pc, " * scale_factor * 1.1"),
        label = "variable"
      ),
      size = 3,
      color = "black",
      hjust = 0.5,
      vjust = 0.5
    )
  
  # Create a dynamic file name based on the PCs
  filename <- paste0("figures/PCA_plots/1_EnvPCA_", x_pc, "-", y_pc, ".tif")
  
  # Save the plot to a file
  ggsave(
    filename = filename,
    plot = pca_plot, device = "tiff", width = 8, height = 6,
    dpi = 300, bg = "white"
  )
}

# 2. plotting PCA in two dimensions for species means with species labels

occs.pca1 <- occs.pca %>% #creating the data.frame with species names and eigenvalue averages
  group_by(tip.labels) %>%
  summarise(
    sections = unique(sections),
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE),
    PC3 = mean(PC3, na.rm = TRUE),
  )

occs.pca1$short_labels <- gsub("_.*", "", occs.pca1$tip.labels) #shortening species names to plot
occs.pca1$short_labels <- sub("Eugenia\\.", "E\\.", occs.pca1$short_labels)
occs.pca1$short_labels <- sub("Myrcianthes\\.", "M\\.", occs.pca1$short_labels)

# Get unique sections to loop through
unique_sections <- unique(occs.pca1$sections)

# Loop through each section
for (current_section in unique_sections) {
  
  # Filter the data for the current section
  occs.pca1_subset <- occs.pca1[occs.pca1$sections == current_section, ]
  
  # Loop through each pair of PCs (inner loop)
  for (pair in pc_pairs) {
    # Extract the x and y PCs
    x_pc <- pair[1]
    y_pc <- pair[2]
    
    # Create a dynamic plot title based on the PCs and the current section
    plot_title <- paste("PCA of Chelsa/Bioclim + SoilGrids (", current_section, "): ", x_pc, " vs ", y_pc, sep = "")
    
    # Create the PCA plot
    # Crucially, we now plot using occs.pca1_subset for points and labels
    pca_plot <- ggplot(occs.pca1_subset, aes_string(x = x_pc, y = y_pc)) +
      geom_point(aes_string(colour = "sections")) +
      geom_text(aes(label = short_labels, colour = sections),
                hjust = -0.1, vjust = 0.5, size = 2.5, show.legend=F) +
      theme_minimal() +
      labs(title = plot_title, x = x_pc, y = y_pc, colour='Section') +
      scale_color_manual(values = custom_colours) +
      
      # Plot arrows for loadings (these remain constant across sections)
      geom_segment(
        data = loadings, aes_string(
          x = "0", y = "0",
          xend = paste0(x_pc, " * scale_factor"),
          yend = paste0(y_pc, " * scale_factor")
        ),
        arrow = arrow(length = unit(0.2, "cm")),
        color = "black"
      ) +
      
      # Annotate with variable names (these also remain constant)
      geom_text(
        data = loadings, aes_string(
          x = paste0(x_pc, " * scale_factor * 1.1"),
          y = paste0(y_pc, " * scale_factor * 1.1"),
          label = "variable"
        ),
        size = 3,
        color = "black",
        hjust = 0.5,
        vjust = 0.5
      )
    
    # Create a dynamic file name based on the PCs and the current section
    # This ensures unique filenames for each plot
    filename <- paste0("figures/PCA_plots/EnvPCAspp/", current_section, "_EnvPCAspp_", x_pc, "-", y_pc, ".tif")
    
    # Save the plot to a file
    ggsave(
      filename = filename,
      plot = pca_plot, device = "tiff", width = 8, height = 6,
      dpi = 300, bg = "white"
    )
  }
}






# Iterate over each pair
for (pair in pc_pairs) {
  # Extract the x and y PCs
  x_pc <- pair[1]
  y_pc <- pair[2]
  
  # Create a dynamic plot title based on the PCs
  plot_title <- paste("PCA of Chelsa/Bioclim + SoilGrids Variables:", x_pc, "vs", y_pc)
  
  # Create the PCA plot
  pca_plot <- ggplot(occs.pca1, aes_string(x = x_pc, y = y_pc)) +
    geom_point(aes_string(colour = "sections")) +
    geom_text(aes(label = short_labels, colour = sections),
              hjust = -0.1, vjust = 0.5, size = 2.5) +
    theme_minimal() +
    labs(title = plot_title, x = x_pc, y = y_pc) +
    scale_color_manual(values = custom_colours) +
    # Plot arrows for loadings
    geom_segment(
      data = loadings, aes_string(
        x = "0", y = "0",
        xend = paste0(x_pc, " * scale_factor"),
        yend = paste0(y_pc, " * scale_factor")
      ),
      arrow = arrow(length = unit(0.2, "cm")),
      color = "black"
    ) +
    # Annotate with variable names
    geom_text(
      data = loadings, aes_string(
        x = paste0(x_pc, " * scale_factor * 1.1"),
        y = paste0(y_pc, " * scale_factor * 1.1"),
        label = "variable"
      ),
      size = 3,
      color = "black",
      hjust = 0.5,
      vjust = 0.5
    )
  
  # Create a dynamic file name based on the PCs
  filename <- paste0("figures/PCA_plots/2_EnvPCAspp_", x_pc, "-", y_pc, ".tif")
  
  # Save the plot to a file
  ggsave(
    filename = filename,
    plot = pca_plot, device = "tiff", width = 8, height = 6,
    dpi = 300, bg = "white"
  )
}






# 3D scatter plot with PC1, PC2, and PC3

p <- plot_ly(
  data = pca.dat,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~sections,
  colors = custom_colours,
  type = 'scatter3d',
  mode = 'markers',
  marker = list(size = 4, opacity = 0.5),
  name = ~sections
) %>%
  layout(
    title = "3D PCA Plot: PC1 vs PC2 vs PC3",
    scene = list(
      xaxis = list(title = 'PC1'),
      yaxis = list(title = 'PC2'),
      zaxis = list(title = 'PC3')
    )
  )

# Add the loading arrows to the plot
for (i in 1:nrow(loadings)) {
  p <- p %>% add_trace(
    x = c(0, loadings$PC1[i] * scale_factor),
    y = c(0, loadings$PC2[i] * scale_factor),
    z = c(0, loadings$PC3[i] * scale_factor),
    type = 'scatter3d',
    mode = 'lines',
    line = list(color = 'black', width = 3),
    showlegend = FALSE, # Don't show individual arrows in legend
    inherit = FALSE, # Important to prevent plotly from trying to map 'sections'
    name = loadings$variable[i] # Name for hover info on the line
  ) %>%
    add_trace(
      x = loadings$PC1[i] * scale_factor,
      y = loadings$PC2[i] * scale_factor,
      z = loadings$PC3[i] * scale_factor,
      type = 'scatter3d',
      mode = 'text',
      text = loadings$variable[i],
      textfont = list(color = 'black', size = 10),
      showlegend = FALSE, # Don't show individual labels in legend
      inherit = FALSE,
      name = loadings$variable[i] # Name for hover info on the text
    )
}

# Print the plot to view it
p

# To save the 3D interactive plot as an HTML file:
htmlwidgets::saveWidget(as_widget(p), "figures/PCA_plots/2_3D_EnvPCA_Biplot.html")

library(colorspace)

custom_colours2 <- darken(custom_colours, amount = 0.3)
custom_colours3 <- darken(custom_colours, amount = 0.7)

# --- VIOLIN PLOT FOR PC1 ---
pc1_violin_plot <- ggplot(pca.dat, aes(x = sections, y = PC1)) +
  
  # 1. Add Jittered Points: Uses custom_colours2 (darker shade)
  geom_jitter(
    aes(color = sections), 
    width = 0.2,
    size = 0.5,
    alpha = 0.8
  ) +
  
  # 2. Add Violin Plot: Fill uses custom_colours, Contour (color) uses custom_colours2
  geom_violin(
    aes(fill = sections, color = sections), # <--- CHANGE: Mapped 'color' to 'sections'
    trim = FALSE, 
    alpha = 0.7 # Keep original alpha
  ) +
  
  # 3. Apply custom_colours2 to the 'color' aesthetic (used by BOTH points and contour)
  scale_color_manual(
    values = custom_colours2, # <--- Both points and contour will use this darker color
    guide = "none" 
  ) +
  
  # 4. Apply original custom_colours to the 'fill' aesthetic (for violins)
  scale_fill_manual(
    values = custom_colours,
    name = "Clades"
  ) +
  labs(
    title = "Distribution of Clades along PC1",
    x = "Clade (Section)",
    y = paste0("PC1 (", round(explained_df$Percentage[1], 1), "%)")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "italic"),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  ) +
  scale_x_discrete(limits = levels(pca.dat$sections))
# End of PC1 plot code

# Display the plot
print(pc1_violin_plot)
