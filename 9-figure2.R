#Project: HISTORICAL BIOGEOGRAPHY AND NICHE EVOLUTION IN EUGENIINAE (MYRTACEAE)
#Author: Paulo Henrique Gaem
#University of Michigan, Ann Arbor

# Setting working directory
setwd("~/Documents/Projetos/3_PhD/1_Project/Eugenia Biogeography")

# PLOTTING FIGURE 2

#loading packages
library(sf)
library(ggplot2)
library(ggspatial)

# FIGURE 2. 

#reading shapefiles

land <- st_read("figures/shapefiles/ne_10m_land.shp")
a.southam <- st_read("figures/shapefiles/a_southamerica.shp")
b.centralam <- st_read("figures/shapefiles/b_centralamerica.shp")
c.caribbean <- st_read("figures/shapefiles/c_caribbean.shp")
d.africa <- st_read("figures/shapefiles/d_africa.shp")
e.madag <- st_read("figures/shapefiles/e_madagascar.shp")
f.reumau <- st_read("figures/shapefiles/f_reumau.shp")
g.southasia <- st_read("figures/shapefiles/g_southasia.shp")
h.malay <- st_read("figures/shapefiles/h_malayetc.shp")
i.newcal <- st_read("figures/shapefiles/i_newcaledonia.shp")
j.pacific <- st_read("figures/shapefiles/j_pacific.shp")

my_crs <- 8857
#my_crs <- "+proj=natearth +lon_0=50"

#transforming projection (Equal Earth, EPSG:8857)

land_eq <- st_transform(land, crs = my_crs)
a.southam <- st_transform(a.southam, crs = my_crs)
b.centralam <- st_transform(b.centralam, crs = my_crs)
c.caribbean <- st_transform(c.caribbean, crs = my_crs)
d.africa <- st_transform(d.africa, crs = my_crs)
e.madag <- st_transform(e.madag, crs = my_crs)
f.reumau <- st_transform(f.reumau, crs = my_crs)
g.southasia <- st_transform(g.southasia, crs = my_crs)
h.malay <- st_transform(h.malay, crs = my_crs)
i.newcal <- st_transform(i.newcal, crs = my_crs)
j.pacific <- st_transform(j.pacific, crs = my_crs)

# Define points as sf objects in longitude/latitude (WGS84)
indon_point <- st_sfc(st_point(c(108.5, 6.8)), crs = 4326)
sa_point     <- st_sfc(st_point(c(-57, -12)), crs = 4326)
congo_point  <- st_sfc(st_point(c(23, -1)), crs = 4326)

# Transform to Equal Earth projection (EPSG:8857)
indon_point_eq <- st_transform(indon_point, crs = my_crs)
sa_point_eq     <- st_transform(sa_point, crs = my_crs)
congo_point_eq  <- st_transform(congo_point, crs = my_crs)

# Extract coordinates
indon_coords <- st_coordinates(indon_point_eq)
sa_coords     <- st_coordinates(sa_point_eq)
congo_coords  <- st_coordinates(congo_point_eq)

#plotting map

ggplot() +
  geom_sf(data = land_eq, fill = "white", color = "grey40", size = 0.1) +
  geom_sf(data = a.southam, fill = "#e31a1c", color = "grey40", size = 0.1) +
  geom_sf(data = b.centralam, fill = "#E3EF3E", color = "grey40", size = 0.1) +
  geom_sf(data = c.caribbean, fill = "gold", color = "grey40", size = 0.1) +
  geom_sf(data = d.africa, fill = "#4daf4a", color = "grey40", size = 0.1) +
  geom_sf(data = e.madag, fill = "darkolivegreen3", color = "grey40", size = 0.1) +
  geom_sf(data = f.reumau, fill = "#90FCE6", color = "grey40", size = 0.1) +
  geom_sf(data = g.southasia, fill = "#1f78b4", color = "grey40", size = 0.1) +
  geom_sf(data = h.malay, fill = "#C766FF", color = "grey40", size = 0.1) +
  geom_sf(data = i.newcal, fill = "#FF66D1", color = "grey40", size = 0.1) +
  geom_sf(data = j.pacific, fill = "#8466FF", color = "grey40", size = 0.1) +
  annotate("text", x = indon_coords[1],  y = indon_coords[2],  label = "103", size = 7, fontface = "bold", color = "grey20") +
  annotate("text", x = sa_coords[1],      y = sa_coords[2],      label = "1018",   size = 7, fontface = "bold", color = "grey20") +
  annotate("text", x = congo_coords[1],   y = congo_coords[2],   label = "157",   size = 7, fontface = "bold", color = "grey20") +
  coord_sf(expand = F) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),  # Remove x axis title (including "x")
    axis.title.y = element_blank(),  # Remove y axis title (including "y")
    panel.grid = element_blank()
  )

ggsave(
  filename = "figures/map_disp_eugeniinae.tiff",
  width = 255,      # width in mm
  height = 140,    # height in mm (adjust as needed for your map's best appearance)
  units = "mm",
  dpi = 600,
  device = "tiff"   # or use "tiff" for some journals, check their submission guidelines
)

# INSERT - Mascarene Islands

ggplot() +
  geom_sf(data = land_eq, fill = "white", color = "grey40", size = 0.1) +
  geom_sf(data = f.reumau, fill = "#90FCE6", color = "grey40", size = 0.1) +
  coord_sf(xlim = c(5000000, 5500000),
           ylim = c(-2750000, -2450000),
           expand = T)+
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    panel.grid = element_blank()
  )

ggsave(
  filename = "figures/map_disp_mascarene.tiff",
  width = 255/8,      # width in mm
  height = 140/8,    # height in mm (adjust as needed for your map's best appearance)
  units = "mm",
  dpi = 600,
  device = "tiff"   # or use "tiff" for some journals, check their submission guidelines
)

# INSERT - New Caledonia

ggplot() +
  geom_sf(data = land_eq, fill = "white", color = "grey40", size = 0.1) +
  geom_sf(data = i.newcal, fill = "#FF66D1", color = "grey40", size = 0.1) +
  coord_sf(xlim = c(15000000, 15900000),
           ylim = c(-2990000, -2400000),
           expand = T)+
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "aliceblue"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    panel.grid = element_blank()
  )

ggsave(
  filename = "figures/map_disp_newcal.tiff",
  width = 255/8,      # width in mm
  height = 140/8,    # height in mm (adjust as needed for your map's best appearance)
  units = "mm",
  dpi = 600,
  device = "tiff"   # or use "tiff" for some journals, check their submission guidelines
)
