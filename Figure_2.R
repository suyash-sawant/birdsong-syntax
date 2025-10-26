# Load libraries
library(sf)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(ggpattern)
library(shadowtext)
library(ggspatial)
library(cowplot)  # For legend extraction

# Load KML
kml_path <- "KODAI WBS AmNat final.kml"
layers <- st_layers(kml_path)$name
territory_layers <- layers[str_detect(layers, "20")]
study_area <- st_read(kml_path, layer = "Study Area", quiet = TRUE)

territories <- bind_rows(lapply(territory_layers, function(layer) {
  st_read(kml_path, layer = layer, quiet = TRUE) %>%
    mutate(Year = str_extract(layer, "20\\d{2}"))
}))
study_area <- st_transform(study_area, st_crs(territories))

territories_cropped <- territories %>%
  mutate(
    Name = str_trim(Name),
    Individual = str_remove(Name, "_20\\d{2}"),
    Description = as.numeric(as.character(Description))
  ) %>%
  st_intersection(study_area)

# Set colors by Description
category_colors <- c(
  "1" = "#1f78b4",  # Shared – light purple
  "0" = "#a6cee3",  # Full – rich medium purple
  "2" = "#d9d9d9"   # Unidentified – light gray (unchanged)
)
category_labels <- c("1" = "Banded Individuals - Songs Analyzed",
                     "0" = "Banded Individuals - No Songs Analyzed", 
                     "2" = "Unanded Individuals")

# Plot function
plot_year <- function(df, year_label, show_legend = FALSE) {
  df_year <- filter(df, Year == year_label)
  df_labels <- df_year %>%
    filter(Description %in% c(0, 1)) %>%
    mutate(pt = suppressWarnings(st_point_on_surface(geometry)),
           lon = st_coordinates(pt)[,1],
           lat = st_coordinates(pt)[,2])
  
  ggplot() +
    geom_sf(data = filter(df_year, Description == 2),
            aes(fill = factor(Description)), color = NA) +
    geom_sf(data = filter(df_year, Description == 1),
            aes(fill = factor(Description)), color = "black") +
    geom_sf(data = filter(df_year, Description == 0),
            aes(fill = factor(Description)), color = "black") +
    geom_shadowtext(data = df_labels,
                    aes(x = lon, y = lat, label = Individual),
                    size = 2, fontface = "bold",
                    color = "black", bg.color = "white", bg.r = 0.15) +
    geom_sf(data = study_area, fill = NA, color = "black", linewidth = 1) +
    scale_fill_manual(
      values = category_colors,
      labels = category_labels,
      guide = if (show_legend) guide_legend(title = "Territory Type") else "none"
    ) +
    annotation_scale(location = "br", width_hint = 0.3) +
    annotation_north_arrow(location = "tl", which_north = "true",
                           style = north_arrow_fancy_orienteering) +
    ggtitle(year_label) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0, vjust = 1.5, face = "bold"),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )
}

# Generate panel plots
p2019 <- plot_year(territories_cropped, "2019")
p2021 <- plot_year(territories_cropped, "2021")
p2022 <- plot_year(territories_cropped, "2022")
p2023 <- plot_year(territories_cropped, "2023")

# Load the required package
library(ggpubr)

# Create a dummy plot for the legend
p_legend <- ggplot() +
  geom_sf(data = filter(territories_cropped, Year == "2023"),
          aes(fill = factor(Description)), color = NA) +
  scale_fill_manual(
    values = category_colors,
    labels = category_labels,
    name = "Territory Type"
  ) +
  theme_void() +
  theme(legend.position = "bottom")

# Extract the legend using ggpubr
legend_patch <- patchwork::wrap_elements(full = ggpubr::get_legend(p_legend))



(p2019 | p2021 | p2022 | p2023) / legend_patch +
  plot_layout(heights = c(10, 1))

