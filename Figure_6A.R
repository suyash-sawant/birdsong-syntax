library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(sf)
library(readr)

# Set working directory
setwd("C:/Users/suyas/Desktop/WBS songs")

# Load and filter notes data (exclude 2019 and unwanted Ind_Years)
df_notes <- read.csv("Final_Classified_Notes_Dec16.csv") %>%
  filter(!Year %in% 2019) %>%
  filter(!paste(Individual, Year, sep = "-") %in% c("BWGG-2022", "BWRU-2021", "GG__-2022", "RU__-2019", "RW__-2019", "YOUU-2021")) %>%
  mutate(Ind_Year = paste(Individual, Year, sep = "-"))

# Only keep Gram1 to Gram3
df_ngrams <- df_notes %>%
  pivot_longer(cols = starts_with("Gram"), names_to = "Gram_Level", values_to = "Ngram_ID") %>%
  filter(Gram_Level %in% c("Gram1", "Gram2", "Gram3"), !is.na(Ngram_ID))

# Frequency table per Ngram per individual/year
ngram_freqs <- df_ngrams %>%
  group_by(Individual, Year, Gram_Level, Ngram_ID) %>%
  summarise(Freq = n(), .groups = "drop") %>%
  mutate(Year = as.integer(as.character(Year)))

# Generate year pairs (consecutive only)
year_pairs <- expand.grid(Year1 = 2021:2022, Year2 = 2022:2023) %>%
  filter(Year2 == Year1 + 1)


# Get all N-grams in Year 1 of the pair
ngrams_y1 <- ngram_freqs %>%
  inner_join(year_pairs, by = c("Year" = "Year1")) %>%
  rename(Year1 = Year, Occurrence = Freq) %>%
  select(Individual, Year1, Year2, Gram_Level, Ngram_ID, Occurrence)

# Mark consistent if present in both years
consistent_ngrams <- ngrams_y1 %>%
  inner_join(ngram_freqs %>% rename(Year2 = Year), 
             by = c("Individual", "Gram_Level", "Ngram_ID", "Year2")) %>%
  mutate(Status = 1)  # consistent

# Mark inconsistent if present in Year1 but *not* in Year2
inconsistent_ngrams <- ngrams_y1 %>%
  anti_join(ngram_freqs %>% rename(Year2 = Year), 
            by = c("Individual", "Gram_Level", "Ngram_ID", "Year2")) %>%
  mutate(Status = 0)  # inconsistent

# Combine
ngram_consistency <- bind_rows(consistent_ngrams, inconsistent_ngrams)
head(ngram_consistency)

# Ngram richness (repertoire size) in Year1
ngram_richness_y1 <- ngram_freqs %>%
  group_by(Individual, Year, Gram_Level) %>%
  summarise(NgramRichness = n_distinct(Ngram_ID), .groups = "drop") %>%
  rename(Year1 = Year)

# Import acoustic parameters
df_params <- read.csv("Ngram_Parameters_All.csv") %>%
  filter(Ngram_Level %in% c("Gram1", "Gram2", "Gram3")) %>%
  rename(Gram_Level = Ngram_Level)


kml_files <- list(
  "2021" = "Summer 2021.kml",
  "2022" = "Summer 2022.kml",
  "2023" = "Summer 2023.kml"
)

territory_sf <- map_dfr(names(kml_files), function(yr) {
  st_read(kml_files[[yr]]) %>%
    mutate(Year = as.integer(yr)) %>%
    select(Individual = Name, Year, geometry)
})

# Calculate territory area
territory_area <- territory_sf %>%
  mutate(Area = as.numeric(st_area(geometry))) %>%
  st_drop_geometry()

# Calculate delta territory area
territory_change <- territory_area %>%
  rename(Year1 = Year, Area1 = Area) %>%
  inner_join(territory_area %>% rename(Year2 = Year, Area2 = Area),
             by = "Individual") %>%
  filter(Year2 == Year1 + 1) %>%
  mutate(DeltaTerritoryArea = (Area2 - Area1) / Area1) %>%
  select(Individual, Year1, Year2, DeltaTerritoryArea)

# Recompute centroids for neighbors
territory_centroids <- territory_sf %>%
  st_centroid() %>%
  mutate(Ind_Year = paste(Individual, Year, sep = "-"))

# Neighbor calculation
neighbor_union <- territory_centroids %>%
  group_split(Year) %>%
  map_dfr(function(df) {
    d <- st_distance(df)
    combs <- t(combn(nrow(df), 2))
    
    tibble(
      Year = unique(df$Year),
      Ind1 = df$Individual[combs[, 1]],
      Ind2 = df$Individual[combs[, 2]],
      Dist_m = as.numeric(d[combs])
    )
  }) %>%
  filter(Dist_m <= 200) %>%
  # Create two directional entries (Ind1 → Ind2 and Ind2 → Ind1)
  mutate(
    Focal1 = Ind1,
    Neighbor1 = Ind2,
    Focal2 = Ind2,
    Neighbor2 = Ind1
  ) %>%
  select(Year,
         Focal1, Neighbor1,
         Focal2, Neighbor2) %>%
  pivot_longer(
    cols = c(Focal1, Focal2, Neighbor1, Neighbor2),
    names_to = c(".value", "group"),
    names_pattern = "(Focal|Neighbor)([12])"
  ) %>%
  group_by(Year, Focal) %>%
  summarise(Neighbors = list(unique(Neighbor)), .groups = "drop")

# Join neighbors across consecutive years
neighbor_combined <- neighbor_union %>%
  rename(Year1 = Year, Neighbors1 = Neighbors, Individual = Focal) %>%
  inner_join(
    neighbor_union %>% rename(Year2 = Year, Neighbors2 = Neighbors, Individual = Focal),
    by = "Individual"
  ) %>%
  filter(Year2 == Year1 + 1) %>%
  mutate(
    AllNeighbors = map2(Neighbors1, Neighbors2, union),
    NeighborCount = map_int(AllNeighbors, length)
  ) %>%
  select(Individual, Year1, Year2, NeighborCount)




# Final merged data
df_consistency <- ngram_consistency %>%
  left_join(ngram_richness_y1, by = c("Individual", "Gram_Level", "Year1")) %>%
  left_join(territory_change, by = c("Individual", "Year1", "Year2")) %>%
  left_join(neighbor_combined, by = c("Individual", "Year1", "Year2")) %>%
  left_join(df_params, by = c("Gram_Level", "Ngram_ID")) %>%
  filter(!is.na(NgramRichness), !is.na(NeighborCount), !is.na(DeltaTerritoryArea))


summary(df_consistency)

library(mgcv)

gam_consistency_models <- list()

for (g in c("Gram1","Gram2","Gram3")) {
  df_g <- df_consistency %>% filter(Gram_Level == g)
  
  gam_consistency_models[[g]] <- gam(
    Status ~ 
      s(NeighborCount, k = 5) +
      s(DeltaTerritoryArea, k = 5) +
      s(NgramRichness, k = 5) +
      s(Occurrence, k = 5) +
      s(Mean_Position, k = 5) +
      s(Mean_Freq, k = 5) +
      s(Delta_Freq, k = 5) +
      s(Duration, k = 5),
    data = df_g,
    family = binomial(link = "logit"),
    method = "REML"
  )
}

summary(gam_consistency_models[[g]])


#---
library(mgcv)
library(dplyr)
library(tidyr)
library(purrr)

# Define predictors to extract
smooth_terms <- c("DeltaTerritoryArea","NeighborCount","Occurrence", "NgramRichness",
                  "Mean_Position", "Mean_Freq", "Delta_Freq", "Duration")


# Function to get prediction data for one model term
get_gam_smooth_data <- function(model, term, gram_label, n = 200) {
  x_seq <- seq(min(model$model[[term]], na.rm = TRUE),
               max(model$model[[term]], na.rm = TRUE),
               length.out = n)
  
  new_data <- model$model[1, , drop = FALSE]
  new_data <- new_data[rep(1, n), ]
  new_data[[term]] <- x_seq
  
  for (col in names(new_data)) {
    if (col != term && is.numeric(new_data[[col]])) {
      new_data[[col]] <- median(model$model[[col]], na.rm = TRUE)
    }
  }
  
  preds <- predict(model, newdata = new_data, se.fit = TRUE, type = "link")
  
  # Convert from logit to probability scale
  fit_prob <- plogis(preds$fit)
  lower_prob <- plogis(preds$fit - 2 * preds$se.fit)
  upper_prob <- plogis(preds$fit + 2 * preds$se.fit)
  
  tibble(
    Ngram = gram_label,
    Predictor = term,
    x_value = x_seq,
    fit = fit_prob,
    lower = lower_prob,
    upper = upper_prob
  )
}



# Order and labeling
order <- c("DeltaTerritoryArea", "NeighborCount", "Occurrence", "NgramRichness",
           "Mean_Position", "Mean_Freq", "Delta_Freq", "Duration")

rowlab <- c("Proportional Change in Territory Area",
            "Total Number of Neighbors",
            "N-gram Occurrence in Year 1",
            "N-gram Richness in Year 1",
            "N-gram Mean Position in Song",
            "N-gram Mean Frequency",
            "N-gram Frequency Bandwidth",
            "N-gram Duration")

collabs <- c("1-Gram", "2-Gram", "3-Gram")

xlabs <- c("Proportional Change", "Number of Neighbors", 
           "N-gram Occurrence", "Ngram Richness", 
           "Mean Position", "Frequency (Hz)", "Bandwidth (Hz)", "Duration (s)")

# Generate prediction data from all models and smooth terms
plot_data_df <- map_dfr(names(gam_consistency_models), function(g) {
  model <- gam_consistency_models[[g]]
  map_dfr(smooth_terms, function(term) {
    get_gam_smooth_data(model, term, gram_label = g)
  })
})


# Add clean Ngram label
plot_data_df <- plot_data_df %>%
  mutate(
    Predictor = factor(Predictor, levels = order),
    Predictor_Label = factor(Predictor, labels = rowlab),
    Ngram = factor(Ngram, levels = c("Gram1", "Gram2", "Gram3"),
                   labels = collabs),
    facet_label = interaction(Predictor_Label, Ngram, sep = " | ", lex.order = TRUE)
  )

library(ggplot2)
library(patchwork)

ngram_colors <- c("1" = "#280b53", "2" = "#9f2a63", "3" = "#f57d15")
xlab_map <- setNames(xlabs, rowlab)

# Split data
split_plots <- split(plot_data_df, plot_data_df$facet_label)

facet_plots <- lapply(names(split_plots), function(facet_label) {
  df <- split_plots[[facet_label]]
  row_label <- strsplit(facet_label, " \\| ")[[1]][1]
  col_label <- strsplit(facet_label, " \\| ")[[1]][2]
  xlab <- xlab_map[[row_label]]
  ngram_num <- gsub("-Gram", "", col_label)
  line_color <- ngram_colors[[ngram_num]]
  
  ggplot(df, aes(x = x_value, y = fit)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = line_color, alpha = 0.2) +
    geom_line(color = line_color, size = 1) +
    labs(x = xlab, y = NULL) +
    theme_bw(base_size = 11) +
    ggtitle(facet_label) +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 9),
      axis.text = element_text(size = 8),
      strip.background = element_blank()
    )
})

final_plot_consistency <- wrap_plots(facet_plots, ncol = 3) +
  plot_annotation(title = "Predicted Probability of Within-Individual N-gram Consistency")

final_plot_consistency

