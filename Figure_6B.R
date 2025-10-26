library(dplyr)
library(tidyr)
library(sf)
library(stringr)
library(purrr)
library(readr)


########### Data preparation

setwd("C:/Users/suyas/Desktop/WBS songs")

df_notes <- read.csv("Final_Classified_Notes_Dec16.csv") %>%
  filter(!paste(Individual, Year, sep = "-") %in% c(
    "BWGG-2022", "BWRU-2021", "GG__-2022", "RU__-2019", "RW__-2019", "YOUU-2021"
  )) %>%
  mutate(
    Ind_Year = paste(Individual, Year, sep = "-"),
    Year = as.character(Year)
  )

ngram_freqs <- df_notes %>%
  pivot_longer(cols = starts_with("Gram"), names_to = "Ngram_Level", values_to = "Ngram_ID") %>%
  filter(!is.na(Ngram_ID)) %>%
  group_by(Ind_Year, Year, Individual, Ngram_Level, Ngram_ID) %>%
  summarise(Freq = n(), .groups = "drop")

ngram_richness <- ngram_freqs %>%
  group_by(Ind_Year, Ngram_Level) %>%
  summarise(Richness = n_distinct(Ngram_ID), .groups = "drop")

# Create pairwise data within year
pair_data <- ngram_freqs %>%
  select(Year, Ngram_Level, Ind_Year) %>%
  distinct() %>%
  group_by(Year, Ngram_Level) %>%
  summarise(
    Pair_Matrix = list(as.data.frame(t(combn(Ind_Year, 2)))),
    .groups = "drop"
  ) %>%
  unnest(Pair_Matrix) %>%
  rename(Ind1_Year = V1, Ind2_Year = V2) %>%
  mutate(
    Pair_ID = map2_chr(Ind1_Year, Ind2_Year, ~ paste(sort(c(.x, .y)), collapse = "-"))
  )


# Join back repertoires and calculate values
pairwise_ngram <- pair_data %>%
  left_join(ngram_freqs, by = c("Ind1_Year" = "Ind_Year", "Ngram_Level", "Year")) %>%
  rename(Freq1 = Freq, Ngram_ID = Ngram_ID) %>%
  left_join(ngram_freqs, by = c("Ind2_Year" = "Ind_Year", "Ngram_Level", "Year", "Ngram_ID")) %>%
  rename(Freq2 = Freq) %>%
  mutate(
    Occurrence_Max = pmax(Freq1, Freq2, na.rm = TRUE),
    Shared = ifelse(!is.na(Freq1) & !is.na(Freq2), 1, 0)
  )


df_params <- read.csv("Ngram_Parameters_All.csv") %>%
  rename(Ngram_ID = Ngram_ID, Ngram_Level = Ngram_Level) %>%
  mutate(Ngram_Level = as.character(Ngram_Level))


pairwise_model_df <- pairwise_ngram %>%
  left_join(df_params, by = c("Ngram_ID", "Ngram_Level")) %>%
  mutate(
    Shared = factor(Shared, levels = c(0, 1)),
    Year = as.character(Year)
  )

# Replace with your actual file names for each year
kml_files <- list(
  "2019" = "Summer 2019.kml",
  "2021" = "Summer 2021.kml",
  "2022" = "Summer 2022.kml",
  "2023" = "Summer 2023.kml"
)

territory_centroids <- map_dfr(names(kml_files), function(yr) {
  sf_data <- st_read(kml_files[[yr]]) %>%
    mutate(Year = yr) %>%
    select(Individual = Name, Year, geometry)  # Keep Year before centroid
  centroids <- st_centroid(sf_data)
  centroids %>% mutate(Ind_Year = paste(Individual, Year, sep = "-"))
})


territory_dists <- territory_centroids %>%
  group_split(Year) %>%
  map_dfr(function(df) {
    d <- st_distance(df)
    combs <- t(combn(nrow(df), 2))
    tibble(
      Year = unique(df$Year),
      Ind1 = df$Ind_Year[combs[, 1]],
      Ind2 = df$Ind_Year[combs[, 2]],
      Dist_m = as.numeric(d[combs])
    )
  }) %>%
  mutate(
    Pair_ID = map2_chr(Ind1, Ind2, ~ paste(sort(c(.x, .y)), collapse = "-")),
    Is_Neighbor = Dist_m <= 200
  )


# Step 1: Get a list of neighbors for each individual-year
neighbor_list <- territory_dists %>%
  filter(Is_Neighbor) %>%
  mutate(
    Focal1 = Ind1,
    Neighbor1 = Ind2,
    Focal2 = Ind2,
    Neighbor2 = Ind1
  ) %>%
  select(Year, Focal1, Neighbor1, Focal2, Neighbor2) %>%
  pivot_longer(
    cols = c(Focal1, Focal2),
    names_to = "Role",
    values_to = "Focal"
  ) %>%
  mutate(
    Neighbor = ifelse(Role == "Focal1", Neighbor1, Neighbor2)
  ) %>%
  select(Year, Focal, Neighbor) %>%
  group_by(Year, Focal) %>%
  summarise(Neighbors = list(unique(Neighbor)), .groups = "drop")


# Step 2: Generate all pairs within a year
pairwise_neighbors <- neighbor_list %>%
  group_by(Year) %>%
  summarise(
    Pairs = list(as_tibble(t(combn(Focal, 2)), .name_repair = ~ c("Ind1", "Ind2"))),
    .groups = "drop"
  ) %>%
  unnest(Pairs)


# Step 3: Join neighbor lists to each pair
pairwise_neighbors <- pairwise_neighbors %>%
  left_join(neighbor_list, by = c("Year", "Ind1" = "Focal")) %>%
  rename(Neighbors1 = Neighbors) %>%
  left_join(neighbor_list, by = c("Year", "Ind2" = "Focal")) %>%
  rename(Neighbors2 = Neighbors)

# Step 4: Calculate union and intersection counts
pairwise_neighbors <- pairwise_neighbors %>%
  mutate(
    Pair_ID = map2_chr(Ind1, Ind2, ~ paste(sort(c(.x, .y)), collapse = "-")),
    UniqueNeighborCount = map2_int(Neighbors1, Neighbors2, ~ length(union(.x, .y))),
    SharedNeighborCount = map2_int(Neighbors1, Neighbors2, ~ length(intersect(.x, .y)))
  ) %>%
  select(Year, Ind1, Ind2, Pair_ID, UniqueNeighborCount, SharedNeighborCount)


# Join all predictors to pairwise_model_df
final_df <- pairwise_model_df %>%
  left_join(territory_dists %>% select(Pair_ID, Dist_m), by = "Pair_ID") %>%
  left_join(ngram_richness, by = c("Ind1_Year" = "Ind_Year", "Ngram_Level")) %>%
  rename(Richness1 = Richness) %>%
  left_join(ngram_richness, by = c("Ind2_Year" = "Ind_Year", "Ngram_Level")) %>%
  rename(Richness2 = Richness) %>%
  mutate(NgramRichness = Richness1 + Richness2) %>%
  left_join(pairwise_neighbors %>% select(Pair_ID, UniqueNeighborCount, SharedNeighborCount), by = "Pair_ID")

# Select relevant variables
final_df_clean <- final_df %>%
  select(
    Year, Ngram_Level, Pair_ID, Occurrence_Max, Shared,
    Mean_Position, Mean_Freq, Delta_Freq, Duration,
    Dist_m, Richness1, Richness2, NgramRichness,
    UniqueNeighborCount, SharedNeighborCount
  )

# Clean up workspace
rm(
  df_notes, df_params, final_df, kml_files, ngram_freqs,
  ngram_richness, pair_data, pairwise_model_df, pairwise_ngram,
  territory_centroids, territory_dists, neighbor_list, pairwise_neighbors
)

# Prepare data for GAM modeling
df_bin <- final_df_clean %>%
  filter(
    !is.na(Shared),
    !is.na(Mean_Position), !is.na(Mean_Freq), !is.na(Delta_Freq),
    !is.na(Duration), !is.na(UniqueNeighborCount), !is.na(SharedNeighborCount)
  ) %>%
  mutate(
    Shared = as.integer(Shared),                  # Ensure integer type
    Occurrence_Max = pmin(Occurrence_Max, 50) 
  )%>%
  mutate(Shared = ifelse(Shared == 2, 1, 0))

####### MODEL CODE
library(mgcv)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rlang)

# Fit GAMs for each N-gram level using only s() terms
gam_sharing_models <- list()

for (g in c("Gram1", "Gram2", "Gram3")) {
  df_g <- df_bin %>% filter(Ngram_Level == g)
  
  gam_sharing_models[[g]] <- gam(
    Shared ~ 
      s(Dist_m, k = 5) +
      s(SharedNeighborCount, k = 5) +
      s(NgramRichness, k = 5) +
      s(Occurrence_Max, k = 5) +
      s(Mean_Position, k = 5) +
      s(Mean_Freq, k = 5) +
      s(Delta_Freq, k = 5) +
      s(Duration, k = 5),
    data = df_g,
    family = binomial(link = "logit"),
    method = "REML"
  )
}


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


# Define your custom order and labels
order <- c("Dist_m", "SharedNeighborCount", "Occurrence_Max", "NgramRichness",
           "Mean_Position", "Mean_Freq", "Delta_Freq", "Duration")

rowlab <- c("Pairwise Distance between Territories",
            "Pairwise Number of Shared Neighbors",
            "Pairwise Maximum N-gram Occurrence",
            "Pairwise Cumulative N-gram Richness",
            "N-gram Mean Position in Song",
            "N-gram Mean Frequency",
            "N-gram Frequency Bandwidth",
            "N-gram Duration")

collabs <- c("1-Gram", "2-Gram", "3-Gram")

xlabs <- c("Distance (m)", "Number of Neighbors", "N-gram Occurrence",
           "N-gram Richness", "Mean Position", "Frequency (Hz)",
           "Frequency (Hz)", "Duration (s)")

# Generate prediction data from all models and smooth terms
plot_data_df <- map_dfr(names(gam_sharing_models), function(g) {
  model <- gam_sharing_models[[g]]
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

# Generate axis labels map
xlab_map <- setNames(xlabs, rowlab)

library(ggplot2)
library(patchwork)

# Define N-gram specific colors
ngram_colors <- c(
  "1" = "#280b53",
  "2" = "#9f2a63",
  "3" = "#f57d15"
)

# Split data by facet_label
split_plots <- split(plot_data_df, plot_data_df$facet_label)

# Generate x-axis labels (already defined earlier)
xlab_map <- setNames(xlabs, rowlab)

# Create individual plots
facet_plots <- lapply(names(split_plots), function(facet_label) {
  df <- split_plots[[facet_label]]
  row_label <- strsplit(facet_label, " \\| ")[[1]][1]
  col_label <- strsplit(facet_label, " \\| ")[[1]][2]
  xlab <- xlab_map[[row_label]]
  
  # Extract N-gram number from col_label
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

# Combine with patchwork (3 columns)
final_plot <- wrap_plots(facet_plots, ncol = 3) +
  plot_annotation(title = "Predicted Probability of N-gram Sharing")

final_plot

