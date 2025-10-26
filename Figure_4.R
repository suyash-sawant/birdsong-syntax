# Load required libraries
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(patchwork)
library(grid)
library(RColorBrewer)
library(ggpubr)
library(vegan)

#--------------------------
# Global Setup
#--------------------------
fill_colors <- c("Neighbor" = "#74c476", "Non-neighbor" = "#ef6548")
years_colors <- c("2019" = "#a6bddb", "2020" = "#74a9cf", "2021" = "#2b8cbe",
                  "2022" = "#045a8d", "2023" = "#023858")

# Load main data
df <- read.csv("Final_Classified_Notes_Dec16.csv") %>%
  mutate(
    Ind_Year = paste(Individual, Year, sep = "-"),
    Year = as.character(Year)
  ) %>%
  filter(!Ind_Year %in% c(
    "BWGG-2022", "BWRU-2021", "GG__-2022", "RU__-2019",
    "RW__-2019", "YOUU-2021"
  ))

# Load pairwise distance data
df1 <- read.csv("WBS_Final_Data - Ter_Dist.csv") %>%
  mutate(
    Ind1_Year = paste(Ind1, Year, sep = "-"),
    Ind2_Year = paste(Ind2, Year, sep = "-"),
    Pair_ID = map2_chr(Ind1_Year, Ind2_Year, ~ paste(sort(c(.x, .y)), collapse = "-"))
  )

#--------------------------
# Function: N-gram Overlap by Pair (Within-Year)
#--------------------------
ngram_overlap_by_pair <- function(df, n) {
  gram_col <- paste0("Gram", n)
  
  df_filtered <- df %>%
    filter(!is.na(.data[[gram_col]])) %>%
    group_by(Ind_Year, Year, Ngram = .data[[gram_col]]) %>%
    filter(n() > 1) %>%
    ungroup()
  
  rep_sets <- df_filtered %>%
    distinct(Ind_Year, Year, Ngram) %>%
    group_by(Ind_Year, Year) %>%
    summarise(Repertoire = list(Ngram), .groups = "drop")
  
  if (nrow(rep_sets) < 2) return(NULL)
  
  pair_df <- rep_sets %>%
    group_by(Year) %>%
    summarise(Combos = list(t(combn(Ind_Year, 2))), .groups = "drop") %>%
    unnest(Combos) %>%
    mutate(
      Ind1_Year = Combos[, 1],
      Ind2_Year = Combos[, 2],
      Pair_ID = map2_chr(Ind1_Year, Ind2_Year, ~ paste(sort(c(.x, .y)), collapse = "-"))
    ) %>%
    select(Pair_ID, Ind1_Year, Ind2_Year, Year) %>%
    left_join(rep_sets, by = c("Ind1_Year" = "Ind_Year", "Year")) %>%
    rename(Rep1 = Repertoire) %>%
    left_join(rep_sets, by = c("Ind2_Year" = "Ind_Year", "Year")) %>%
    rename(Rep2 = Repertoire) %>%
    rowwise() %>%
    mutate(
      Overlap = length(intersect(Rep1, Rep2)) / length(union(Rep1, Rep2)),
      Ngram = n
    ) %>%
    ungroup() %>%
    select(Pair_ID, Year, Ngram, Overlap)
}

# Calculate N-gram overlaps
overlap_df <- map_dfr(1:3, ~ngram_overlap_by_pair(df, .x))

# Merge with spatial distance
df_combined <- overlap_df %>%
  left_join(df1 %>% select(Pair_ID, Distance_m), by = "Pair_ID") %>%
  filter(!is.na(Distance_m)) %>%
  mutate(
    Neighbor_Status = factor(ifelse(Distance_m <= 200, "Neighbor", "Non-neighbor"),
                             levels = c("Neighbor", "Non-neighbor"))
  )

#--------------------------
# ROW A: Boxplots
#--------------------------
rowA_boxplots <- map(1:3, function(n) {
  df_n <- filter(df_combined, Ngram == n)
  
  p_val <- wilcox.test(Overlap ~ Neighbor_Status, data = df_n)$p.value
  stars <- case_when(
    p_val < 0.001 ~ "***",
    p_val < 0.01 ~ "**",
    p_val < 0.05 ~ "*",
    TRUE ~ "n.s."
  )
  
  y_max <- max(df_n$Overlap, na.rm = TRUE)
  
  ggplot(df_n, aes(x = Neighbor_Status, y = Overlap)) +
    geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = 16, outlier.size = 1.5) +
    #geom_jitter(width = 0.15, alpha = 0.4, size = 1.5, color = "gray30") +
    annotate("text", x = 1.5, y = y_max * 1.05, label = stars, size = 5) +
    labs(
      title = paste0(n, "-Gram"),
      x = "", y = ifelse(n == 1, "N-gram Overlap", "")
    ) +
    theme_bw() +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))
})

rowA <- wrap_plots(rowA_boxplots, ncol = 3)

#--------------------------
# ROW B: Cumulative Richness
#--------------------------
get_accum_by_year_ngram <- function(df, gram_col, year) {
  df_sub <- df %>%
    filter(Year == year, !is.na(.data[[gram_col]])) %>%
    mutate(Ind_Year = paste(Individual, Year, sep = "-")) %>%
    distinct(Ind_Year, Token = .data[[gram_col]]) %>%
    mutate(presence = 1) %>%
    pivot_wider(names_from = Token, values_from = presence, values_fill = 0) %>%
    column_to_rownames("Ind_Year")
  
  if (nrow(df_sub) < 2) return(NULL)
  
  acc <- vegan::specaccum(df_sub, method = "random", permutations = 100)
  tibble(
    N_Individuals = acc$sites,
    Richness = acc$richness,
    SD = acc$sd,
    Year = year,
    Ngram = gram_col
  )
}

rich_df <- map_dfr(c("Gram1", "Gram2", "Gram3"), function(gram) {
  map_dfr(c("2021", "2022"), function(yr) {
    get_accum_by_year_ngram(df, gram, yr)
  })
}) %>%
  mutate(Ngram_Label = recode(Ngram, Gram1 = "1-Gram", Gram2 = "2-Gram", Gram3 = "3-Gram"))

rowB_rich <- map(unique(rich_df$Ngram_Label), function(label) {
  ggplot(filter(rich_df, Ngram_Label == label),
         aes(x = N_Individuals, y = Richness, color = Year, fill = Year)) +
    geom_line(linewidth = 0.9) +
    geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), alpha = 0.25, color = NA) +
    scale_color_manual(values = years_colors) +
    scale_fill_manual(values = years_colors) +
    labs(title = label, x = "Number of Individuals",
         y = ifelse(label == "1-Gram", "Cumulative N-gram Richness", "")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
})

legend_B <- get_legend(
  ggplot(rich_df, aes(x = N_Individuals, y = Richness, color = Year, fill = Year)) +
    geom_line() + geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), alpha = 0.3, color = NA) +
    scale_color_manual(values = years_colors) +
    scale_fill_manual(values = years_colors) +
    theme(legend.position = "top")
)

rowB <- wrap_plots(rowB_rich, ncol = 3)

#--------------------------
# ROW C: N-gram Sharing
#--------------------------
ngram_share_df <- df %>%
  filter(Year %in% c("2021", "2022")) %>%
  pivot_longer(cols = all_of(c("Gram1", "Gram2", "Gram3")),
               names_to = "Ngram", values_to = "Token") %>%
  filter(!is.na(Token)) %>%
  distinct(Year, Individual, Ngram, Token) %>%
  group_by(Year, Ngram, Token) %>%
  summarise(Num_Individuals = n(), .groups = "drop") %>%
  group_by(Year, Ngram, Num_Individuals) %>%
  summarise(Ngram_Count = n(), .groups = "drop") %>%
  group_by(Year, Ngram) %>%
  mutate(Percent = Ngram_Count / sum(Ngram_Count) * 100) %>%
  ungroup() %>%
  mutate(Ngram_Label = recode(Ngram,
                              Gram1 = "1-Gram",
                              Gram2 = "2-Gram",
                              Gram3 = "3-Gram"))


rowC_share <- map(unique(ngram_share_df$Ngram_Label), function(label) {
  ggplot(filter(ngram_share_df, Ngram_Label == label),
         aes(x = Num_Individuals, y = Percent, color = Year)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.5) +
    scale_color_manual(values = years_colors) +
    ylim(0, 100) +
    labs(title = label, x = "Individuals Sharing",
         y = ifelse(label == "1-Gram", "Percentage of N-grams", "")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
})

legend_C <- get_legend(
  ggplot(ngram_share_df, aes(x = Num_Individuals, y = Percent, color = Year)) +
    geom_line() + geom_point() +
    scale_color_manual(values = years_colors) +
    theme(legend.position = "top")
)

rowC <- wrap_plots(rowC_share, ncol = 3)

#--------------------------
# Final Assembly
#--------------------------
rowA_title <- wrap_elements(grid::textGrob("(C) Neighbor vs Non-neighbor", x = 0, hjust = 0, gp = gpar(fontsize = 12, fontface = "bold")))
rowB_title <- wrap_elements(grid::textGrob("(A) Cumulative N-gram Richness", x = 0, hjust = 0, gp = gpar(fontsize = 12, fontface = "bold")))
rowC_title <- wrap_elements(grid::textGrob("(B) N-gram Sharing Across Individuals", x = 0, hjust = 0, gp = gpar(fontsize = 12, fontface = "bold")))

final_plot <- (
  rowB_title / wrap_elements(full = legend_B) / rowB /
    rowC_title / wrap_elements(full = legend_C) / rowC /
    rowA_title / rowA
) +
  plot_layout(heights = c(0.15, 0.15, 1.2, 0.15, 0.15, 1.2, 0.15, 1.2))

# Display
final_plot

