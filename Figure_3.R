library(dplyr)
library(stringdist)
library(tidyr)
library(purrr)
library(lubridate)
library(ggplot2)
library(patchwork)

# --- Individual filter ---
individual_order <- c("YYXX", "YOXX", "RUGG", "YOGG")

# --- Load and clean data ---
setwd("C:/Users/suyas/Desktop/WBS songs")
df <- read.csv("Final_Classified_Notes_Dec16.csv") %>%
  mutate(
    Individual = gsub("__", "XX", Individual),
    Year = as.character(Year),
    Note_DateTime = dmy_hms(paste(Note_Date, Note_Time), tz = "UTC")
  ) %>%
  filter(Individual %in% individual_order)

# --- Song start times ---
song_times <- df %>%
  group_by(Individual, Year, Note_Date, Song_Num) %>%
  summarise(Song_Time = min(Note_DateTime), .groups = "drop")

# --- Faster short-term pairing ---
get_song_pairs_fast <- function(sub_df, max_songs = 30) {
  sub_df %>%
    group_by(Note_Date, Song_Num, Begin_File) %>%
    slice_min(Note_Num, n = 1) %>%
    ungroup() %>%
    group_by(Note_Date, Individual, Year) %>%  # include grouping vars here
    group_modify(~{
      songs <- .x
      if (nrow(songs) < 2) {
        return(tibble(
          Song1 = numeric(0),
          Song2 = numeric(0),
          Begin_File1 = character(0),
          Begin_File2 = character(0)
        ))
      }
      
      songs <- slice_sample(songs, n = min(max_songs, nrow(songs)))
      combs <- t(combn(songs$Song_Num, 2))
      tibble(
        Song1 = combs[, 1],
        Song2 = combs[, 2],
        Begin_File1 = songs$Begin_File[match(combs[, 1], songs$Song_Num)],
        Begin_File2 = songs$Begin_File[match(combs[, 2], songs$Song_Num)]
      )
    }) %>%
    ungroup()
}




# --- Faster long-term pairing ---
get_long_term_pairs_fast <- function(sub_df, max_songs = 100) {
  songs <- sub_df %>%
    group_by(Song_Num, Note_Date, Begin_File) %>%
    slice_min(Note_Num, n = 1) %>%
    ungroup() %>%
    distinct(Song_Num, Note_Date, Begin_File)
  
  if (nrow(songs) < 2) return(tibble())
  
  songs <- slice_sample(songs, n = min(max_songs, nrow(songs)))
  combs <- t(combn(songs$Song_Num, 2))
  
  song_time_lookup <- song_times %>%
    filter(Individual == sub_df$Individual[1], Year == sub_df$Year[1]) %>%
    select(Song_Num, Song_Time)
  
  left <- tibble(Song1 = combs[,1]) %>%
    left_join(song_time_lookup, by = c("Song1" = "Song_Num")) %>%
    rename(Song1_Time = Song_Time)
  
  right <- tibble(Song2 = combs[,2]) %>%
    left_join(song_time_lookup, by = c("Song2" = "Song_Num")) %>%
    rename(Song2_Time = Song_Time)
  
  tibble(
    Individual = sub_df$Individual[1],
    Year = sub_df$Year[1],
    Song1 = combs[,1],
    Song2 = combs[,2],
    Song1_Time = left$Song1_Time,
    Song2_Time = right$Song2_Time
  ) %>%
    mutate(Time_Diff_Day = abs(as.numeric(difftime(Song1_Time, Song2_Time, units = "days")))) %>%
    drop_na()
}

# --- Levenshtein Distance Calculator ---
calc_lev <- function(s1, s2, meta) {
  seq1 <- df %>%
    filter(Song_Num == s1, Individual == meta$Individual, Year == meta$Year) %>%
    arrange(Note_Num) %>%
    pull(Gram1)
  seq2 <- df %>%
    filter(Song_Num == s2, Individual == meta$Individual, Year == meta$Year) %>%
    arrange(Note_Num) %>%
    pull(Gram1)
  if (length(seq1) < 1 || length(seq2) < 1) return(NA_real_)
  all_notes <- unique(c(seq1, seq2))
  note_pool <- c(LETTERS, letters, as.character(0:9), "!", "@", "#", "$", "%", "&", "*", "+", "=", "~")
  if (length(all_notes) > length(note_pool)) return(NA_real_)
  note_map <- setNames(note_pool[1:length(all_notes)], all_notes)
  s1_str <- paste(note_map[seq1], collapse = "")
  s2_str <- paste(note_map[seq2], collapse = "")
  if (nchar(s1_str) == 0 || nchar(s2_str) == 0) return(NA_real_)
  stringdist(s1_str, s2_str, method = "lv") / max(nchar(s1_str), nchar(s2_str))
}

# --- Wrappers ---
get_dist_df_long <- function(df, time_col, scale_label) {
  df %>%
    mutate(Time = .data[[time_col]], Scale = scale_label) %>%
    pmap_dfr(function(Song1, Song2, Individual, Year, Time, Scale) {
      dist <- tryCatch(
        calc_lev(Song1, Song2, list(Individual = Individual, Year = Year)),
        error = function(e) NA_real_
      )
      tibble(Levenshtein_Dist = dist, Time = Time, Scale = Scale, Individual = Individual)
    }) %>%
    drop_na(Levenshtein_Dist)
}

# --- Sampling function ---
sample_intervals <- function(df, time_var, breaks, n_per_bin) {
  df %>%
    mutate(Bin = cut(.data[[time_var]], breaks = breaks)) %>%
    group_by(Bin) %>%
    group_modify(~ slice_sample(.x, n = min(nrow(.x), n_per_bin))) %>%
    ungroup()
}

# --- Generate short/long pairs ---
short_pairs <- df %>%
  group_by(Individual, Year) %>%
  group_split() %>%
  map_dfr(~get_song_pairs_fast(.x, max_songs = 30)) %>%
  left_join(song_times, by = c("Individual", "Year", "Note_Date", "Song1" = "Song_Num")) %>%
  rename(Song1_Time = Song_Time) %>%
  left_join(song_times, by = c("Individual", "Year", "Note_Date", "Song2" = "Song_Num")) %>%
  rename(Song2_Time = Song_Time) %>%
  mutate(Time_Diff_Sec = abs(as.numeric(difftime(Song1_Time, Song2_Time, units = "secs")))) %>%
  drop_na()


long_pairs <- df %>%
  group_by(Individual, Year) %>%
  group_split() %>%
  map_dfr(~get_long_term_pairs_fast(.x, max_songs = 100))

short_sampled <- sample_intervals(short_pairs, "Time_Diff_Sec", c(0,4,15,60,240,900,3600,14400), 1000)
long_sampled  <- sample_intervals(long_pairs, "Time_Diff_Day", c(1, 2, 4, 7, 14, 30, 60, 120), 1000)

# --- Distance calculation ---
lev_short <- short_sampled %>%
  mutate(Time = Time_Diff_Sec) %>%
  select(Song1, Song2, Individual, Year, Note_Date, Begin_File1, Begin_File2, Song1_Time, Song2_Time, Time) %>%
  pmap_dfr(function(Song1, Song2, Individual, Year, Note_Date, Begin_File1, Begin_File2, Song1_Time, Song2_Time, Time) {
    dist <- tryCatch(
      calc_lev(Song1, Song2, list(
        Individual = Individual,
        Year = Year,
        Note_Date = Note_Date,
        Begin_File1 = Begin_File1,
        Begin_File2 = Begin_File2
      )),
      error = function(e) NA_real_
    )
    tibble(Levenshtein_Dist = dist, Time = Time, Individual = Individual)
  }) %>%
  drop_na(Levenshtein_Dist)


lev_long <- long_sampled %>%
  rename(Time = Time_Diff_Day) %>%
  select(Song1, Song2, Individual, Year, Song1_Time, Song2_Time, Time) %>%
  pmap_dfr(function(Song1, Song2, Individual, Year, Song1_Time, Song2_Time, Time) {
    dist <- tryCatch(
      calc_lev(Song1, Song2, list(
        Individual = Individual,
        Year = Year
      )),
      error = function(e) NA_real_
    )
    tibble(Levenshtein_Dist = dist, Time = Time, Individual = Individual)
  }) %>%
  drop_na(Levenshtein_Dist)

# --- Across year ---
get_across_year_pairs <- function(df, year_gap, max_per_gap = 100) {
  map_dfr(unique(df$Individual), function(ind) {
    sub <- df %>% filter(Individual == ind)
    yrs <- sort(unique(as.numeric(sub$Year)))
    pairs <- list()
    for (y in yrs) {
      if ((y + year_gap) %in% yrs) {
        s1 <- sub %>% filter(Year == y) %>% distinct(Song_Num) %>% pull()
        s2 <- sub %>% filter(Year == as.character(y + year_gap)) %>% distinct(Song_Num) %>% pull()
        if (length(s1) > 0 && length(s2) > 0) {
          sampled <- expand.grid(Song1 = s1, Song2 = s2) %>%
            slice_sample(n = min(nrow(.), max_per_gap))
          sampled$Individual <- ind
          sampled$Year1 <- as.character(y)
        }
        pairs[[length(pairs) + 1]] <- sampled
      }
    }
    bind_rows(pairs)
  })
}

get_dist_df_yeargap <- function(df_pairs, gap) {
  df_pairs %>%
    mutate(Time = gap, Scale = paste0(gap, " year")) %>%
    rowwise() %>%
    mutate(Levenshtein_Dist = tryCatch(
      calc_lev(Song1, Song2, list(Individual = Individual, Year = Year1)),
      error = function(e) NA_real_
    )) %>%
    ungroup() %>%
    drop_na(Levenshtein_Dist)
}

yeargap_data <- map_dfr(1:4, function(gap) {
  across_pairs <- get_across_year_pairs(df, gap, max_per_gap = 100)
  get_dist_df_yeargap(across_pairs, gap)
})

# --- Plot ---
library(RColorBrewer)
individual_order <- c("YYXX", "YOXX", "RUGG", "YOGG")
colors <- setNames(brewer.pal(n = length(individual_order), name = "Set2"), individual_order)
shapes <- c("YYXX" = 16, "YOXX" = 17, "RUGG" = 15, "YOGG" = 8)

p1 <- ggplot(lev_short, aes(x = Time, y = Levenshtein_Dist)) +
  geom_point(size = 1.5, aes(color = Individual, shape = Individual)) +
  geom_smooth(se = T, color = "black", fill = "gray20") +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  scale_x_log10(
    limits = c(4, 14400),
    breaks = c(4,15,60,240,900,3600,14400),
    labels = c("4 sec","15 sec", "1 min", "4 min", "15 min", "1 hour", "4 hour")
  ) +
  labs(x = "Time gap", y = "Normalized Levenshtein Distance (NLD)", color = "Individual") +
  theme_bw() + ylim(0.3, 1) +
  theme(legend.position = "",
        axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(lev_long, aes(x = Time, y = Levenshtein_Dist)) +
  geom_point(size = 1.5, aes(color = Individual, shape = Individual)) +
  geom_smooth(se = T, color = "black", fill = "gray20") +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  scale_x_log10(
    limits = c(1, 120),
    breaks = c(1, 2, 4, 7, 14, 30, 60, 120),
    labels = c("1 day","2 day", "4 day", "1 week", "2 week", "1 month", "2 month", "4 month")
  ) +
  labs(x = "Time gap", y = NULL) +
  theme_bw() + ylim(0.3, 1) +
  theme(legend.position = "",
        axis.text.x = element_text(angle = 45, hjust = 1))

p3 <- ggplot(yeargap_data, aes(x = Time, y = Levenshtein_Dist)) +
  geom_point(position = position_jitter(width = 0.1, height = 0),
             size = 1.5, aes(color = Individual, shape = Individual))+
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = TRUE,
              color = "black", fill = "gray20") +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  scale_x_continuous(breaks = 1:4, labels = paste0(1:4, " year")) +
  labs(x = "Time gap", y = NULL) +
  theme_bw() + ylim(0.3, 1) +
  theme(legend.position = "",
        axis.text.x = element_text(angle = 45, hjust = 1))


#-------------------------------------------------------------------------------


library(tidyverse)
library(lubridate)
library(patchwork)

# N-gram color palette
ngram_colors <- c("Gram1" = "#000004", "Gram2" = "#8c2981", "Gram3" = "#fe9f6d")

# Individual color palette
individual_order <- c("YYXX", "YOXX", "RUGG", "YOGG")
colors <- setNames(RColorBrewer::brewer.pal(4, "Set2"), individual_order)
shapes <- c("YYXX" = 16, "YOXX" = 17, "RUGG" = 15, "YOGG" = 8)

# Load and clean
df <- read.csv("Final_Classified_Notes_Dec16.csv") %>%
  mutate(
    Individual = gsub("__", "XX", Individual),
    Year = as.character(Year),
    Note_DateTime = dmy_hms(paste(Note_Date, Note_Time), tz = "UTC")
  ) %>%
  filter(Individual %in% individual_order)

df$Individual <- factor(df$Individual, levels = individual_order)
df$Year <- factor(df$Year, levels = as.character(2019:2023))

# Function to compute cumulative unique ngrams
get_cumulative_unique <- function(df, gram_var) {
  df %>%
    filter(!is.na(.data[[gram_var]])) %>%
    distinct(Individual, Year, .data[[gram_var]]) %>%
    group_by(Individual, Year) %>%
    summarise(Gram_Set = list(unique(.data[[gram_var]])), .groups = "drop") %>%
    arrange(Individual, Year) %>%
    group_by(Individual) %>%
    mutate(
      Cumulative_Set = accumulate(Gram_Set, union),
      Cumulative_Count = lengths(Cumulative_Set),
      Ngram = gram_var
    ) %>%
    select(Individual, Year, Ngram, Cumulative_Count)
}

# Calculate for Gram1, Gram2, Gram3
cumulative_df <- bind_rows(
  get_cumulative_unique(df, "Gram1"),
  get_cumulative_unique(df, "Gram2"),
  get_cumulative_unique(df, "Gram3")
) %>%
  mutate(Year_num = as.integer(as.character(Year)))

cumulative_segments <- cumulative_df %>%
  arrange(Individual, Ngram, Year_num) %>%
  group_by(Individual, Ngram) %>%
  mutate(
    Year_from = lag(Year_num),
    Count_from = lag(Cumulative_Count)
  ) %>%
  filter(!is.na(Year_from)) %>%
  mutate(
    Line_Type = if_else(Year_from == 2019 & Year_num == 2021, "dashed", "solid")
  ) %>%
  ungroup()

plot_ngram_panel <- function(ngram_level) {
  df_points <- cumulative_df    %>% filter(Ngram == ngram_level)
  df_lines  <- cumulative_segments %>% filter(Ngram == ngram_level)
  
  ggplot() +
    geom_segment(
      data = df_lines,
      aes(x    = Year_from,  xend = Year_num,
          y    = Count_from, yend = Cumulative_Count,
          color = Individual,
          linetype = Line_Type),   # keep linetype mapping
      size = 1
    ) +
    geom_point(
      data  = df_points,
      aes(x = Year_num, y = Cumulative_Count, color = Individual, shape = Individual),
      size = 2, alpha = 0.6
    ) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    scale_linetype_manual(values = c(solid = "solid", dashed = "dashed"),
                          guide  = "none")+
    scale_x_continuous(breaks = 2019:2023) +
    labs(
      x = "Year",
      y = "Cumulative Unique N-gram Types",
      color  = "Individual"          # only color appears in legend
    ) +
    theme_bw() +
    theme(
      plot.title  = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
}

# --- Generate plots ---
g1 <- plot_ngram_panel("Gram1")
g2 <- plot_ngram_panel("Gram2") +
  theme(legend.position = "none", axis.title.y = element_blank())
g3 <- plot_ngram_panel("Gram3") +
  theme(legend.position = "none", axis.title.y = element_blank())


#--------------------------------------------------------------------------
library(grid)
library(patchwork)

# --- Headers for NLD row ---
n1 <- wrap_elements(grid::textGrob("Within Day", gp = gpar(fontsize = 12, fontface = "bold")))
n2 <- wrap_elements(grid::textGrob("Across Days", gp = gpar(fontsize = 12, fontface = "bold")))
n3 <- wrap_elements(grid::textGrob("Across Years", gp = gpar(fontsize = 12, fontface = "bold")))
nld_header <- n1 + n2 + n3

# --- Headers for N-gram row ---
h1 <- wrap_elements(grid::textGrob("1-Gram", gp = gpar(fontsize = 12, fontface = "bold")))
h2 <- wrap_elements(grid::textGrob("2-Gram", gp = gpar(fontsize = 12, fontface = "bold")))
h3 <- wrap_elements(grid::textGrob("3-Gram", gp = gpar(fontsize = 12, fontface = "bold")))
ngram_header <- h1 + h2 + h3

# --- Row labels ---
label_A <- wrap_elements(grid::textGrob("(A)", x = unit(0, "npc"), hjust = 0, gp = gpar(fontsize = 14, fontface = "bold")))
label_B <- wrap_elements(grid::textGrob("(B)", x = unit(0, "npc"), hjust = 0, gp = gpar(fontsize = 14, fontface = "bold")))

# --- Row plots ---
nld_row   <- p1 + p2 + p3
ngram_row <- g1 + g2 + g3

# --- Spacer to increase row gap ---
gap <- plot_spacer() + plot_spacer() + plot_spacer()




#_______________


library(tidyverse)
library(vegan)

# Use consistent colors and shapes
individual_order <- c("YYXX", "YOXX", "RUGG", "YOGG")
colors <- setNames(RColorBrewer::brewer.pal(4, "Set2"), individual_order)
shapes <- setNames(c(16, 17, 15, 3), individual_order)

# Filter 2022 data
df_2022 <- df %>%
  filter(Year == "2022", Individual %in% individual_order)

# Function to get accumulation data for one individual and one N-gram level
get_accum_curve <- function(ind, gram_col) {
  df_ind <- df_2022 %>%
    filter(Individual == ind, !is.na(.data[[gram_col]])) %>%
    group_by(Song_Num) %>%
    summarise(across(all_of(gram_col), unique), .groups = "drop") %>%
    unnest(cols = all_of(gram_col)) %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = all_of(gram_col), values_from = value, values_fill = 0) %>%
    arrange(Song_Num) %>%
    select(-Song_Num)
  
  if (nrow(df_ind) < 2 || ncol(df_ind) < 1) return(NULL)
  
  acc <- specaccum(df_ind, method = "random", permutations = 100)
  tibble(
    Songs = acc$sites,
    Richness = acc$richness,
    SD = acc$sd,
    Individual = ind,
    Ngram = gram_col
  )
}

# Get accumulation data for all combinations
accum_all <- map_dfr(c("Gram1", "Gram2", "Gram3"), function(gram) {
  map_dfr(individual_order, function(ind) get_accum_curve(ind, gram))
}) %>%
  mutate(Ngram = recode(Ngram,
                        "Gram1" = "1-Gram",
                        "Gram2" = "2-Gram",
                        "Gram3" = "3-Gram"))

# --- Define accumulation curve plotter ---
plot_accum_panel <- function(ngram_label) {
  df_plot <- accum_all %>% filter(Ngram == ngram_label)
  
  ggplot(df_plot, aes(x = Songs, y = Richness, color = Individual, shape = Individual)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2, alpha = 0.6) +
    geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD, fill = Individual), alpha = 0.2, color = NA) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    labs(
      x = "Number of Songs",
      y = "Accumulated Unique N-gram Types",
      color = "Individual"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = ""
    )
}

b1 <- plot_accum_panel("1-Gram")
b2 <- plot_accum_panel("2-Gram") +
  theme(legend.position = "none", axis.title.y = element_blank())
b3 <- plot_accum_panel("3-Gram") +
  theme(legend.position = "none", axis.title.y = element_blank())


#----


# Update headers
label_B <- wrap_elements(grid::textGrob("(B)", x = unit(0, "npc"), hjust = 0, gp = gpar(fontsize = 14, fontface = "bold")))
label_C <- wrap_elements(grid::textGrob("(C)", x = unit(0, "npc"), hjust = 0, gp = gpar(fontsize = 14, fontface = "bold")))

accum_header <- h1 + h2 + h3  # reuse from earlier

accum_row <- b1 + b2 + b3

# Update final layout
final_plot <- (
  label_A / nld_header / nld_row /
    gap /
    label_B / accum_header / accum_row /
    gap /
    label_C / ngram_header / ngram_row
) +
  plot_layout(
    heights = c(0.1, 0.1, 1, 0.01, 0.1, 0.1, 1, 0.01, 0.1, 0.1, 1),
    guides = "collect"
  ) +
  plot_annotation(
    theme = theme(
      legend.position = "top",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.margin = margin(t = 20, b = 10, l = 10, r = 10)
    )
  )


final_plot

