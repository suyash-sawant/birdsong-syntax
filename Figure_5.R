# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(patchwork)
library(lubridate)

setwd("C:/Users/suyas/Desktop/WBS songs")

# Constants
focus_inds <- c("YYXX", "YOXX", "RUGG", "YOGG")
df <- read.csv("Final_Classified_Notes_Dec16.csv") %>%
  mutate(
    Individual = gsub("__", "XX", Individual),
    Year = as.character(Year),
    Note_DateTime = dmy_hms(paste(Note_Date, Note_Time), tz = "UTC")
  )

# Function for plotting unique N-grams (individual signatures)
get_signature_plot <- function(gram_col, title_label) {
  df_sub <- df %>%
    filter(Year %in% c("2021", "2022"), !is.na(.data[[gram_col]])) %>%
    select(Individual, Year, Ngram_Type = all_of(gram_col))
  
  unique_ngrams <- df_sub %>%
    filter(!is.na(Individual)) %>%
    distinct(Individual, Year, Ngram_Type) %>%
    group_by(Year, Ngram_Type) %>%
    filter(n() == 1) %>%
    ungroup()
  
  unique_counts <- df_sub %>%
    inner_join(unique_ngrams, by = c("Individual", "Year", "Ngram_Type")) %>%
    count(Individual, Year, Ngram_Type, name = "Count") %>%
    group_by(Individual, Ngram_Type) %>%
    mutate(Total = sum(Count)) %>%
    ungroup()
  
  top_unique <- unique_counts %>%
    group_by(Individual) %>%
    slice_max(Total, n = 3, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      Ngram_Type = fct_reorder(Ngram_Type, Total, .desc = TRUE),
      Individual = factor(Individual, levels = focus_inds)
    ) %>% 
    filter(!is.na(Individual))
  
  # Add dummy rows for empty individuals
  missing_inds <- setdiff(focus_inds, unique(top_unique$Individual))
  if (length(missing_inds) > 0) {
    dummy <- expand.grid(Individual = missing_inds, Year = c("2021", "2022"), Ngram_Type = "(none)")
    top_unique <- bind_rows(top_unique, dummy %>% mutate(Count = 0, Total = 0))
  }
  
  ggplot(top_unique, aes(x = Ngram_Type, y = Count, fill = Year)) +
    geom_col(position = "stack") +
    facet_wrap(~Individual, nrow = 1, scales = "free_x") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.6))) +  # More space on top
    scale_fill_manual(values = c("2021" = "#2b8cbe", "2022" = "#045a8d")) +
    labs(title = title_label, x = "Unique N-gram Type", y = "Occurrences", fill = "Year") +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "top"
    )
}

# Function for plotting shared N-grams (population dialect)
get_shared_plot <- function(gram_col) {
  df_sub <- df %>%
    filter(Year %in% c("2021", "2022"), !is.na(.data[[gram_col]])) %>%
    select(Individual, Year, Ngram_Type = all_of(gram_col)) %>%
    mutate(Year = as.factor(Year))
  
  # Get shared N-grams (present in all individuals that year)
  inds_per_year <- df_sub %>% distinct(Year, Individual) %>% count(Year, name = "n_inds")
  
  shared_ngrams <- df_sub %>%
    distinct(Year, Individual, Ngram_Type) %>%
    count(Year, Ngram_Type, name = "n_present") %>%
    left_join(inds_per_year, by = "Year") %>%
    filter(n_present == n_inds) %>%
    select(Year, Ngram_Type)
  
  # Count total usage and get top 3 per year
  top_shared <- df_sub %>%
    semi_join(shared_ngrams, by = c("Year", "Ngram_Type")) %>%
    count(Year, Ngram_Type, name = "Total") %>%
    group_by(Year) %>%
    slice_max(Total, n = 3, with_ties = FALSE) %>%
    ungroup()
  
  # Count usage for those N-grams in focal individuals
  shared_counts <- df_sub %>%
    filter(Individual %in% focus_inds) %>%
    inner_join(top_shared, by = c("Year", "Ngram_Type")) %>%
    count(Individual, Year, Ngram_Type, name = "Count") %>%
    group_by(Ngram_Type) %>%
    mutate(Total = sum(Count)) %>%
    ungroup() %>%
    mutate(
      Individual = factor(Individual, levels = focus_inds),
      Ngram_Type = fct_reorder(Ngram_Type, Total, .desc = TRUE)
    )
  
  # Plot
  ggplot(shared_counts, aes(x = Individual, y = Count, fill = Year)) +
    geom_col(position = "stack") +
    facet_wrap(~Ngram_Type, nrow = 1, scales = "free_x") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.6))) +
    scale_fill_manual(values = c("2021" = "#2b8cbe", "2022" = "#045a8d")) +
    labs(
      title = "Population Dialect: Top 3 Shared Gram1 Types per Year",
      x = "Shared N-gram Type",
      y = "Occurrences",
      fill = "Year"
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )
}


# Generate plots
sig1 <- get_signature_plot("Gram1", "Individual Signatures: Top Unique 1-grams")
sig2 <- get_signature_plot("Gram2", "Individual Signatures: Top Unique 2-grams")
sig3 <- get_signature_plot("Gram3", "Individual Signatures: Top Unique 3-grams")
shared1 <- get_shared_plot("Gram1")
# Combine plots
(sig1 / sig2 / sig3 / shared1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

#_____________


setwd("C:/Users/suyas/Desktop/WBS songs/Singnatures")

list.files()


library(ggplot2)
library(ggimage)
library(tuneR)
library(seewave)
library(dplyr)
library(ggimage)

wav_files <- list.files("C:/Users/suyas/Desktop/WBS songs/Singnatures",
                        pattern = "\\.wav$", full.names = TRUE)

# Create a directory for images
dir.create("Spectrograms", showWarnings = FALSE)


for (file in wav_files) {
  # 1. Read in
  wav <- readWave(file)
  sr  <- wav@samp.rate
  orig  <- wav@left
  n_orig <- length(orig)
  
  # 2. Pad to 1.5 seconds if needed
  n_target <- sr * 1.5
  if (n_orig < n_target) {
    pad_total  <- n_target - n_orig
    pad_left   <- floor(pad_total / 2)
    pad_right  <- pad_total - pad_left
    new_left   <- c(
      rep(0, pad_left),
      orig,
      rep(0, pad_right)
    )
    wav <- Wave(left = new_left, samp.rate = sr, bit = wav@bit)
  }
  
  # 3. Plot spectrogram with x and y scales, no color scale
  png_filename <- file.path("Spectrograms",
                            paste0(tools::file_path_sans_ext(basename(file)), ".png"))
  png(png_filename, width = 1200, height = 1200, res = 300)
  
  par(mar = c(4, 4, 1, 1))  # Room for axis labels
  
  spectro(
    wav,
    f     = sr,
    flim  = c(1, 10),     # Y-axis: 1 to 10 kHz
    tlim  = c(0, 1.5),    # X-axis: 0 to 1.5 sec
    scale = FALSE,        # Disable color scale
    axisX = TRUE,         # Show X axis (time)
    axisY = TRUE,         # Show Y axis (frequency)
    osc   = FALSE,
    grid  = FALSE
  )
  
  dev.off()
}
