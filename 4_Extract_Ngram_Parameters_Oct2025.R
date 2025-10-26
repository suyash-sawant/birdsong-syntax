library(dplyr)
library(tidyr)
library(purrr)
library(zoo)

df <- read.csv("Final_Classified_Notes_Dec16.csv") %>%
  filter(!paste(Individual, Year, sep = "-") %in% c(
    "BWGG-2022", "BWRU-2021", "GG__-2022", "RU__-2019", 
    "RW__-2019", "YOUU-2021", 
    "YO__-2019", "YY__-2019","GG__-2019"))

# First, add Mean_Freq_Hz column
df <- df %>%
  mutate(
    Mean_Freq_Hz = (Low_Freq_Hz + High_Freq_Hz) / 2
  )

# ---------- Gram1 traits ----------
gram1_traits <- df %>%
  group_by(Ngram_ID = Gram1) %>%
  summarise(
    Mean_Position = mean(Note_Num, na.rm = TRUE),
    Mean_Freq = mean(Mean_Freq_Hz, na.rm = TRUE),
    Delta_Freq = mean(Delta_Freq_Hz, na.rm = TRUE),
    Duration = mean(Delta_Time_s, na.rm = TRUE),
    Ngram_Level = "Gram1",
    .groups = "drop"
  )

# ---------- Function to compute Gram2 or Gram3 traits ----------
get_ngrams <- function(df, n = 2, trim = 0.1) {
  df %>%
    arrange(Individual, Year, Song_Num, Note_Num) %>%
    group_by(Individual, Year, Song_Num) %>%
    mutate(
      Ngram_ID = zoo::rollapply(Gram1, width = n, FUN = function(x) paste(x, collapse = ""), fill = NA, align = "left"),
      Pos = zoo::rollapply(Note_Num, width = n, FUN = mean, fill = NA, align = "left"),
      Begin_s = zoo::rollapply(Begin_Time_s, width = n, FUN = min, fill = NA, align = "left"),
      End_s = zoo::rollapply(End_Time_s, width = n, FUN = max, fill = NA, align = "left"),
      Low_Min = zoo::rollapply(Low_Freq_Hz, width = n, FUN = min, fill = NA, align = "left"),
      High_Max = zoo::rollapply(High_Freq_Hz, width = n, FUN = max, fill = NA, align = "left")
    ) %>%
    ungroup() %>%
    filter(!is.na(Ngram_ID)) %>%
    mutate(
      Delta_F = High_Max - Low_Min,
      Mean_F = (Low_Min + High_Max) / 2,
      Dur = End_s - Begin_s
    ) %>%
    group_by(Ngram_ID) %>%
    summarise(
      Mean_Position = mean(Pos, trim = trim, na.rm = TRUE),
      Mean_Freq = mean(Mean_F, trim = trim, na.rm = TRUE),
      Delta_Freq = mean(Delta_F, trim = trim, na.rm = TRUE),
      Duration = mean(Dur, trim = trim, na.rm = TRUE),
      Ngram_Level = paste0("Gram", n),
      .groups = "drop"
    )
}


# ---------- Apply to Gram2 and Gram3 ----------
gram2_traits <- get_ngrams(df, n = 2)
gram3_traits <- get_ngrams(df, n = 3)

# ---------- Combine all ----------
ngram_traits_all <- bind_rows(gram1_traits, gram2_traits, gram3_traits)

# View result
head(ngram_traits_all)
write.csv(ngram_traits_all, "Ngram_Parameters_All.csv", row.names = FALSE)
