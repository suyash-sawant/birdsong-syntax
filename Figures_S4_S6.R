library(dplyr)
library(ggplot2)
library(factoextra)
library(tidyr)
library(patchwork)

# Set working directory
setwd("C:/Users/suyas/OneDrive/Desktop/WBS_songs")

# Load and prepare data
df <- read.csv("WBS_Final_Data - Songs_Params.csv", header = TRUE)
df$Individual <- as.factor(df$Individual)
df$Year <- as.factor(df$Year)
df$Month <- as.factor(df$Month)

# Custom color palettes
years_colors <- c(
  "2019" = "#a6bddb", "2021" = "#2b8cbe",
  "2022" = "#045a8d", "2023" = "#023858"
)

group1_colors <- c(
  "BWXX" = "#fee8c8", "RWXX" = "#fdd49e", "RUXX" = "#fdbb84",
  "GGXX" = "#fc8d59", "YYXX" = "#ef6548", "BWRW" = "#d7301f",
  "RUGG" = "#b30000", "BWGG" = "#7f0000", "YOXX" = "#993404", "YORW" = "#662506"
)

group2_colors <- c(
  "UUXX" = "#e5f5e0", "YOUU" = "#c7e9c0", "YOGG" = "#a1d99b",
  "BWRU" = "#74c476", "BWGY" = "#41ab5d", "BWRR" = "#238b45", "PWGG" = "#006d2c"
)

individual_colors <- c(group1_colors, group2_colors)

# Normalize names like "YY__" → "YYXX"
df$Individual <- gsub("__", "XX", df$Individual)

# PCA variable groups
pca_vars1 <- c("Low_Freq_Hz", "High_Freq_Hz", "Delta_Freq_Hz", "Mean_Freq_Hz", "Delta_Time_s")
pca_vars2 <- c("Note_Pace", "Note_Count", "Note_Types", "Note_Richness", "Note_Diversity",
               "Note_Order")

run_pca_plot_grid <- function(df, vars, label) {
  # Prepare data
  pca_df <- df %>%
    select(Individual, Year, all_of(vars)) %>%
    drop_na() %>%
    mutate(across(all_of(vars), scale))
  
  # Show year distribution to check for missing 2023
  message("Included years for ", label, ":")
  print(table(pca_df$Year))
  
  # PCA
  pca_result <- prcomp(pca_df %>% select(all_of(vars)), center = TRUE, scale. = TRUE)
  scores_df <- cbind(pca_df[, c("Individual", "Year")], pca_result$x)
  
  # Filter key individuals
  key_inds <- c("YYXX", "YOXX", "YOGG", "RUGG")
  facet_df <- scores_df %>% filter(Individual %in% key_inds)
  
  # Plot 1: Color by Individual, facet by Year
  p1 <- ggplot(scores_df, aes(x = PC1, y = PC2, color = Individual)) +
    geom_point(alpha = 0.7) +
    facet_wrap(~Year, ncol = 4) +
    scale_color_manual(values = individual_colors) +
    labs(title = NULL) +
    theme_bw()+
    theme(
      legend.position = "top",
      legend.box = "horizontal"
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))
  
  # Plot 2: Color by Year, facet by Individual
  p2 <- ggplot(facet_df, aes(x = PC1, y = PC2, color = Year)) +
    geom_point(alpha = 0.7) +
    facet_wrap(~Individual, ncol = 4) +
    scale_color_manual(values = years_colors) +
    labs(title = NULL) +
    theme_bw()+
    theme(
      legend.position = "top",
      legend.box = "horizontal"
    )
  
  
  # Plot 3a: Scree plot
  explained_var <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  explained_df <- data.frame(
    PC = paste0("PC", 1:length(explained_var)),
    Variance_Explained = explained_var
  )
  
  p3a <- ggplot(explained_df, aes(x = PC, y = Variance_Explained)) +
    geom_col(fill = "grey50") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = "Variance Explained", y = "Variance Explained") +
    theme_bw()
  
  # Plot 3b: Arrow biplot (using factoextra)
  p3b <- fviz_pca_biplot(
    pca_result,
    repel = TRUE,
    col.var = "black",
    col.ind = "gray60",
    label = "var"
  ) +
    theme_bw() +
    labs(title = "PCA Loadings")
  
  # Combine and return
  return(p1 / p2 / (p3a | p3b) + plot_layout(ncol = 1, heights = c(1, 1, 1))+
           plot_annotation(title = label))
}


# Run for both PCA sets
run_pca_plot_grid(df, pca_vars1, "Spectro-Temporal Parameters")
run_pca_plot_grid(df, pca_vars2, "Song Complexity Measurements")
run_pca_plot_grid(df, c(pca_vars1, pca_vars2),
                  "Spectro-Temporal Parameters & Song Complexity Measurements")
