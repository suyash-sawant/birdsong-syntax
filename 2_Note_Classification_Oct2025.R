rm(list = ls())
# Packages required
library(readr)
library(progress)
library(Rraven)
library(data.table)
library(warbleR) # Ensure this is loaded for cross_correlation


# Set dir
#setwd("D:/CA/GRASS-notes")
Dir <- "C:/Users/suyas/OneDrive/Desktop/WBS_songs/Set_07/manually_classified_data_Individual-Notes"
setwd(Dir)

# Training and test data file names
train_data <- "Subsampled_Training_Data.txt"
test_data <- "manually_classified_data_Individual-Notes.txt"

# Read Selection Tables
train_data_warbler <- Rraven::imp_raven(path = Dir, warbler.format = TRUE,
                                        files = train_data)
test_data_warbler <- Rraven::imp_raven(path = Dir, warbler.format = TRUE,
                                       files = test_data)

file_path_test <- file.path(Dir, test_data)
file_path_train <- file.path(Dir, train_data)

train_data_complete <- fread(file_path_train)
#train_data_complete$`Note Type` <- NA
train_data_complete$`Correlation` <- NA
test_data_complete <- fread(file_path_test)
test_data_complete$`Note types` <- NA
test_data_complete$`Correlation` <- NA

common_columns <- intersect(colnames(train_data_complete), colnames(test_data_complete))
test_data_complete <- test_data_complete[, ..common_columns]


# Note Types
train_types <- train_data_complete$`Note types`
unique_train_types <- unique(train_types)

# Create combinations of possible note type names - A01, A02, A03,...
combinations <- expand.grid(LETTERS, sprintf("%02d", 1:30))
combinations_vec <- paste(combinations$Var1, combinations$Var2, sep = "")

unique_new_types <- setdiff(combinations_vec, unique_train_types)
rm(combinations_vec, combinations)

# Set Parameters
lowlim <- 0.8 # Threshold for new note type
highlim <- 0.97 # Threshold to stop look and assign note type

Test_N <- nrow(test_data_warbler[c(1:50),]) # Test for a subset

# Progress bar
pb <- progress::progress_bar$new(
  format = "  [:bar] :current/:total (:percent) elapsed: :elapsed eta: :eta",
  total = Test_N, clear = FALSE, width = 60
)

# Function to save data regularly
save_data_csv <- function(test_data, iteration) {
  fwrite(test_data, file = sprintf("test_data_complete_%04d.csv", iteration))
}

# FOR loop
for (Test_n in 1:Test_N) {
  
  corr_vec <- c()
  Train_N <- nrow(train_data_warbler)
  
  #Randomization of rows in tarining data
  train_data_warbler <- train_data_warbler[sample(nrow(train_data_warbler)),]
  train_data_complete <- 
    train_data_complete[match(train_data_warbler$selec,
                              train_data_complete$Selection),]
  
  for (Train_n in 1:Train_N) {
    correlation <- suppressMessages(suppressWarnings(
      warbleR::cross_correlation(
        X = rbind(train_data_warbler[Train_n,],
                  test_data_warbler[Test_n,]), pb = FALSE,
        wl = 512, ovlp = 70, bp = c(1, 10))))
    
    corr <- correlation[2,1]
    
    if (corr >= highlim) {
      corr_vec <- c(corr_vec, corr)
      
      break
    }
    corr_vec <- c(corr_vec, corr)
  }
  
  if (max(corr_vec) >= lowlim) {
    test_data_complete$`Note types`[Test_n] <- 
      train_data_complete$`Note types`[which(corr_vec == max(corr_vec))]
    test_data_complete$Correlation[Test_n] <- max(corr_vec)
    print(sum(train_data_complete$`Note types` == 
                train_data_complete$`Note types`[which(corr_vec == max(corr_vec))]))
    if (sum(train_data_complete$`Note types` == 
            train_data_complete$`Note types`[which(corr_vec == max(corr_vec))]) 
        < 5){  #add replicate of note type in the training data if less than 5 replicates
      
      train_data_warbler <- rbind(train_data_warbler, 
                                  test_data_warbler[Test_n,])
      train_data_warbler$selec[nrow(train_data_warbler)] <- 
        nrow(train_data_warbler)
      train_data_complete <- 
        rbind(train_data_complete, 
              test_data_complete[Test_n,])
      train_data_complete$Selection[nrow(train_data_complete)] <- 
        nrow(train_data_complete)
    }
    
  } else {
    test_data_complete$`Note types`[Test_n] <- unique_new_types[1]
    test_data_complete$Correlation[Test_n] <- 1
    
    train_data_warbler <- rbind(train_data_warbler, test_data_warbler[Test_n,])
    train_data_warbler$selec[nrow(train_data_warbler)] <- 
      nrow(train_data_warbler)
    train_data_complete <- 
      rbind(train_data_complete, 
            test_data_complete[Test_n, 1:ncol(test_data_complete)])
    train_data_complete$Selection[nrow(train_data_complete)] <- 
      nrow(train_data_complete)
    
    unique_new_types <- unique_new_types[-1]
  }
  
  pb$tick()
  
  # Save data every 100 iterations
  if (Test_n %% 100 == 0) {
    save_data_csv(test_data_complete, Test_n)
  }
}

# Final save
save_data_csv(test_data_complete, Test_N)

#write.table(train_data_complete, "Subsampled_Training_Data.txt", sep = "\t", quote = F,row.names = FALSE)
