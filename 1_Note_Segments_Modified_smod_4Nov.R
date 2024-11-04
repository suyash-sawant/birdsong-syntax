################################################################################
## this script is for extracting individual note level WAV files from a recordings
## this script will only run for one WAV recordings file and its corresponding 
## selection table with note level annotations


#JUST ONE SELECTION TABLE - MANUAL - FILE.CHOOSE(READ.TABLE())
#install.packages("Rraven")
#install.packages("tuneR")
#install.packages("seewave")

library(Rraven)
library(seewave)
library(tuneR)

rm(list = ls())
setwd("C:/Users/suyas/OneDrive/Desktop/WBS_songs/Sholicola songs selection tables/BWRR-2022") 

dir()


##below is the function to extract note segments
Note_segment <- function(File, Input.Path, Output.Path, tbuffer=0.05, fbuffer=500) ##set tbuffer (sec) and fbuffer(Hz) and keep it consistent across all recordings
  {
  
  if (is.null(Input.Path)){
    Input.Path = getwd()
  }
  setwd(Input.Path)
  
  if (is.null(Output.Path)){
    print("Output Path not provided")
    stop()
  }
  
  if (!dir.exists(Output.Path)) {
    dir.create(Output.Path, recursive = TRUE)
  }
  
  # Inspect the data
  #names(Selection_table)
  
  Selection_table$`Delta Time (s)` <- Selection_table$`End Time (s)` - Selection_table$`Begin Time (s)`
  
  #Create Sound_files vector
  read_sound_files <- Selection_table$`Begin File`
  Sound_files <- file.path(Input.Path, read_sound_files)
  

  Start_time <- Selection_table$`File Offset (s)`
  End_time <- Selection_table$`File Offset (s)` + Selection_table$`Delta Time (s)`
  Total_notes <- nrow(Selection_table)
  
  pb <- txtProgressBar(min = 0, max = Total_notes, style = 3) #progress bar
  
  #print(Total_notes)
  output_file_names <- c()
  for (i in 1:Total_notes) {
    
    input_sound_file_name <- Sound_files[i]  
    
    #print(input_sound_file_name)  # Debug: Check file path
    
    # Debug: Check if the file exists before trying to read it
    if (!file.exists(input_sound_file_name)) {
      stop(paste("Audio file not found:", input_sound_file_name))
    }
    
    # Use tuneR::readWave() to read the file
    input_sound_file <- tuneR::readWave(
      filename = input_sound_file_name, 
      from = Start_time[i] - tbuffer, 
      to = End_time[i] + tbuffer, 
      units = "seconds"
    )
    
    # Filtering Sound
    low.freq <- Selection_table$`Low Freq (Hz)`[i] - fbuffer
    low.freq <- base::pmax(low.freq, 0)
    high.freq <- Selection_table$`High Freq (Hz)`[i] + fbuffer
    #print(paste("Low frequency:", low.freq))   # Debug: Check low frequency
    #print(paste("High frequency:", high.freq)) # Debug: Check high frequency
    
    output_sound_file <- seewave::ffilter(
      input_sound_file,
      f = input_sound_file@samp.rate,
      from = low.freq,
      to = high.freq,
      output = "Wave",
      rescale = TRUE
    )
    
    # Optional: Visualize Spectrogram
    #par(mar = c(2, 7, 7, 2))
    #seewave::spectro(output_sound_file, f = 44100)
    
    output_file_name <- paste(
      paste(tools::file_path_sans_ext(basename(input_sound_file_name)),  # Remove extension
      formatC(Selection_table$Selection[i], width = 4, flag = "0"),  ##add note number extension with 4 char to the filename while saving a note
      sep = "_"),
    ".WAV", sep=""
    )# Add the ".WAV" extension

    # Save the filtered sound file
    seewave::savewav(output_sound_file,
                     f = input_sound_file@samp.rate,
                     channel = 1,
                     filename = paste(Output.Path, output_file_name, sep = "/"))
    
    output_file_names <- c(output_file_names,output_file_name)
    
    setTxtProgressBar(pb, i)
  }
  
  Output_selections <- Selection_table
  
  for (i in 1:Total_notes){
    Output_selections$`End Time (s)`[i] <-
      sum(Output_selections$`Delta Time (s)`[c(1:i)]) +
      tbuffer +(2*tbuffer*(i-1))
    Output_selections$`Begin Time (s)`[i] <-
      Output_selections$`End Time (s)`[i]-
      Output_selections$`Delta Time (s)`[i]
    Output_selections$`File Offset (s)`[i] <-
      tbuffer
    
    input_sound_file_name <- Sound_files[i]
    output_file_name <-
      paste(paste(tools::file_path_sans_ext(input_sound_file_name),
                  Selection_table$`Song No.`[i],
                  Selection_table$Selection[i], sep = "_"),
            ".WAV", sep = "")
  }
  Output_selections$`Begin File` <- output_file_names
  
  Output_sel_table <- paste(tools::file_path_sans_ext(File),
                            "Note-segments",sep="-")
  write.table(Output_selections,
              file = paste(Output.Path,
                           paste(Output_sel_table,"txt",sep="."),
                           sep = "/"),
              sep = "\t", quote = F,row.names = FALSE)
}


Selection_table <- as.data.frame(readr::read_tsv(file.choose(),col_names = T,
                                                 show_col_types = FALSE)) #select the selection table in the pop up generated
#Calling the function
Note_segment(File="Combine-notes-BWRR-2022.txt",#name of the file to be segmented
             Input.Path="C:/Users/suyas/OneDrive/Desktop/WBS_songs/Sholicola songs selection tables/BWRR-2022",#path to the file to be segmented
             Output.Path="C:/Users/suyas/OneDrive/Desktop/WBS_songs/Sholicola songs segments/BWRR-2022", #save to wd, this is not incorporated in the script yet
             tbuffer=0.05, #time buffer set for each note
             fbuffer=500) #frequency buffer set for each note
  
