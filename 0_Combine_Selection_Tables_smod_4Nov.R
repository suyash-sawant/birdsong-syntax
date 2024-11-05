File_Dur <- function(Directory, No_Digit=3) {
  List_wav <- list.files(path = Directory, pattern = ".*.WAV")
  File_count <- as.numeric(length(List_wav))
  File_lengths <-list()
  for ( i in 1:File_count){
    File_name <- List_wav[i]
    File_name <- base::gsub("WAV", "wav", File_name)
    wav <- tuneR::readWave(paste(Directory,File_name, sep = "/"))
    sound_length <- base::round(length(wav@left) / wav@samp.rate, No_Digit)
    File_lengths[i] <- sound_length
  }
  return(File_lengths)
}

Combine_Tables <- function(Directory, File_Dur = c()) {
  
  notefilename <- list.files(path = Directory, pattern = ".*.txt")
  Number_of_files <- length(notefilename)
  
  File_Dur <- as.array(File_Dur)
  
  #reading first file in the list
  notefile1 <- notefilename[1]
  Note_Selections1 <-
    as.data.frame(readr::read_tsv(notefile1, col_names = T, show_col_types = FALSE))
  Note_Selections1$`Song No.` <- as.numeric(Note_Selections1$`Song No.`)
  
  Time_diff1 <- 0
  
  #Clubbing all the files into one raven importable table
  for (x in c(2:(Number_of_files))){
    
    notefile2 <- notefilename[x]
    Previous_songs <- length(unique(Note_Selections1$`Song No.`))
    Time_diff2 <- Time_diff1 + as.numeric(File_Dur[x-1])
    
    Note_Selections2 <-
      as.data.frame(readr::read_tsv(notefile2, col_names = T, show_col_types = FALSE))
    Note_Selections2$`Song No.` <-
      as.numeric(Note_Selections2$`Song No.`) + Previous_songs
    Note_Selections2$`Begin Time (s)` <-
      Note_Selections2$`Begin Time (s)` + Time_diff2
    Note_Selections2$`End Time (s)` <-
      Note_Selections2$`End Time (s)` + Time_diff2
    
    clubfile <- rbind(Note_Selections1,Note_Selections2)
    
    Note_Selections1 <- clubfile
    Time_diff1 <- Time_diff2
  }
  
  #set the selection numbers
  clubfile[,1] <- c(1:nrow(clubfile))
  
  #convert song number format from '1' to '001'
  clubfile$`Song No.` <- as.factor(clubfile$`Song No.`)
  clubfile$`Song No.` <- formatC(clubfile$`Song No.`,x,flag="0",width=3)
  
  return(clubfile)
}

Directory <- "C:/Users/suyas/OneDrive/Desktop/WBS_songs/Sholicola songs selection tables/BWXX-2021"
File.Dur <- File_Dur(Directory, No_Digit=3)
Combine.Table <- Combine_Tables(Directory, File_Dur = File.Dur)

write.table(Combine.Table, "Combine-notes-BW__-2021.txt",
            sep = "\t", quote = F,row.names = FALSE)
