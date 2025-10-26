
Extract_Syllables <- function(Raven_Selections, max.ngram.length = 3,
                              SongNo, NoteType, Notes = c()){
  
  ##---------------------------------------------------------------------------------------------------------------------------
  
  #function for inserting spaces in a character string
  fun_insert <- function(x, pos, insert) {
    base::gsub(paste0("^(.{", pos, "})(.*)$"),
               paste0("\\1", insert, "\\2"),
               x)
  }
  
  ##---------------------------------------------------------------------------------------------------------------------------
  
  #Import the song sequences to be analyzed
  Raven_Selections <- Raven_Selections[Notes,c(1:NoteType)]
  Raven_Selections <- stats::na.omit(Raven_Selections)
  
  max_n <- max.ngram.length
  
  ##---------------------------------------------------------------------------------------------------------------------------
  
  Classified_notes <-
    as.data.frame(Raven_Selections[,NoteType], stringsAsFactors=FALSE)
  colnames(Classified_notes) <- c("Note Types")
  
  nt_len <- as.numeric(nchar(Classified_notes[1,1]))
  
  Total_notes <- nrow(Raven_Selections)
  
  #Extract length of each song in terms of number of notes
  Note_Count <-
    stats::aggregate(Raven_Selections$Channel ,
                     by = list(Category = Raven_Selections[,SongNo]), FUN = sum )
  colnames(Note_Count) <- c("song_no","note_count")
  
  ##---------------------------------------------------------------------------------------------------------------------------
  
  #Create datasheet for the start and end of each song
  start_end <- data.frame(matrix(nrow = nrow(Note_Count), ncol = 2))
  start_end <- cbind(Note_Count$note_count,start_end)
  colnames(start_end) <- c("Note_Count","Start_note", "end_note")
  for (l in 2:nrow(start_end)){
    start_end[1,2] <- 1
    start_end[1,3] <- as.numeric(start_end[1,1])
    start_end[l,2] <- 1 + as.numeric(start_end[l-1,3])
    start_end[l,3] <- as.numeric(start_end[l,1]) + as.numeric(start_end[l-1,3])
  }
  rm(l)
  
  ##---------------------------------------------------------------------------------------------------------------------------
  
  Song_seq <- data.frame(matrix(nrow = nrow(Note_Count)))
  
  
  for (m in 1:nrow(start_end)){
    x = start_end[m,2]
    y = start_end[m,3]
    z = start_end[m,1]
    
    Song_seq[m,] <- t(as.data.frame(do.call(paste, quote = F, as.list(Classified_notes[c(x:y),1]))))
    Song_seq[m,] <- gsub(" ","",Song_seq[m,])
    colnames(Song_seq) <- c("Note_Types_in_song")
  }
  rm(m,x,y,z, Note_Count)
  
  ##---------------------------------------------------------------------------------------------------------------------------
  
  Old_New <- data.frame(matrix(nrow = Total_notes, ncol = 2))
  colnames(Old_New) <- c("Old","New")
  
  #Fill the data frame for transition of note types
  for (l in 1:Total_notes){
    Old_New[l,1] <- Classified_notes[l,1]
    Old_New[l,2] <- Classified_notes[l+1,1]
  }
  
  Old_New <- Old_New[-start_end$end_note,]
  
  #Extract transition probabilities for each transition
  transTable <- base::prop.table(with(Old_New, table(Old,New)), 1)
  
  Trans_prob_long <- as.data.frame(as.table(transTable))
  Trans_prob_long$Freq[is.nan(Trans_prob_long$Freq)] <- 0
  
  rm(Old_New,l,transTable,Total_notes,Classified_notes,start_end)
  
  ##---------------------------------------------------------------------------------------------------------------------------
  
  #Create data frame for the output
  N_grams <- data.frame(matrix(ncol=2, nrow=nrow(Song_seq)))
  colnames(N_grams) <- c("Song","N_grams")
  
  ##---------------------------------------------------------------------------------------------------------------------------
  
  for (i1 in 1:nrow(Song_seq)){
    
    Song_string <- as.character(Song_seq[i1,1])
    song_length <- length((strsplit(Song_string, "")[[1]]))/nt_len
    song_gaps <- (length((strsplit(Song_string, "")[[1]]))/nt_len)-1
    
    gap_combinations <- unlist(Map(utils::combn, list(c(1:song_gaps)),
                                   seq_along(c(1:song_gaps)),
                                   simplify = FALSE),recursive = FALSE)
    
    waysofwriting = data.frame(matrix())
    
    
    for (i2 in 1:length(gap_combinations)){
      gap_loc <- unlist(gap_combinations[i2])
      
      Song_string1 <- Song_string
      for ( i3 in 1:length(gap_loc)){
        
        Song_string2 <- fun_insert(x = Song_string1,
                                   pos = nt_len*(gap_loc[i3])+i3-1,
                                   insert = " ")
        Song_string1 <- Song_string2
        
      }
      waysofwriting = rbind(waysofwriting, Song_string2)
    }
    waysofwriting[1,] <- Song_string
    
    waysofwriting1 <-
      waysofwriting %>%
      tidyr::separate(matrix.., sep = " ",
                      into = as.character(c(1:song_length)), fill = "right")
    
    waysofwriting1[is.na(waysofwriting1)] <- 0
    
    waysofwriting2 <- data.frame(matrix(ncol = ncol(waysofwriting1), nrow = nrow(waysofwriting1)))
    
    for (j2 in 1:nrow(waysofwriting1)){
      for (j3 in 1:ncol(waysofwriting1)){
        if (nchar(waysofwriting1[j2,j3]) <= max_n*nt_len){
          waysofwriting2[j2,j3] <- waysofwriting1[j2,j3]
        } else {
          waysofwriting2[j2,j3] <- NA
        }
      }
    }
    waysofwriting2 <- stats::na.omit(waysofwriting2)
    waysofwriting2[waysofwriting2 == 0] <- NA
    
    waysofwriting <- data.frame(matrix(ncol =1, nrow = nrow(waysofwriting2)))
    
    for(i in 1:nrow(waysofwriting2)){
      waysofwriting[i,1]  <- paste0(waysofwriting2[i,1:sum(!is.na(waysofwriting2[i,]))], collapse = " ")
    }
    
    intra_tran_prob = data.frame(matrix(nrow=nrow(waysofwriting2),
                                        ncol=ncol(waysofwriting2)))
    
    for (i4 in 1:nrow(waysofwriting2)){
      for (i5 in 1:ncol(waysofwriting2)){
        prob0A = 0
        way <- waysofwriting2[i4,i5]
        if ( length(way) != 0){
          way1 <- unlist(strsplit(way, split = ""))
          way2 <- character(length = length(way1)/nt_len)
          for (j1 in 0:(length(way1)/nt_len)-1) {
            str1 <- paste(way1[nt_len*j1+1],way1[nt_len*j1+2],way1[nt_len*j1+3], sep="")
            way2[j1+1] <- str1
          }
          way2 <- as.character(way2)
        }
        
        #print(way2)
        if (length(way2)==1){
          prob0A = 0
        }else{
          for (i6 in 2:length(way2)){
            note1A<- as.character(way2[i6-1])
            note2A<- as.character(way2[i6])
            prob1A <- Trans_prob_long$Freq[which(Trans_prob_long$Old == note1A &
                                                   Trans_prob_long$New == note2A)]
            if ( is.numeric(prob1A) == TRUE){
              prob1A <- prob1A
            } else{
              prob1A = 0
            }
            prob0A = (prob0A + prob1A)/(length(way2)-1)
          }
          
        }
        #print(paste(way, prob0A, i6, i5, sep = ","))
        if ( length(prob0A) == 0){
          intra_tran_prob[i4,i5] <- 0
        } else{
          intra_tran_prob[i4,i5] <- prob0A
        }
        
      }
    }
    
    
    intra_prob = data.frame(matrix(nrow=nrow(intra_tran_prob)))
    for (i7 in 1:nrow(intra_tran_prob)){
      intra_prob[i7,] <- sum(intra_tran_prob[i7,])/sum(intra_tran_prob[i7,]!=0)
    }
    intra_prob[nrow(intra_prob),] <- 0
    DeltaP <- cbind(waysofwriting, intra_prob)
    
    inter_tran_prob = data.frame(matrix(ncol=1, nrow=nrow(waysofwriting2)))
    
    for (i8 in 1:nrow(waysofwriting2)){
      prob0B = 0
      len1 <- sum(!is.na(waysofwriting2[i8,]))
      for (i9 in 2:len1){
        
        string1 <- waysofwriting2[i8,i9-1]
        string2 <- waysofwriting2[i8,i9]
        
        note1B<- as.character(substr(string1, nchar(string1)-nt_len+1, nchar(string1)))
        note2B<- as.character(substr(string2, 1, nt_len))
        
        prob1B <- Trans_prob_long$Freq[which(Trans_prob_long$Old == note1B &
                                               Trans_prob_long$New == note2B)]
        
        if ( is.numeric(prob1B) == TRUE){
          prob1B = prob1B
        } else{
          prob1B = 0
        }
        prob0B = prob0B + prob1B
        #print(prob0B)
      }
      if ( length(prob0B) == 0){
        inter_tran_prob[i8,] <- 0
      } else{
        inter_tran_prob[i8,] <- prob0B/len1
      }
      #print(prob0B)
    }
    
    DeltaP[,3]<-inter_tran_prob
    DeltaP[,4]<-DeltaP[,2]-DeltaP[,3]
    colnames(DeltaP)<- c("Way", "Intra_Prob(P)", "Inter_Prob(p)", "Delta_P")
    
    Bestway <- DeltaP$Way[which(DeltaP$Delta_P == max(DeltaP$Delta_P))]
    N_grams[i1,1]<- as.character(Song_string)
    N_grams[i1,2]<- as.character(Bestway[1])
    
    rm(DeltaP,gap_combinations,inter_tran_prob,intra_prob,intra_tran_prob,
       waysofwriting, waysofwriting1,waysofwriting2,Bestway,gap_loc,
       i,i1,i2,i3,i4,i5,i6,i7,i8,i9,j1,j2,j3,len1,note1A,note1B,note2A,note2B,
       prob0A,prob0B,prob1A,prob1B,song_gaps,song_length,str1,
       Song_string,Song_string1,Song_string2,string1,string2,way,way1,way2)
  }
  ##___________________________________________________________________________________________________________________________
  
  return(N_grams)
}
