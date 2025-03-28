#
# Kiem Lab
# Fred Hutchinson Cancer Center
# November 2024
# 1. Script to calculate the background error rate in the protospacer sequence in the pre-treatment samples.
# 
#
# Load the required packages.
library('stringr')
library('stringi')
library('ggplot2')
library('plotly')
library('pheatmap')
library('dplyr')
library('GMCM')

#A21051 <- c("A21051_pre_treatment","A21051_12DPT","A21051_26DPT","A21051_40DPT","A21051_54DPT","A21051_68DPT","A21051_83DPT","A21051_118DPT")
#A21033 <- c("A21033_13DPT","A21033_26DPT","A21033_40DPT","A21033_54DPT")
#A20134 <- c("A20134_pre_treatment","A20134_26DPT","A20134_33DPT","A20134_47DPT")

sample <- "Sample to be run from the above 3."
testDir = 'File path to results for each timepoint.'
setwd(testDir)
dir.create('File path to create the results_directory.')
resdir = 'File path to the results_directory.'

# Filtering the files based on read count for each off-target. 
# We chose 1000 as the cut-off based on the read count curve (log scale).
# The tsv files contain the alignment results of the fastq sequences with the reference and their frequencies.
# Create a dataframe with file names and labels.
fs <- list.files(testDir, pattern = "*.tsv", full.names = T, recursive = F)
fs <- basename(fs)
fs
ts = data.frame(fs,stringsAsFactors = FALSE)
ts 
colnames(ts) = c("file")

lt = lapply(ts[,1],function(s) {
  l = substr(s,1,nchar(s)-4)
})
ts$label = unlist(lt)
ts

# Create a dataframe to read each file and save the filename, the strand to which the sequence aligned and the off-target seqeunce. 
df.filter <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(df.filter) <- c("label","Strand","OfSeq")

for (i in 1:length(ts$file)){
  print(ts$file[i])
  
  if (ts$file[i] == "out.tsv" | ts$file[i] == "proto_spacer_mean.tsv"){
    next
  }
  else{
    
    ff <- read.table(ts$file[i], sep = "\t", header = F)
    # Retrieve the filename and the strand.
    label <- ts$label[i]
    Strand <- ff[3,1]
    # Also add the off-target sequence to the file.
    OfSeq <- str_split((str_trim(ff[5,2])), " ")[[1]][1]
    
    tmp <- cbind.data.frame(label, Strand, OfSeq)
    df.filter <- rbind.data.frame(df.filter, tmp)
  }
}
# Save this dataframe with the file name, strand alignment and off-target sequence. 
write.csv(df.filter, paste0(resdir,sample,"_df_filter.csv"), row.names = F)

# If df.filter already exists, read the saved file.
df.filter <- read.csv(paste0(resdir,sample,"_df_filter.csv"), stringsAsFactors = T, header = T)

# Create a dataframe with filenames and labels for the results files.
# These files contain the counts for all bases at each position in the sequence.
# Eg: position 20: A-10, C-0, G-180, T-0
files <- list.files(testDir, pattern = "^r_*", full.names = T,recursive = F)
files <- basename(files)
files
tests = data.frame(files,stringsAsFactors = FALSE)
tests 
colnames(tests) = c("file")

lt = lapply(tests[,1],function(s) {
  l = substr(s,3,nchar(s)-4)
})
tests$label = unlist(lt)
tests

# Merge the above 2 dataframes, so that each filename is now associated with the base counts, strand and off-target sequence.
df.all <- merge.data.frame(tests, df.filter, by = "label")

# Routine to process each file and calculate for each position in the sequence the following:
# Original base, base change if any (substitutions) and the counts for all bases at that position.
df.guide <- data.frame(data.frame(matrix(nrow = 0,ncol = 9)))
colnames(df.guide) <- c("Lab", "Position", "Frequency", "Base_change", "Off_target", "A", "C", "G", "T")

for (i in 1:length(df.all$file)){

  # Read the first line of the file to check for the total sequence count.
  f <- readLines(df.all$file[i])[1]
  Count <- as.numeric(strsplit(f,",")[[1]][2])
  
  # Read the rest of the file to get the base information.
  fh <- read.csv(df.all$file[i], stringsAsFactors = T, header = T, comment.char = "#")
  colnames(fh) <- c("Base", "Protospacer_loc", "Ref", "A", "C", "G", "T", "Gap")
  
  # Retrieve the position where the protospacer sequence starts to align.
  guide_pos <- as.numeric(fh$Protospacer_loc[1])
  
  # Check if the total sequence count meets the cut-off criteria.
  if(Count >= 1000){
    
    # Make a new data frame to include just the protospacer base pair frequencies.
    # Find the index which matches the guide position.
    gp <- which(fh$Base == guide_pos)
    
    # Depending on which strand the off-target aligns to, the numbering of the positions will change.
    if (df.all$Strand[i] == "pro.f"){
      
      # Make a small data frame with the protospacer positions. [ACGT]
      # The length of the protospacer sequence is 20.
      pro.df <- fh[c(gp:(gp + 19)),]
      Frequency <- c()
      Base_change <- c()
      # Retrieve the base counts.
      ct.df <- pro.df[,c(4:7)]
      # We will also need to calculate the base errors and the base change.
      # Base error = mismatches/total sequences.
      for (j in 1:length(pro.df$Base)){
        # Reference base.   
        ref <- pro.df$Ref[j]
        # All the base counts.
        bp.freq <- pro.df[j,4:7]
        # Which of the bases matches the reference, exclude it from calculating the base error rate.
        b.col <- which(colnames(bp.freq) == ref)
        b.error <- (sum(bp.freq[-c(b.col)]))/Count
        Frequency <- c(Frequency, round(b.error,5))
        
        # Look at the base change.
        if(b.error > 0){
          bc <- bp.freq[-c(b.col)]
          m <- as.numeric(which.max(bc))
          be_change <- paste0(ref,"-",colnames(bc[m]))
          Base_change <- c(Base_change, be_change)
        }
        else{
          be_change <- paste0("NoChange")
          Base_change <- c(Base_change, be_change)
        }
      }
      # Add the protospacer frequencies and base change to df.guide dataframe.
      Lab <- rep(df.all$label[i],20)
      Position <- c(1:20)
      Off_target <- rep(df.all$OfSeq[i],20)
      tmp1 <- data.frame(Lab, Position, Frequency, Base_change, Off_target)
      colnames(tmp1) <- c("Lab", "Position", "Frequency", "Base_change", "Off_target")
      tmp1 <- cbind.data.frame(tmp1, ct.df)
      df.guide <- rbind.data.frame(df.guide, tmp1)
      
    }
    else{
      # If the off-target aligns to the reverse strand. fh[c((gp + 3):(gp + 22)),]
      print(df.all$file[i])
      pro.df <- fh[c(gp:(gp + 19)),]
      Frequency <- c()
      Base_change <- c()
      ct.df <- rev(pro.df[,c(4:7)])
      colnames(ct.df) <- c("A", "C", "G", "T")
      # Since the sequence aligns to the reverse strand, we need to look at the reverse complement.
      rev_ct.df <- cbind.data.frame(rev(ct.df$A), rev(ct.df$C), rev(ct.df$G), rev(ct.df$T))
      colnames(rev_ct.df) <- c("A", "C", "G", "T")
      
      for (j in 1:length(pro.df$Base)){
        
        ref <- pro.df$Ref[j]
        bp.freq <- pro.df[j,4:7]
        b.col <- which(colnames(bp.freq) == ref)
        b.error <- (sum(bp.freq[-c(b.col)]))/Count
        Frequency <- c(Frequency, round(b.error,5))
        
        # Look at the base change.
        if(b.error > 0){
          bc <- bp.freq[-c(b.col)]
          m <- as.numeric(which.max(bc))
          
          #Reverse complement the reference base since it is on the reverse strand.
          if(ref == "A"){
            ref.rev = "T"
          }
          else if(ref == "T"){
            ref.rev = "A"
          }
          else if(ref == "G"){
            ref.rev = "C"
          }
          else if(ref == "C"){
            ref.rev = "G"
          }
          
          # Reverse complement the changed base.
          mut <- colnames(bc[m])
          if(mut == "A"){
            mut.rev = "T"
          }
          else if(mut == "T"){
            mut.rev = "A"
          }
          else if(mut == "G"){
            mut.rev = "C"
          }
          else if(mut == "C"){
            mut.rev = "G"
          }
          
          be_change <- paste0(ref.rev,"-",mut.rev)
          Base_change <- c(Base_change, be_change)
        }
        else{
          be_change <- paste0("NoChange")
          Base_change <- c(Base_change, be_change)
        }
      }
      
      # Add the protospacer freqeuncies and base change to df.guide dataframe.
      Lab <- rep(df.all$label[i],20)
      Position <- c(1:20)
      Off_target <- rep(df.all$OfSeq[i],20)
      tmp1 <- data.frame(Lab, Position, rev(Frequency), rev(Base_change), Off_target)
      colnames(tmp1) <- c("Lab", "Position", "Frequency", "Base_change", "Off_target")
      tmp1 <- cbind.data.frame(tmp1, rev_ct.df)
      df.guide <- rbind.data.frame(df.guide, tmp1)
    }
  }
  else {
    next
  }
}

# The result files contains the per base error rate for each of the 20 bases in the protospacer sequence.
# We have 900 different off-targets for each timepoint.
write.csv(df.guide, paste0(resdir,sample,"_protospacer_error_rates_reverse_comp.csv"), row.names = F)


