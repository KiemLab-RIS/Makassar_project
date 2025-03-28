#
# Kiem Lab
# Fred Hutchinson Cancer Center
# November 2024
#
# 3. For each animal, identify the SNPs from the pre-treatment sample in all the edited time points.
#
library(ggplot2)

baseDir = 'File path to results for each timepoint.'
setwd(baseDir)

edited <- 'Read the protospacer base error rate file for the edited/later timepoint.'
snp <- 'Read the corresponding SNP file for the animal.'

# Identify substitutions in the edited samples that are not SNPs and are greater than the background error rate.
df.edit_snp <- as.data.frame(matrix(nrow = 0, ncol = 4))
colnames(df.edit_snp) <- c("File","Position","Frequency","Base_change")

# Remove the off-targets with frequency more than 10% (0.1)
#
for (i in seq(1, length(edited$Lab), by = 20)){
  
  edit.freq <- edited[c(i:(i+19)),]
  # Check if the edited off-target is present in the snp file.
  if(edited$Lab[i] %in% snp$Lab == "TRUE"){
    
    # Find the index of the off-target in the snp file.
    snp.bp <- which(snp$Lab == edited$Lab[i])
    snp.pos <- snp$Position[snp.bp]
    print(edited$Lab[i])
    
    # Check if the base error rate is greater than background error rate (0.004) in the edited sample.
    poi <- which(edit.freq$Frequency > 0.004)
    print(poi)
    
    if(length(poi) > 0){
      for(j in 1:length(poi)){
        # Check if the bp with frequency > 0.004 is a snp.
        if(poi[j] %in% snp.pos == "FALSE"){
          
          tmp <- cbind.data.frame(edited$Lab[i], poi[j], (edit.freq$Frequency[poi[j]]), edit.freq$Base_change[poi[j]])
          colnames(tmp) <- c("File","Position","Frequency","Base_change")
          df.edit_snp <- rbind.data.frame(df.edit_snp, tmp)
        }
      }
    }
  }
  else{
    # Check if the base error rate is greater than background error rate (0.004) in the edited sample.
    poi <- which(edit.freq$Frequency > 0.004)
    
    if(length(poi) > 0){
      for(j in 1:length(poi)){
        
        tmp <- cbind.data.frame(edited$Lab[i], poi[j], (edit.freq$Frequency[poi[j]]), edit.freq$Base_change[poi[j]])
        colnames(tmp) <- c("File","Position","Frequency","Base_change")
        df.edit_snp <- rbind.data.frame(df.edit_snp, tmp)
      }
    }
  }
}

# Save all the substitutions in all the samples that are not SNPs and are above the background error rate.
write.csv(df.edit_snp, paste0(baseDir,sample,"_mutations_after_filtering.csv"), row.names = F)
