#
# Kiem Lab
# Fred Hutchinson Cancer Center
# November 2024
#
# 2. Calculate the background error rate for the pre-treatment samples.
# Compare the pre-treatment samples with the edited samples from the same animal to find potential SNPs. 
# Identify the SNPs and annotate them.
#
# USE THIS SCRIPT ONLY ON PRE-TREATMENT SAMPLES.
#

baseDir = 'File path to results for each timepoint.'
setwd(baseDir)

sample1 = 'Pre-treatment sample'
sample2 = "Choose a later timepoint from the same animal to identify SNPs."

mock <- 'Read the protospacer base error rate file for the pre-treatment sample.'
edited <- 'Read the protospacer base error rate file for the edited/later timepoint for the same animal.'

# Make a dataframe with the filtered snps.
df.all <- as.data.frame(matrix(nrow = 0, ncol = 21))
colnames(df.all) <- c("Lab",c(1:20))

df.snp <- as.data.frame(matrix(nrow = 0, ncol = 8))
colnames(df.snp) <- c("Lab", "Position", "Base_change", "Frequency", "A", "C", "G", "T")

for(i in seq(1, length(mock$Lab), by = 20)){
 
  # Retrieve all the 20 positions for the off-target.
  mock.freq <- mock[c(i:(i+19)),]
  
  # Check if the off-target is present in the edited sample so that it can be a probable SNP.
  if(mock$Lab[i] %in% edited$Lab == "TRUE"){
    # Identify the index of the off-target in the edited sample.
    edit.bp <- (which(edited$Lab == mock$Lab[i]))[1]
    edit.freq <- edited[c(edit.bp:(edit.bp+19)),]
    
    # Assert that the positions are in the order 1 to 20.
    if(all(edit.freq$Position) == all(c(1:20))){
      # Filter the SNP if the frequency is greater than 10% (0.1).
      snp <- which(edit.freq$Frequency >= 0.1)
      
      if(length(snp) > 0){
        # Check if the frequency of the SNP is also greater than 10% in the pre-treatment.
        t <- as.numeric(which(mock.freq$Frequency[c(snp)] >= 0.1))
        print(mock$Lab[i])
        print(mock.freq$Frequency[c(snp)])
        print(t)
        
        if(length(t) > 0){
          # Check if the base changes are also the same at these positions.
          if(edit.freq$Base_change[snp] == mock.freq$Base_change[snp[t]]){
            
            # Save the snp locations and the base changes in the mock.
            Lab <- rep(mock$Lab[i], length(t))
            s.tmp <- cbind.data.frame(Lab, mock.freq$Position[snp[t]], mock.freq$Base_change[snp[t]], 
                     mock.freq$Frequency[snp[t]], mock.freq$A[snp[t]], mock.freq$C[snp[t]], 
                     mock.freq$G[snp[t]], mock.freq$T[snp[t]])
            colnames(s.tmp) <- c("Lab","Position","Base_change","Frequency","A","C","G","T")
            df.snp <- rbind.data.frame(df.snp, s.tmp)
            
            # Include the off-target in the analysis but don't include the SNP bases.
            mock.freq$Frequency[snp[t]] <- "snp"
            
            tmp <- cbind.data.frame(mock$Lab[i],t(mock.freq$Frequency))
            colnames(tmp) <- c("Lab",c(1:20))
            df.all <- rbind.data.frame(df.all, tmp)
          }
          else{
            tmp4 <- cbind.data.frame(mock$Lab[i],t(mock.freq$Frequency))
            colnames(tmp4) <- c("Lab",c(1:20))
            print(tmp4)
            df.all <- rbind.data.frame(df.all, tmp4)
          }
        }
        else{
          tmp1 <- cbind.data.frame(mock$Lab[i],t(mock.freq$Frequency))
          colnames(tmp1) <- c("Lab",c(1:20))
          df.all <- rbind.data.frame(df.all, tmp1)
        }
        
      }
      else{
        tmp2 <- cbind.data.frame(mock$Lab[i],t(mock.freq$Frequency))
        colnames(tmp2) <- c("Lab",c(1:20))
        df.all <- rbind.data.frame(df.all, tmp2)
      }
    }
    else{
      print("Position error")
    }
  }
  else{
    tmp3 <- cbind.data.frame(mock$Lab[i],t(mock.freq$Frequency))
    colnames(tmp3) <- c("Lab",c(1:20))
    df.all <- rbind.data.frame(df.all, tmp3)
  }
}

# Save the SNPs and their frequencies for each animal.
write.csv(df.snp, paste0(baseDir,sample1,"_snps_freq.csv"), row.names = F)

# Calculate the mean for each protospacer position by excluding the SNP postions (na).
cmean <- c()
csd <- c()
df.stat <- df.all[,2:21]

for(i in 1:20){
  na <- as.numeric(which(df.stat[,i] == "snp"))
  
  if(length(na) > 0){
    print(df.stat[na,i])
    var <- as.numeric(df.stat[-c(na),i])
    assign(paste("p", i, sep = "_"), var)  
    pmean <- mean(var)
    cmean <- c(cmean, pmean)
    psd <- sd(var)
    csd <- c(csd, psd)
  }
  else{
    var <- as.numeric(df.stat[,i])
    assign(paste("p", i, sep = "_"), var)  
    pmean <- mean(var)
    cmean <- c(cmean, pmean)
    psd <- sd(var)
    csd <- c(csd, psd)
  }
}

# Calculate the cumulative mean and the standard deviation for each animal based on the pre-treatment.
mean(cmean)
mean(csd)



