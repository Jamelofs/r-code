#calling the necessary packages:
suppressWarnings({library(ggplot2)
library(stringr)
library(data.table)
library(readr)
library(dplyr)
library(tidyr)
library(reshape)     
library(gplots)
library(pheatmap)
library(heatmap3)
library(magrittr) 
library(tibble) 
library(RColorBrewer)
library(paletteer)
library(clr)
library(tidyverse) #used for data science. The eight core packages inside this library are: ggplot2 (data visualisation), dplyr (data manipulation), tidyr, readr, purrr, tibble, stringr, and forcats
library(KODAMA) # to use the normalisation function
library(ggrepel) #mainly used to repel overlapping text labels in ggplots
library(vegan) #popular library for analysing ecological diversity and for multivariate analysis of community data. Here, we use it for PCoA
library(svglite) #to save the plots in support vector graphics (svg) format
library(factoextra) #for extracting and visualizing outputs of multivariate analyses such as PCA, k-means
library(ggsci) #provides color palettes for ggplot2 that can be used for scientific journals
library(matrixStats) #contains highly optimized functions to perform statistics on matrix data
library(cowplot)}) #efficient functions to arrange several plots
####################################Data clean-up steps###########################################

#upload feature table and metadata
setwd("\\\\SERVER\\Users\\james\\Desktop\\DESKTOP\\NAUSeurm")
ft <- read.csv("mzmine_pos_test1_iimn_gnps_quant_no_wash.csv", header = T, check.names = F, sep = ",")
md <- read.csv("Cope_NAU_Metadata_serum_pos - Copy.tsv", header = T, check.names = F, sep = "\t")

#Function: InsideLevels - checking data
InsideLevels <- function(metatable){
  lev <- c()
  typ<-c()
  for(i in 1:ncol(metatable)){
    x <- levels(droplevels(as.factor(metatable[,i])))
    if(is.double(metatable[,i])==T){x=round(as.double(x),2)}
    x <-toString(x)
    lev <- rbind(lev,x)
    
    y <- class(metatable[,i])
    typ <- rbind(typ,y)
  }
  out <- data.frame(INDEX=c(1:ncol(metatable)),ATTRIBUTES=colnames(metatable),LEVELS=lev,TYPE=typ,row.names=NULL)
  return(out)
}

InsideLevels(md)
new_ft <- ft
new_md <- md

#Removing Peak area extensions
colnames(new_ft) <- gsub(' Peak area','',colnames(new_ft))
rownames(new_md) <- gsub(' Peak area','',rownames(new_md))

new_md <- new_md[,colSums(is.na(new_md))<nrow(new_md)] #Removing if any NA columns present in the md file
rownames(new_md) <- trimws(rownames(new_md), which = c("both")) #remove the (front & tail) spaces, if any present, from the rownames of md

#Removing the spaces (if any) from each column of md and converting them all to UPPERCASE
for(i in 1:ncol(new_md)){
  if(is.factor(new_md[,i]) | is.character(new_md[,i]) == T){
    new_md[,i] <- trimws(new_md[,i], which = c("both")) #First remove spaces in the front and end of each column of md
    new_md[,i] <- gsub(' ','_', new_md[,i]) # Replace the spaces (in the middle) to underscore
    #new_md[,i] <- factor(casefold(new_md[,i], upper=T)) #convert all to UPPERCASE
  } else if (is.numeric(new_md[,i]) | is.integer(new_md[,i]) | is.double(new_md[,i]) == T){
    new_md[,i] <- new_md[,i]
  }
}

#Changing the row names of the files
#new_ft <- data.frame(new_ft)
rownames(new_ft) <- paste(new_ft$'row ID',round(new_ft$'row m/z',digits = 3),round(new_ft$'row retention time',digits = 3), sep = '_')

#Picking only the files with column names containing 'mzML'
new_ft <- new_ft[,grep('mzML',colnames(new_ft))]

#check again
InsideLevels(new_md)
new_ft<- new_ft[,order(colnames(new_ft))] #ordering the ft by its column names
new_md <-new_md[order(rownames(new_md)),] #ordering the md by its row names
write.csv(new_md, "~/NAU_pos_new_md.csv",row.names = TRUE)

new_md <- read_csv('~/NAU_pos_new_md.csv') %>% column_to_rownames("filename")
#new_md <- new_md[, -1]

#lists the colnames(ft) which are not present in md
{unmatched_ft <- colnames(new_ft)[which(is.na(match(colnames(new_ft),rownames(new_md))))] 
cat("These", length(unmatched_ft),"columns of feature table are not present in metadata:")
if((length(unmatched_ft) %% 2) ==0)
{matrix(data=unmatched_ft,nrow=length(unmatched_ft)/2,ncol=2)}else
{matrix(data=unmatched_ft,nrow=(length(unmatched_ft)+1)/2,ncol=2)}

flush.console()
Sys.sleep(0.2)

#lists the rownames of md which are not present in ft
unmatched_md <- rownames(new_md)[which(is.na(match(rownames(new_md),colnames(new_ft))))] 
cat("These", length(unmatched_md),"rows of metadata are not present in feature table:")
if((length(unmatched_md) %% 2) ==0)
{matrix(data=unmatched_md,nrow=length(unmatched_md)/2,ncol=2)}else
{matrix(data=unmatched_md,nrow=(length(unmatched_md)+1)/2,ncol=2)}}

#checking the dimensions of our new ft and md:
cat("The number of rows and columns in our original ft is:",dim(ft),"\n")
cat("The number of rows and columns in our new ft is:",dim(new_ft),"\n")
cat("The number of rows and columns in our new md is:",dim(new_md))

#If subset_data exists, it will take it as "data", else take new_md as "data"
if(exists("subset_data")==T){data <-subset_data}else{data <-new_md}
InsideLevels(data)

Condition <- as.double(unlist(readline("Enter the index number of the attribute to split sample and blank:")))
Levels_Cdtn <- levels(droplevels(as.factor(data[,Condition[1]])))

###############################Selecting Blanks#########################################
#Among the shown levels of an attribute, select the ones to keep
Blk_id <- as.double(unlist(readline("Enter the index number of your BLANK:")))
paste0('You chosen blank is:',Levels_Cdtn[Blk_id])

#Splitting the data into blanks and samples based on the metadata
md_Blank <- data[(data[,Condition] == Levels_Cdtn[Blk_id]),]
Blank <- new_ft[,which(colnames(new_ft)%in%rownames(md_Blank)),drop=F] 
md_Samples <- data[(data[,Condition] != Levels_Cdtn[Blk_id]),]
Samples <- new_ft[,which(colnames(new_ft)%in%rownames(md_Samples)),drop=F] 

#write.csv(new_md, 'md_new_2.csv', row.names=TRUE)

###################################Step 1: Blank removal###################################
if(casefold(readline('Do you want to perform Blank Removal- Y/N:'),upper=T)=='Y'){
  
  #When cutoff is low, more noise (or background) detected; With higher cutoff, less background detected, thus more features observed
  Cutoff <- as.numeric(readline('Enter Cutoff value between 0.1 & 1:')) # (i.e. 10% - 100%). Ideal cutoff range: 0.1-0.3
  
  #Getting mean for every feature in blank and Samples
  Avg_blank <- rowMeans(Blank, na.rm= FALSE, dims = 1) # set na.rm = FALSE to check if there are NA values. When set as TRUE, NA values are changed to 0
  Avg_samples <- rowMeans(Samples, na.rm= FALSE, dims = 1)
  
  #Getting the ratio of blank vs Sample
  Ratio_blank_Sample <- (Avg_blank+1)/(Avg_samples+1)
  
  # Creating a bin with 1s when the ratio>Cutoff, else put 0s
  Bg_bin <- ifelse(Ratio_blank_Sample > Cutoff, 1, 0 )
  Blank_removal <- cbind(Samples,Bg_bin)
  
  # Checking if there are any NA values present. Having NA values in the 4 variables will affect the final dataset to be created
  temp_NA_Count <-cbind(Avg_blank ,Avg_samples,Ratio_blank_Sample,Bg_bin)
  
  print('No of NA values in the following columns:')
  print(colSums(is.na(temp_NA_Count)))
  
  #Calculating the number of background features and features present
  print(paste("No.of Background or noise features:",sum(Bg_bin ==1,na.rm = TRUE)))
  print(paste("No.of features after excluding noise:",(nrow(Samples) - sum(Bg_bin ==1,na.rm = TRUE)))) 
  
  Blank_removal <- Blank_removal %>% filter(Bg_bin == 0) # Taking only the feature signals
  Blank_removal <- as.matrix(Blank_removal[,-ncol(Blank_removal)]) # removing the last column Bg_bin 
}

#check output - blank removed
dim(Blank_removal)

write.csv(Blank_removal, "~\\Blank_removed.csv", row.names =TRUE)

#########################################Step 2 - LOD Imputation#########################
#creating bins from -1 to 10^10 using sequence function seq()
bins <- c(-1,0,(1 * 10^(seq(0,10,1)))) 

#cut function cuts the give table into its appropriate bins
scores_gapfilled <- cut(as.matrix(Blank_removal),bins,labels = c('0','1','10','1E2','1E3','1E4','1E5','1E6','1E7','1E8','1E9','1E10')) 

#transform function convert the tables into a column format: easy for visualization 
FreqTable<-transform(table(scores_gapfilled)) #contains 2 columns: "scores_x1", "Freq"
FreqTable$Log_Freq <- log(FreqTable$Freq+1) #Log scaling the frequency values

colnames(FreqTable)[1] <- 'Range_Bins'
#FreqTable #Uncomment the line if you want to see the FreqTable used for the following ggplot.

## GGPLOT2
ggplot(FreqTable, aes(Range_Bins, Log_Freq))+ 
  geom_bar(stat="identity",position = "dodge", width=0.3) + 
  scale_fill_brewer(palette = "Set1") +
  ggtitle(label="Frequency plot - Gap Filled") +
  xlab("Range") + ylab("(Log)Frequency") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +   # setting the angle for the x label
  theme(axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1)) +   # setting the angle for the y label
  theme(plot.title = element_text(hjust = 0.5))

Cutoff_LOD <- round(min(Blank_removal[Blank_removal!=min(Blank_removal)]))
print(paste0("The minimum value greater than 0 in gap-filled table: ",Cutoff_LOD))

##############################################################
#random number imputation
#set.seed(141222)
#ran_values <-round(runif(length(Blank_removal),0,Cutoff_LOD),digits=1)
#ran_values
#blk_rem <- data.frame(Blank_removal)
#imp <- blk_rem  %>% mutate(across(everything(), ~replace(., . == 0 , sample(ran_values, size=1))))
#if(casefold(readline('Do you want to perform Imputation? - Y/N:'),upper=T) == 'Y'){
  #Imputed <- blk_rem  %>%
    #mutate(blk_rem, ~replace(., . == 0 , sample(ran_values, size=1)))
  #head(Imputed)
#}
################################################################
#minimum LOD imputation
if(casefold(readline('Do you want to perform Imputation? - Y/N:'),upper=T) == 'Y'){
  Imputed <- Blank_removal
  Imputed[Imputed <Cutoff_LOD] <- Cutoff_LOD
  head(Imputed)
}

Ydim(Imputed)

write.csv(Imputed,"Imputed_QuantTable_filled_with_116_20240124.csv",row.names =TRUE)

###############################Step 3 - Normalization - do if step 4 is not performed########################
#library(KODAMA) # to use the normalization function

Normalized_data_TIC <- t(normalization(t(Imputed), method = "sum")$newXtrain) 
head(Normalized_data_TIC,n=3)
dim(Normalized_data_TIC)
print(paste('No.of NA values in Normalized data:',sum(is.na(Normalized_data_TIC)== T)))
write.csv(Normalized_data_TIC, "Normalised_Quant_table.csv",row.names =TRUE)

#################################Step 4 - scaling with clr transform function - do if step 3 is not performed##############
Imp_t <- as.data.frame(t(Imputed)) #transposing the imputed table
# center and scale data with scale() function
library(vegan)
Imp_clrt <- Imp_t %>% decostand(method = "rclr")
write.csv(Imp_clrt, "IMP_clr.csv",row.names = TRUE)


############################UNIVARIATE#################

