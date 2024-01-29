#read in libraries
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(dplyr)
library(stringr)
# Read the data into R
setwd("\\\\SERVER\\Users\\james\\Desktop\\DESKTOP\\NAUFecal\\")


isotope = "16O"
week = "24"



fileNameNQT = paste(isotope,"W",week, sep = "")
NQTFileName = "Normalised_Quant_table.csv"
fileNameMetadata = "Cope_NAU_FecalMetadata_pos_no_wash.tsv"
metadata <- read_tsv(fileNameMetadata,num_threads =1,progress = show_progress())
dset<-read.csv("Normalised_Quant_table.csv", header = TRUE)


# metadata <- metadata[metadata$ATTRIBUTE_Isotope %in% c(isotope), ]
metadata <- metadata[metadata$ATTRIBUTE_Week %in% c(week), ]

wtMetadata <- metadata[metadata$ATTRIBUTE_Strain %in% c("WT"), ]

wtFileNames <- c(wtMetadata$filename)
wtDSet <- dset[, wtFileNames, drop = FALSE]
wtData<- as.matrix(wtDSet %>% select(-1))



tgMetadata <- metadata[metadata$ATTRIBUTE_Strain %in% c("Tg"), ]

tgFileNames <- c(tgMetadata$filename)
tgDSet <- dset[, tgFileNames, drop = FALSE]
tgData<- as.matrix(tgDSet %>% select(-1))

data <- as.matrix(cbind(wtDSet,tgDSet))
featureIDs <- dset[,1]

# Check the loaded dataset
# Dimension of the dataset
dim(data) 
head(data)

# Separate the two conditions into two smaller data frames
wt = wtData
ko = tgData

# Compute the means of the samples of each condition
wt.mean = apply(wt, 1, mean)
ko.mean = apply(ko, 1, mean)

# Just get the maximum of all the means for plotting
limit = max(wt.mean, ko.mean)
means <- cbind.data.frame(wt.mean, ko.mean)

###############################################################
######### WT VS KO ANALYSIS FROM THIS POINT ONWARD ###########
DEres.df <- as.data.frame(data)
DEres.df <- cbind.data.frame(featureIDs, DEres.df)
#compute fold
fold = wt.mean / ko.mean
log2fold <- log2(fold)
#add to results dataframe
DEres.df$logFC <- log2fold

# Generate histogram of the fold differences
foldhist <- ggplot(DEres.df, aes(x=logFC)) +
  geom_histogram( binwidth=0.01, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("Fold Change Distribution") +
  labs(x= "Log fold change in abundance between groups", y = "Count") +
  theme_minimal() 
foldhist

#compute stat sig - t-test
pvalue = NULL # Empty list for the p-values
tstat = NULL # Empty list of the t test statistics
for(i in 1 : nrow(data)) { # For each gene : 
  x = wt[i,] # WT of gene number i
  y = ko[i,] # KO of gene number i
  
  # Compute t-test between the two conditions
  t = t.test(x, y, na.rm=TRUE, paired=FALSE, exact=FALSE)
  
  # Put the current p-value in the pvalues list
  pvalue[i] = t$p.value
  # Put the current t-statistic in the tstats list
  tstat[i] = t$statistic
}
# add results to results dataframe
DEres.df$pvalue <- pvalue

# Generate -log10(pvalues)
pvallogged <- -log10(pvalue)
# add results to results dataframe
DEres.df$logpval <- pvallogged
# Generate histogram of -log10(pvalues)
pval_hist <- ggplot(DEres.df, aes(x=logpval)) +
  geom_histogram( binwidth=0.05, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
  ggtitle("-log10(pvalue) Distribution") +
  labs(x= "Significance levels occuring in D.E. results", y = "Count") +
  xlim(0,4) +
  theme_minimal() 
pval_hist
# Volcano Plots:
# put the biological significance (fold changes)
# and statistical significance (p-value) in one plot
# Generate the volcano plot
volcanoPlot <- function(results.df, fold_cutoff=0.5,
                        pvalue_cutoff=0.01, title = 
                          "Volcano Plot of Features"){
  # Inputs: ##############################################
  
  # results.df: dataframe containing columns: 
  #             1) "featureIDs": feature IDs
  #             2) "logFC": foldchange values from DE analysis
  #             3) "logpval": log-transformed p-values from DE analysis
  # fold cutoff: significance threshold to render visually on plot;
  #               denotes fold difference between mutant and wildtype
  #               also referred to as "biologial signal"
  # p_value_cutoff: significance threshold to render visually on plot;
  #                 denotes statistical significance
  ########################################################
  # create factor denoting differential expression status
  stats.df <- results.df %>% mutate(diffexpressed = 
                                      ifelse(logFC < -(fold_cutoff) & 
                                               logpval >= -log10(pvalue_cutoff),
                                             "Significant, downregulated", 
                                             ifelse(logFC > fold_cutoff & 
                                                      logpval >= -log10(pvalue_cutoff), 
                                                    "Significant, upregulated", "Normally expressed")))
  
  stats.df$diffexpressed <- as.factor(stats.df$diffexpressed)
  DEres.df <- stats.df %>%
    mutate(featureIDs = sub("_.*", "", featureIDs)) %>% filter(!(diffexpressed == "Normally expressed"))
  write.csv(DEres.df, "chr_output.csv", row.names = FALSE)
  
  # Base plot
  plot <- ggplot(stats.df, aes(x = logFC, y = logpval, col = diffexpressed, text = featureIDs)) +
    geom_point() +
    ggtitle(title) + 
    ylab("-log10(pvalue)") +
    geom_vline(xintercept=c(-fold_cutoff, fold_cutoff), col="red") +
    geom_hline(yintercept=-log10(pvalue_cutoff), col="darkturquoise") +
    scale_color_manual(values=c("black", "blue", "red"), 
                       name = "Differential expression\nstatus") +
    theme_minimal()
  return(plot)
  
}

library("plotly")
# Volcano: put the biological significance (fold-change)
# and statistical significance (p-value) in one plot
plot <- volcanoPlot(DEres.df, fold_cutoff=0.5, pvalue_cutoff=0.01)
ggplotly(plot)

ggsave(file=paste(getwd(),"\\Week_Isotope\\Output\\Volcano\\",paste(fileNameNQT, "Volcano.png", sep = ""), sep = ""), plot=plot, width=10, height=8)
write.csv(DEres.df, paste(getwd(),"\\Week_Isotope\\Output\\Volcano\\CSV\\",paste(fileNameNQT, "_volcano_t-test_pos.csv", sep = ""), sep = ""), row.names=FALSE)

############################################################################
#####box plot with wilcoxon and krustal-wallis stats for feature m/z 797####
library("ggpubr")
library("datatable")
library("dplyr")
library("tidyverse")
library("readr")
setwd("//SERVER/Users/james/Desktop/DESKTOP/NAUFecal/Week_Isotope/output/Volcano/CSV")
file_list <- list.files(getwd())
#for (file in file_list) {
  # Read the CSV file
  file = file_list[3]
  dataW <- read.csv(file, header = TRUE)  # Adjust this based on your CSV structure
  dataW <- t(dataW)
  
  dataW <- head(dataW, -3)
  
  colnames(dataW) <- unlist(dataW[1,])
  dataW <- dataW[-1,]
  
  
  row_names_as_column <- data.frame(group = rownames(dataW), row.names = NULL)
  dataW <- cbind(row_names_as_column, dataW)  
  rownames(dataW) <- NULL
  row_names_to_replace <- c(wtFileNames)

  for (row_name in row_names_to_replace) {
    if (row_name %in% dataW$group) {
      dataW$group[dataW$group == row_name] <- "wt"
    }
  }
  row_names_to_replace <- c(tgFileNames)
  
  for (row_name in row_names_to_replace) {
    if (row_name %in% dataW$group) {
      dataW$group[dataW$group == row_name] <- "tg"
    }
  }
  
  
  num_columns <- ncol(dataW)
  
  firstCol <- dataW[,1]
  for (i in 2:num_columns) {
    iCol <- dataW[, i]
    print(iCol)
    my_data2 <- data.frame(group = firstCol,  norm_peak_area= iCol)

p <- ggboxplot(my_data2,x="group", y="norm_peak_area", color="group", add = "jitter", shape = "group", yscale="log10", title = )
p

my_comparisons <- list( c("wt", "tg") )
fig <- p + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  stat_compare_means(label.y = .0012)
fig
BH_corr <- pairwise.wilcox.test(my_data2$norm_peak_area, my_data2$group, p.adjust.method = "BH", paired = FALSE)
m1m2sup <- BH_corr[["p.value"]]
p_m1m2sup <- c(m1m2sup[1,1], m1m2sup[2,1], m1m2sup[2,2])
paste("//SERVER/Users/james/Desktop/DESKTOP/NAUFecal/Week_Isotope/output/Volcano/Wilcoxnon/")
ggsave(file="~/Dropbox/2023-Lanthanide_Nature/R-Code_Positive/Volcano_Inputs/wilcox_boxplot.svg", plot=fig, width=10, height=8)

  }
  
  