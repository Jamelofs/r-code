library(mixOmics) 
library(lme4) 
library(stringr)
library(data.table)
library(readr)
library(dplyr)
library(tidyr)
library(reshape)     
library(ggplot2)
library(gplots)
library(magrittr) 
library(tibble) 
library(RColorBrewer)
library(paletteer)
library(caret)
library(tidyverse)
library(Biobase)
library(gridExtra)

library(ggrepel)

# setwd("//Server/Users/james/Desktop/DESKTOP/NAUFecal")
# getwd() #to see the working directory 
# metadata_full<- read_tsv("Cope_NAU_FecalMetadata_pos_no_wash.tsv")
# scaled_table <- read_csv("Normalised_Quant_table.csv")
# GNPS <- read_tsv("GNPS\\FEATURE-BASED-MOLECULAR-NETWORKING-16d6ec53-view_all_annotations_DB-main.tsv")

setwd("//Server/Users/james/Desktop/DESKTOP/NAUSerum")
getwd() #to see the working directory 

metadata_full<- read_tsv("Cope_NAU_Metadata_serum_pos - Copy.tsv")
scaled_table <- read_csv("Normalised_Quant_table.csv")
GNPS <- read_tsv("GNPS\\FEATURE-BASED-MOLECULAR-NETWORKING-8b50a94c-view_all_annotations_DB-main.tsv")

filtered_metadata <- metadata_full %>% filter(!(ATTRIBUTE_Type == "Blank" | ATTRIBUTE_Isotope == "18O"))


scaled_table <- rbind(colnames(scaled_table), scaled_table)
colnames(scaled_table) <- NULL



transposed_data <- t(scaled_table)
colnames(transposed_data) <- transposed_data[1, ]
transposed_data <- transposed_data[-1, ]
transposed_data <- as.data.table(transposed_data)


scaled_table <- transposed_data

data_merge <- filtered_metadata %>% dplyr::select("filename", "ATTRIBUTE_Week", "ATTRIBUTE_Strain") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("filename")

column_names <- colnames(data_merge)
print(column_names)
# Rename columns starting from the third column
for (i in 3:length(column_names)) {
  split_names <- strsplit(column_names[i], "_")[[1]]
  new_name <- split_names[1]
  names(data_merge)[i] <- new_name
}



plots_list <- list( )

scans_cholic_acid <- GNPS$`#Scan#`[grepl("Tryptophan", GNPS$Compound_Name, ignore.case = TRUE)]
new_element <- "ATTRIBUTE_Strain"
scans_cholic_acid <- c(scans_cholic_acid, new_element)
new_element <- "ATTRIBUTE_Week"
scans_cholic_acid <- c(scans_cholic_acid, new_element)
pattern <- paste(scans_cholic_acid, collapse = "|")
#SingleGraphs(TRUE, compound_name)
MultiGraphs(TRUE, 4, "")


names_cholic_acid <- GNPS$Compound_Name[grepl("tryptophan", GNPS$Compound_Name, ignore.case = TRUE)]
unique_values <- unique(names_cholic_acid)
plots_list <- list( )
compound_name <- 'L-Tryptophan'
for (compound_name in unique_values) {
  print(compound_name)
  scans_cholic_acid <- GNPS$`#Scan#`[GNPS$Compound_Name == compound_name]
  new_element <- "ATTRIBUTE_Strain"
  scans_cholic_acid <- c(scans_cholic_acid, new_element)
  new_element <- "ATTRIBUTE_Week"
  scans_cholic_acid <- c(scans_cholic_acid, new_element)
  pattern <- paste(scans_cholic_acid, collapse = "|")
  #SingleGraphs(TRUE, compound_name)
  MultiGraphs(TRUE, 4, compound_name)
}










SingleGraphs <- function(save = TRUE, compound) { 
  scatter_ND <- data_merge %>% select(matches(pattern))  #before ND
  scatter_ND$ATTRIBUTE_Week<- as.numeric(scatter_ND$ATTRIBUTE_Week)
  my_plot <- scatter_ND %>%
    pivot_longer(cols = -c(ATTRIBUTE_Week, ATTRIBUTE_Strain), names_to = "Variable") %>%
    ggplot(aes(x = ATTRIBUTE_Week, y = value, color = ATTRIBUTE_Strain)) +
    ggtitle(compound_name) +
    geom_point() +
    scale_color_manual(values = c("royalblue1", "red3", "green4")) +
    labs(x = "Week", y = "Data", color = "Strain") +
    theme_bw() + 
    scale_x_continuous(breaks = unique(scatter_ND$ATTRIBUTE_Week)) +
    theme(panel.grid = element_blank()) +
    annotate("rect", xmin = 47, xmax = Inf, ymin = -Inf, ymax = Inf, 
             fill = "#C8E6C9", alpha = 0.2, color = NA, inherit.aes = FALSE) 
  plots_list <- c(plots_list, list(my_plot))
  for (i in 1:length(plots_list)) {
    print(plots_list[[i]])
  }
  if(save){
    ggsave(paste(compound,".png",sep=""), plot = my_plot, width = 6, height = 4, dpi = 300)
  }
}

MultiGraphs <- function(save = TRUE, numOfGraphs, compound = "") { 
  startNum <- 3
  scatter_ND <- data_merge %>% select(matches(pattern)) 
  scatter_ND$ATTRIBUTE_Week<- as.numeric(scatter_ND$ATTRIBUTE_Week)
  while (startNum <= ncol(scatter_ND)) {
    scatter_ND[[startNum]] <- as.numeric(scatter_ND[[startNum]]) * 10000000000
    startNum = startNum + 1
  }
  startNum <- 3
  while (startNum <= ncol(scatter_ND)) {
    
    endNum <- min(startNum + numOfGraphs-1, ncol(scatter_ND))  # Calculate end column number
    my_plot <- scatter_ND %>%
      select(ATTRIBUTE_Week, startNum:endNum, ATTRIBUTE_Strain) %>%
      pivot_longer(cols = -c(ATTRIBUTE_Week, ATTRIBUTE_Strain), names_to = "Variable") %>%
      ggplot(aes(x = ATTRIBUTE_Week, y = value, color = ATTRIBUTE_Strain)) +
      facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
      geom_point() +
      scale_color_manual(values = c("royalblue1", "red3", "green4")) +
      labs(x = "Week", y = "Data", color = "Strain") +
      theme_bw() + 
      scale_x_continuous(breaks = unique(scatter_ND$ATTRIBUTE_Week)) +
      theme(panel.grid = element_blank()) +
      annotate("rect", xmin = 47, xmax = Inf, ymin = -Inf, ymax = Inf, 
               fill = "#C8E6C9", alpha = 0.2, color = NA, inherit.aes = FALSE) 
    plots_list <- c(plots_list, list(my_plot))
    for (i in 1:length(plots_list)) {
      print(plots_list[[i]])
    }
      ggsave(paste(compound, " Multigraph",startNum-2,"-",endNum-2,".png",sep=""), plot = my_plot, width = 6, height = 4, dpi = 300)

    startNum <- endNum + 1
  }
}
