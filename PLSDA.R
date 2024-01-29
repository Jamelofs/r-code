library(stringr)
library(data.table)
library(readr)
library(dplyr)
library(tidyr)
library(reshape)     
library(ggplot2)
library(gplots)
library(pheatmap)
library(heatmap3)
library(magrittr) 
library(tibble) 
library(RColorBrewer)
library(paletteer)
library(tidyverse)
options(install.packages.compile.from.source="never")
if (!require("pacman")) install.packages("pacman") #Installing pacman if not present
pacman::p_load("tidyverse","factoextra","KODAMA","vegan","IRdisplay","svglite")
#Global settings for plot size in the output cell:
options(repr.plot.width=10, repr.plot.height=10,res=600) # the parameters: width, height & resolution can be changed
pacman::p_load("matrixStats","ggsci","FSA","cowplot") #installing & loading the necessary libraries

setwd("\\\\SERVER\\Users\\james\\Desktop\\DESKTOP\\NAUFecal")
getwd() #to see the working directory 

metadata <- read_tsv("Cope_NAU_FecalMetadata_pos_no_wash.tsv")
scaled_table <- read_csv("IMP_clr.csv")
data_merge <- metadata %>% dplyr::select("filename", "ATTRIBUTE_Week_Strain") %>% 
  left_join(scaled_table) %>% 
  column_to_rownames("filename")



## First example
X <- data_merge %>% dplyr::select(-ATTRIBUTE_Week_Strain)# predictor variable
Y <- data_merge$ATTRIBUTE_Week_Strain

plsda_after <- mixOmics::plsda(X, Y, ncomp = 7)
indiv <- mixOmics::plotIndiv(plsda_after, ind.names = TRUE, ellipse = TRUE, legend = TRUE)

loadings <- mixOmics::plotLoadings(plsda_after, comp = 1, title = 'Loadings on comp 1', legend.color = c("royalblue", "red3"),
             contrib = 'max', method = 'mean', ndisplay = 30, size.name = 0.6, size.legend = 1, size.title = 1,xlim = c(-0.1, 0.1), ylim = c(-2, 2))
ggsave("TEST.SVG")

table2ppt(x=indiv,file=paste(filename,"_plots.pptx"),digits=4,digitspvals=1, append=TRUE) 

table2ppt(x=loadings,file=paste(Z,"_plots.pptx"),digits=4,digitspvals=1, append=TRUE) 