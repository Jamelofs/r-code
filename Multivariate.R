#calling the necessary packages:
library(ggplot2)
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
library(cowplot) #efficient functions to arrange several plots

#read in clr-transformed feature table and metadata
setwd("//Server/Users/james/Desktop/DESKTOP/NAUFrontalCortex_LowerNoise")
Imp_clr <- read.csv("IMP_clr.csv", header = T, check.names = F, sep = ",")
metadata <- read_tsv("Cope_NAU_Metadata_FrontalCortex_pos.tsv")
week <- 24
isotope <- 16

filtered_metadata <- metadata %>% filter(ATTRIBUTE_Isotope == isotope & ATTRIBUTE_Week == week)


#merge metadata column with clr-transformed feature table
data_merge <- filtered_metadata %>% dplyr::select("filename", "ATTRIBUTE_Strain") %>% 
  left_join(Imp_clr) %>% 
  column_to_rownames("filename")

column_names <- colnames(data_merge)

for (i in 2:length(column_names)) {
  split_names <- strsplit(column_names[i], "_")[[1]]
  new_name <- split_names[1]
  names(data_merge)[i] <- new_name
}


#make scree plot 
res.pca <- data_merge %>%
  dplyr::select(-ATTRIBUTE_Strain) %>% 
  prcomp(scale = FALSE) #TRUE if scaling was not performed before
fviz_eig(res.pca)

#plot pca plot
pca <- fviz_pca_ind(res.pca, col.ind = (data_merge)$ATTRIBUTE_Strain, addEllipses = TRUE, label = "none")
pca
ggsave(file="W8O16.svg", plot=pca, width=10, height=8)

#output loadings plot
pc_loadings <- res.pca$rotation
pc_loadings <- as_tibble(pc_loadings, rownames = "feature")
write.csv(pc_loadings, "pca_loadings_all_timeopoints_all_isotopes.csv",row.names = TRUE)

#plot driving features
tmp <- fviz_pca_var(res.pca, select.var = list(contrib = 25), col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
write.csv(tmp$data['name'], "chr_output.csv", row.names = FALSE)

# PERMANOVA
metadata <- data_merge[, 1]
metadata_df <- data.frame(metadata = metadata)
metabolites <- data_merge[, 2:ncol(data_merge)]
dist_metabolites <- vegdist(metabolites, method = "euclidean", na.rm = TRUE)
permanova_all <- adonis2(dist_metabolites ~ metadata_df$metadata, metadata_df, na.action = na.omit)

