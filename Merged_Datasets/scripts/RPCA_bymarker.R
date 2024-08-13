# Plot RPCA Results (Gemelli) - by Analyzing Institute
# Katie Pitz
# 08/12/24

# Load Libraries -----------------------------------------------------------------
library(readr) #read csv files
library(lubridate) #for date modifications
library(dplyr)
library(ggplot2)
library(ggthemes)
library(magrittr)
library(tidyr)
library(RColorBrewer) #colors for plotting
library(forcats) 
library(stringr)
library(viridis) #color maps
library(wesanderson) #colors

# Set Constants -----------------------------------------------------------------

marker = sym("18S")
prefix = "GLOMICON"
# Set directory to save plots
plot_directory <- './Merged_Datasets/figures/RPCA/'
# Set directory to retrieve data
data_directory = "./Merged_Datasets/RPCA/"

# Loop over Analyzing Institutes ----------------------------------------------
Analyzing_Institutes <- c('AWI', 'SBR', 'MBARI')
for (val in Analyzing_Institutes) {
  Analyzing_Institute = sym(val)
  
  ## Import RPCA Data -------------------------------------------------------------
  
  file = "_scores.csv"
  filepath = paste(data_directory,Analyzing_Institute,'/',Analyzing_Institute, file, sep='')
  print(filepath)
  scores <- read_csv(filepath) %>% rename('sample_name' = 1)
  
  file = "_loadings.csv"
  filepath = paste(data_directory,Analyzing_Institute,'/',Analyzing_Institute, file, sep='')
  print(filepath)
  loadings <- read_csv(filepath) %>% rename('ASV_ID' = 1)
  
  file = "_prop_explained.csv"
  filepath = paste(data_directory,Analyzing_Institute,'/',Analyzing_Institute, file, sep='')
  print(filepath)
  prop_explained <- read_csv(filepath)
  
  
  file = "_distance.csv"
  filepath = paste(data_directory,Analyzing_Institute,'/',Analyzing_Institute, file, sep='')
  print(filepath)
  distance <- read_csv(filepath) %>% rename('sample_name' = 1)
  #need to get that sample_name column as rownames in matrix...
  aitmat1 <- as.matrix(distance)
  aitmat <- aitmat1[,-1]
  rownames(aitmat) <- aitmat1[,1]
  
  aitdist <- as.dist(aitmat)
  
  ## Bin Metadata. ----------------------------------------------------------
  
  # create replicate numbers for plotting
  scores_meta <- scores %>%
    group_by(Analyzing_Institute, Collecting_Institute) %>%
    mutate(replicateID = row_number()) %>%
    ungroup()
  
  ## Cluster ----------------
  
  aitclust <- hclust(aitdist, method= "ward.D2")
  plot(aitclust)
  aitdend <- as.dendrogram(aitclust, hang = -1, lwd = 3, lty = 3, sub = "")
  #dend_cols <- select(samp.c,c('SampleID','date'))
  plot(aitdend)
  
  #Can plot with ggplot2
  library(ggdendro)
  dendr <- dendro_data(aitclust, type="rectangle") 
  
  #Define cluster based on similarity percent - good for plotting
  clust <- cutree(aitclust, h = 1)               # find 'cut' clusters (k) or choose level (h) numeric scalar or vector with heights where the tree should be cut.
  clust2 <- cutree(aitclust, h = 2)
  clust3 <- cutree(aitclust, h = 3)  
  clust4 <- cutree(aitclust, h = 5)  
  
  clust.df <- data.frame(label = names(clust), cluster_6 = clust, cluster_7 = clust2, cluster_10 = clust3, cluster_16 = clust4)
  sapply(clust.df, mode)
  clust.df$label <- as.character(clust.df$label)
  clust.df <- clust.df[order(clust.df$label),]
  
  #join sample data with cluster df and simprofCLUSTERS df:
  tree_scores <- as_tibble(clust.df) %>% rename(sample_name = label) %>% left_join(scores_meta, by='sample_name')
  tree_data <- as_tibble(dendr$labels) %>% 
    rename(sample_name = label) %>% 
    left_join(tree_scores, by='sample_name') %>% 
    arrange(Collecting_Institute)
  
  
  ggplot() + 
    geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), size=0.2, color='darkgrey') +
    #geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=year_char),alpha=1, size=1)+
    geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-1,ymax=y-2, color=Collecting_Institute),alpha=1, size=1)+
    #geom_rect(data=tree_data, aes(xmin=x,xmax = x+.01, ymin=y-2,ymax=y-3, color=ESP),alpha=1, size=1)+
    geom_text(data=tree_data, aes(x=x, y=y, label=Collecting_Institute, hjust=-.4), size=2) +
    geom_text(data=tree_data, aes(x=x, y=y, label=sample_name, hjust=-2.5), size=1) +
    #geom_text(data=tree_data, aes(x=x, y=y, label=libraryID, hjust=-6), size=2) +
    coord_flip() + 
    scale_y_reverse(expand=c(.5, 0)) + 
    # scale_colour_manual(values= c( "#6baed6", "#3182bd", "#08519c", 
    #                                "#a1d99b", "#74c476", "#31a354", "#006d2c", 
    #                                "#fdae6b", "#fd8d3c", "#e6550d", "#a63603","#9ecae1"))+
    theme(axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank())
  
  filename = paste(plot_directory,Analyzing_Institute, '_dendrogram_Collecting_Institute.png', sep='')
  print(filename)
  ggsave(filename,height = 5, width =8, units = 'in')
  
  
  ## Plot RPCA --------------------------------------------------------------------
  
  scores_meta %>%
    ggplot(aes(x=PC1,y=PC2, color=Collecting_Institute)) +
    geom_point() +
    #scale_color_manual(values = wes_palette("Darjeeling1")) +
    scale_color_brewer(palette='Dark2') +
    theme_minimal()
  
  filename = paste(plot_directory,Analyzing_Institute, '_RPCA_PC1PC2_Collecting_Institute.png', sep='')
  print(filename)
  ggsave(filename,height = 5, width =8, units = 'in')
  
  scores_meta %>%
    ggplot(aes(x=PC2,y=PC3, color=Collecting_Institute)) +
    geom_point() +
    #scale_color_manual(values = wes_palette("Darjeeling1")) +
    scale_color_brewer(palette='Dark2') +
    theme_minimal()
  
  filename = paste(plot_directory,Analyzing_Institute, '_RPCA_PC2PC3_Collecting_Institute.png', sep='')
  print(filename)
  ggsave(filename,height = 5, width =8, units = 'in')
  
  
  
}

