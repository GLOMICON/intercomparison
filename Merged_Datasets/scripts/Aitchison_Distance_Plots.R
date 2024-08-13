#Aitchison_distance_plots
# Katie Pitz
# 08/12/24

# Load Libraries -----------------------------------------------------------------
library(qiime2R)
library(tidyverse)
library(lubridate) #for date modifications
library(vegan)
library(ggthemes)
library(magrittr)
#for heatmap section
library(gridExtra)
library(RColorBrewer) #colors for plotting
library(forcats) #working with factor data (reordering)
library(patchwork) #putting graphs together on same plot

library(ggdendro)

# Set Constants -----------------------------------------------------------------

marker = sym("18S")
prefix = "GLOMICON"
# Set directory to save plots
plot_directory <- './Merged_Datasets/figures/Aitchison_distance/'
# Set directory to retrieve data
data_directory = "./Merged_Datasets/RPCA/"

Analyzing_Institutes <- c('AWI', 'SBR', 'MBARI')


# Loop over Analyzing Institutes ----------------------------------------------

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
  
  ## Distance between Collecting Institute Samples  --------------------------------------------------------------------

  #create df for plotting dist between points
  df <- as_tibble(aitmat1) %>%
    pivot_longer(-sample_name, names_to = "sample_name.y", values_to = "distance") %>%
    left_join(scores_meta %>% select(Collecting_Institute, sample_name)) %>%
    left_join(scores_meta %>% select(Collecting_Institute, sample_name), by = c("sample_name.y" = "sample_name")) %>%
    rename('sample_name.x' = 'sample_name') %>%
    filter(sample_name.x != sample_name.y) %>%
    filter(Collecting_Institute.x %in% c('BLOOMMOCK', 'EVENMOCK')==FALSE) %>%
    filter(Collecting_Institute.y %in% c('BLOOMMOCK', 'EVENMOCK')==FALSE) %>%
    mutate(distance = as.numeric(distance)) %>%
    drop_na()
    
  # Plot
  p<- ggplot(data=df, aes(x=Collecting_Institute.x, y=distance,fill = Collecting_Institute.y)) +
    geom_boxplot() + 
    facet_wrap(~ Collecting_Institute.x, scales = "free_x") +
    theme_minimal()
  p
  filename = paste(plot_directory, 'Aitdist_',Analyzing_Institute,'_boxplot.png', sep='')
  filename
  ggsave(filename,height = 8, width =8, units = 'in')
  
}

