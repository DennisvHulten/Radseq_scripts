library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(cowplot)
library(phytools)
library(dendextend)
library(readr)
library("stringr")
library("adegenet")
library("vcfR")
library("treeio")
library("reshape2")
library("ggstance")
library("tidyverse")
library("gridExtra")
library("readr")
library("dplyr")

#### change for denovo or reference or dart
setwd("/Users/dvan216/Documents/Inkfish-Phd/Project_Mir/MADR_dart/D1b_building_initial_nj_trees/")

# Load and root the tree
tree <- read.nexus("reference/MADR_reference_d1.tre")
tree <- midpoint.root(tree)

# Create a ggtree object
tree_plot <- ggtree(tree, layout = "rectangular", size = 0.15, open.angle = 5) 


# Add color information to tree_plot$data
tree_plot$data <- tree_plot$data %>%
  mutate(color = case_when(
    grepl("SNA$", label) ~ "firebrick",                         
    grepl("SEA$", label) ~ "steelblue",                        
    grepl("KAL$", label) ~ "mediumpurple",                       
    grepl("SNA_SHA$", label) ~ "coral3",              
    grepl("SEA_SHA$", label) ~ "royalblue2",               
    grepl("KAL_SHA$", label) ~ "plum",
    grepl("HEL$", label) ~ "orange",
    TRUE ~ "black"                                        
  ))

# Add color to the tip labels
p <- tree_plot +
  geom_tippoint(size=1,aes(color = color)) +
  geom_tiplab(size=1,align =FALSE, linetype = "dotted", linesize =0.1, aes(color=color)) + # Map color column to the tip labels
  scale_color_identity()             # Use exact colors without mapping

# Display the tree
p

MADR_reference_het <- read_delim("reference/MADR_reference_d1.het", delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
MADR_reference_miss <- read_delim("reference/MADR_reference_d1.imiss", delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
MADR_reference_miss$GENOTYPED <- MADR_reference_miss$N_DATA - MADR_reference_miss$N_MISS

#MADR_reference_frq <- read_delim("reference/MADR_reference_d1.frq", delim = "\t", escape_double = FALSE,  trim_ws = TRUE)
MADR_reference_depth <- read_delim("reference/MADR_reference_d1.idepth", delim = "\t", escape_double = FALSE,  trim_ws = TRUE)

d1 <- data.frame(INDV=factor(tree$tip.label))
d1$F <- with(MADR_reference_het, F [match(d1$INDV, INDV)])
d1$het <- with(MADR_reference_het, ((N_SITES - `O(HOM)`)/N_SITES) [match(d1$INDV, INDV)])
d1$miss <- with(MADR_reference_miss, F_MISS [match(d1$INDV, INDV)])
d1$genotyped <- with(MADR_reference_miss, GENOTYPED [match(d1$INDV, INDV)])
d1$depth <- with(MADR_reference_depth, MEAN_DEPTH [match(d1$INDV, INDV)])
tree_tips <- tree_plot$data$label[tree_plot$data$isTip]
d1 <- d1[match(tree_tips, d1$INDV), ]

p1 <- p + 
  geom_facet(
    panel = "Genotyped", 
    data = d1, 
    geom = geom_barh, 
    mapping = aes(x = genotyped),
    width = 1, stat = 'identity'
  ) + 
  scale_x_continuous(expand = c(0, 0))  # Remove extra padding on the x-axis
p1

p2 <- p1 + 
  geom_facet(
    panel = "Heterozygosity", 
    data = d1, 
    geom = geom_barh, 
    mapping = aes(x = het),
    width = 1, stat = 'identity')

p2

p3 <- p2 + geom_facet(
  panel = "Average sequencing depth", 
  data = d1, 
  geom = geom_barh, 
  mapping = aes(x = depth),
  width = 1, stat = 'identity'
)
p3
