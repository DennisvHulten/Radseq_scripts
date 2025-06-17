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
library("adegenet")
library("vcfR")
library("poppr")
library("dplyr")
library("reshape2")
library("tidyr")
library("hierfstat")
library("igraph")
library("treeio")
library("phytools")
library("RColorBrewer")
library("scales")

#### change for denovo or reference or dart
setwd("/Users/dvan216/Documents/Inkfish-Phd/Project_Mir/MADR_dart/D2_clone_detection_analysis/reference/")

# Load and root the tree
tree <- read.nexus("reference/MADR_denovo_d1.tre")
tree <- midpoint.root(tree)

# Create a ggtree object
tree_plot <- ggtree(tree, layout = "rectangular", size = 0.15, open.angle = 5) 


#looking at highest matches and missing data
#for detect_clones_vcf output
#highest_matches <- read_delim("~/Documents/Inkfish-Phd/Project_Mir/MADR_dart/D2_clone_detection_analysis/dart/MADR_dart_d1_highest_matches.txt", 
                                             #delim = "\t", escape_double = FALSE, 
                                             #trim_ws = TRUE)
#pairs <- highest_matches %>% separate(samples, into = c("indv_1", "indv_2"), sep = " vs ")
#for vcf_clone_detect output
pairs <- read_csv("MADR_reference_d2_genetic_similarity.txt")

png("genetic_similarity_histogram.png", width = 800, height = 600)
# Create a histogram and suppress default x-axis
# Create a histogram and capture its breaks
hist_data <- hist(pairs$match_perc, 
                  breaks = seq(floor(min(pairs$match_perc)), ceiling(max(pairs$match_perc)), by = 1), 
                  col = "skyblue", 
                  main = "Reference d2 Genetic Similarity", 
                  xlab = "Genetic Similarity", 
                  ylab = "Frequency", 
                  xaxt = "n",
                  ylim = c(0, 5000)) # Suppress default x-axis

# Use the breaks of the histogram as tick positions
tick_positions <- hist_data$breaks  # Breaks used in the histogram
axis(1, at = tick_positions, labels = tick_positions)  # Add custom x-axis with matching labels
dev.off()

#1. Filter the data above the average of similarity in technical replicates
threshold <- 99
#could be issues if match_perc is not numeric
#pairs$match_perc <- as.numeric(as.character(pairs$match_perc))
above_99 <- pairs %>% filter(match_perc >= threshold)
#write.table(above_99.5, here("clones_overall_99.5.txt"), quote = FALSE, col.names = TRUE, sep = "\t")

#Visualize the similarity
similarity <- graph_from_data_frame(
  d = above_99[, c("ind1", "ind2")], 
  directed = FALSE)

similarity

#All sites
png("MADR_reference_d2_clonality_across_locations.png", width = 800, height = 800)

par(mfrow = c(1, 1))
plot(similarity,
     vertex.size = 5,          # Size of nodes
     vertex.label.cex = 0.4,    # Size of node labels
     vertex.color = "lightblue",# Color of nodes
     edge.color = "gray",       # Color of edges
     edge.width = 2,            # Width of edges
     main = "Clonality overall")

dev.off()

#2. Visualize by location
#Keep in mind there are some clone pairs across sites that are not visible in
#these graphs by location
#Save
png("MADR_reference_d2_clonality_within_locations.png", width = 1200, height = 400)
par(mfrow = c(1, 3))

# Split the "pair" column into two populations and filter
above_99_KAL <- above_99 %>%
  mutate(
    pop1 = sub("-.*", "", pop),      # Extract the part before "-"
    pop2 = sub(".*-", "", pop)       # Extract the part after "-"
  ) %>%
  filter(substr(pop1, 1, 3) == "KAL" & substr(pop2, 1, 3) == "KAL" ) 

sim_KAL <- graph_from_data_frame(
  d = above_99_KAL[, c("ind1", "ind2")], 
  directed = FALSE)

clon_KAL <- plot(sim_KAL,
                 vertex.size = 5,          # Size of nodes
                 vertex.label.cex = 0.4,    # Size of node labels
                 vertex.color = "lightblue",# Color of nodes
                 edge.color = "gray",       # Color of edges
                 edge.width = 2,            # Width of edges
                 main = "Kalki")      # Width of edges

#Snakebay
above_99_SNA <- above_99 %>%
  mutate(
    pop1 = sub("-.*", "", pop),      # Extract the part before "-"
    pop2 = sub(".*-", "", pop)       # Extract the part after "-"
  ) %>%
  filter(substr(pop1, 1, 3) == "SNA" & substr(pop2, 1, 3) == "SNA" ) 

sim_SNA <- graph_from_data_frame(
  d = above_99_SNA[, c("ind1", "ind2")], 
  directed = FALSE)

clon_SNA <- plot(sim_SNA,
                 vertex.size = 5,          # Size of nodes
                 vertex.label.cex = 0.4,    # Size of node labels
                 vertex.color = "lightblue",# Color of nodes
                 edge.color = "gray",       # Color of edges
                 edge.width = 2,            # Width of edges
                 main = "Snakebay")      # Width of edges

#Seaquarium
above_99_SEA <- above_99 %>%
  mutate(
    pop1 = sub("-.*", "", pop),      # Extract the part before "-"
    pop2 = sub(".*-", "", pop)       # Extract the part after "-"
  ) %>%
  filter(substr(pop1, 1, 3) == "SEA" & substr(pop2, 1, 3) == "SEA" ) 

sim_SEA <- graph_from_data_frame(
  d = above_99_SEA[, c("ind1", "ind2")], 
  directed = FALSE)

clon_SEA <- plot(sim_SEA,
                 vertex.size = 5,          # Size of nodes
                 vertex.label.cex = 0.4,    # Size of node labels
                 vertex.color = "lightblue",# Color of nodes
                 edge.color = "gray",       # Color of edges
                 edge.width = 2,            # Width of edges
                 main = "Seaquarium")      # Width of edges
dev.off()


# clones between locations
# Split the "pair" column into two populations and filter
above_99_between_sites <- above_99 %>%
  mutate(
    pop1 = sub("-.*", "", pop),      # Extract the part before "-"
    pop2 = sub(".*-", "", pop)       # Extract the part after "-"
  ) %>%
  filter(substr(pop1, 1, 3) != substr(pop2, 1, 3)) 

#Visualize the similarity
similarity <- graph_from_data_frame(
  d = above_99_between_sites[, c("ind1", "ind2")], 
  directed = FALSE)

similarity

#All sites
png("MADR_reference_d2_clonality_between_locations.png", width = 800, height = 800)

par(mfrow = c(1, 1))
plot(similarity,
     vertex.size = 5,          # Size of nodes
     vertex.label.cex = 0.4,    # Size of node labels
     vertex.color = "lightblue",# Color of nodes
     edge.color = "gray",       # Color of edges
     edge.width = 2,            # Width of edges
     main = "Clonality overall")

dev.off()

# Combine both columns and get unique values
unique_individuals <- unique(c(above_99_between_sites$ind1, above_99_between_sites$ind2))
between_sites_clones
#######identifying MLG's using poppr
## Import files and convert to poppr
#Individuals
ind <- read.delim("MADR_reference_popfile_d2.txt", sep = "\t", header = FALSE)
colnames(ind)[1] <- "Sample"
colnames(ind)[2] <- "Location"


genind <- vcfR2genind(read.vcfR("MADR_denovo_d2.vcf"),return.alleles = TRUE)
genind
#head(pop(genind))
tibble <- as_tibble(genind@tab)

#Create a matrix to add strata and add pop and strata
strata <- ind
strata <- strata[,-1]
# Create a data frame from the strata vector
strata_df <- data.frame(pop = strata)
genind@pop <- as.factor(ind$Location)
strata(genind) <- strata_df
genind

head(pop(genind))
summary(genind)
names(genind)

#Check number of alleles per locus
# Open a PNG device
png("number_of_alleles_per_locus_dart.png", width = 800, height = 600)

# Create your plot
barplot(genind$loc.n.all, 
        ylab = "Number of alleles",
        main = "Number of alleles per locus")

# Close the device
dev.off()

#Convert to genclone
genclone <- poppr::as.genclone(genind)


#Calculate genetic diversity measures
#Naive ("original") multilocus genotype calculation. 

mll(genclone) <- "original"      #function displays the assignment of the multilocus lineages in the sample.
nmll(genclone)    #number of multilocus lineages in the sample
#[1] 453
genclone

poppr(genind, 
      #sample=999
      missing = "ignore")


#Filtered ("contracted")
#Check the distribution of the MLG before filtering
mlg.table(genclone)

#Check the distance between individuals and using the previous knowledge about clones, identify the threshold
#(If your maps are already annotated, currently not annotated for Scolymia or Mycetophyllia - 25/11/24)
#This calculates Ecuclidean and dissimilarity matrices
dist <- as.data.frame(bitwise.dist(genclone, percent=TRUE, mat=TRUE, missing_match = FALSE, scale_missing=FALSE,euclidean = TRUE))

#Check the values of technical replicates
max(dist)
#[1] 0.2927653
min(dist)
#[1] 0


#####matchin per genotyped
max_sites_den <- 98416
max_sites_ref <- 121930
max_sites_dar <- 43903

for (i in nrow(clones)) {
  max_genotyped <- min(clones$total_snps_indv_1
  
}


