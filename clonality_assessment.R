#Clonality assessment based on pairwise comparisons and poppr
#Alejandra Hernandez

## Dependencies ==============================================================
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


#############################################################
# Set the main working directory:
setwd("/Users/alejhernandez/Dropbox/Rare species/scolymia/dart/C - Filtering and QC of loci and SNPs/C4 - Detect and eliminate clones")



#############################################################
#Open data
clones <- read.delim("sco_clones_c4_80.txt", sep = ",", header = TRUE)
colnames(clones)[1] <- "ind1"


################################################################################################################
#Evaluate the clonality in pairwise comparisons above the average of similarity in technical replicates
str(clones)

#1. Filter the data above the average of similarity in technical replicates
threshold <- 99.5
above_99.5 <- clones %>% filter(match_perc >= threshold)
#write.table(above_99.5, here("clones_overall_99.5.txt"), quote = FALSE, col.names = TRUE, sep = "\t")

#Visualize the similarity
similarity <- graph_from_data_frame(
  d = above_99.5[, c("ind1", "ind2")], 
  directed = FALSE)

similarity

#All sites
png("clonality_across_locations.png", width = 800, height = 800)

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
png("clonality_within_locations.png", width = 1200, height = 800)
par(mfrow = c(2, 3))

#Kalki
above_99.5_KAL <- subset(above_99.5, pop=="KAL")
sim_KAL <- graph_from_data_frame(
  d = above_99.5_KAL[, c("ind1", "ind2")], 
  directed = FALSE)

clon_KAL <- plot(sim_KAL,
                 vertex.size = 5,          # Size of nodes
                 vertex.label.cex = 0.4,    # Size of node labels
                 vertex.color = "lightblue",# Color of nodes
                 edge.color = "gray",       # Color of edges
                 edge.width = 2,            # Width of edges
                 main = "Kalki")      # Width of edges

#Boka Hulu
above_99.5_HUL <- subset(above_99.5, pop=="HUL")
sim_HUL <- graph_from_data_frame(
  d = above_99.5_HUL[, c("ind1", "ind2")], 
  directed = FALSE)

clon_HUL <- plot(sim_HUL,
                 vertex.size = 5,          # Size of nodes
                 vertex.label.cex = 0.4,    # Size of node labels
                 vertex.color = "lightblue",# Color of nodes
                 edge.color = "gray",       # Color of edges
                 edge.width = 2,            # Width of edges
                 main = "Boka Hulu")    

#Coral Estate
above_99.5_EST <- subset(above_99.5, pop=="EST")
sim_EST <- graph_from_data_frame(
  d = above_99.5_EST[, c("ind1", "ind2")], 
  directed = FALSE)

clon_EST <- plot(sim_EST,
                 vertex.size = 5,          # Size of nodes
                 vertex.label.cex = 0.4,    # Size of node labels
                 vertex.color = "lightblue",# Color of nodes
                 edge.color = "gray",       # Color of edges
                 edge.width = 2,            # Width of edges
                 main = "Coral Estate")   

#Snake Bay
above_99.5_SNA <- subset(above_99.5, pop=="SNA")
sim_SNA <- graph_from_data_frame(
  d = above_99.5_SNA[, c("ind1", "ind2")], 
  directed = FALSE)

clon_SNA <- plot(sim_SNA,
                 vertex.size = 5,          # Size of nodes
                 vertex.label.cex = 0.4,    # Size of node labels
                 vertex.color = "lightblue",# Color of nodes
                 edge.color = "gray",       # Color of edges
                 edge.width = 2,            # Width of edges
                 main = "Snake Bay")   

#Marie Pampoen
above_99.5_MAR <- subset(above_99.5, pop=="MAR")
sim_MAR <- graph_from_data_frame(
  d = above_99.5_MAR[, c("ind1", "ind2")], 
  directed = FALSE)

clon_MAR <- plot(sim_MAR,
                 vertex.size = 5,          # Size of nodes
                 vertex.label.cex = 0.4,    # Size of node labels
                 vertex.color = "lightblue",# Color of nodes
                 edge.color = "gray",       # Color of edges
                 edge.width = 2,            # Width of edges
                 main = "Marie Pampoen")  


#Director's Bay
above_99.5_DIR <- subset(above_99.5, pop=="DIR")
sim_DIR <- graph_from_data_frame(
  d = above_99.5_DIR[, c("ind1", "ind2")], 
  directed = FALSE)

clon_DIR <- plot(sim_DIR,
                 vertex.size = 5,          # Size of nodes
                 vertex.label.cex = 0.4,    # Size of node labels
                 vertex.color = "lightblue",# Color of nodes
                 edge.color = "gray",       # Color of edges
                 edge.width = 2,            # Width of edges
                 main = "Director's Bay")  

# Close the device to save the file
dev.off()



#3. Find groups (connected components)
clone_groups <- components(similarity)$membership

# Create a data frame with group information
clone_groups_data <- data.frame(
  sample = names(clone_groups),
  group = clone_groups)


#4. Calculate the average similarity for each group
# Merge the group info back to the filtered data
above_99.5 <- above_99.5 %>%
  mutate(
    group1 = clone_groups[as.character(ind1)],
    group2 = clone_groups[as.character(ind2)]
  ) %>%
  filter(group1 == group2)  # Ensure the pairs are in the same group

# Calculate the average similarity for each group
group_averages <- above_99.5 %>%
  group_by(group1) %>%
  summarise(average_similarity = mean(match_perc)) %>%
  rename(group = group1)

# View the groups and their average similarity
print(group_averages)



################################################################################################################
#Asses clonality with poppr


## Import files and convert to poppr
#Individuals
ind <- read.delim("sco_popfile_c3_80.txt", sep = "\t", header = FALSE)
colnames(ind)[1] <- "Sample"
colnames(ind)[2] <- "Location"


genind <- vcfR2genind(read.vcfR("sco_c3_80.vcf"),return.alleles = TRUE)
genind
#head(pop(genind))
tibble <- as_tibble(genind@tab)

#Create a matrix to add strata and add pop and strata
strata <- ind
strata <- strata[,-1]


genind@pop <- as.factor(ind$Location)
strata(genind) <- strata
genind

head(pop(genind))
summary(genind)
names(genind)

#Check number of alleles per locus
barplot(genind$loc.n.all, ylab="Number of alleles",
        main="Number of alleles per locus")


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

#SCCESTS01418
dist["SCCESTS01418","SCCESTSX1418"]
#[1] 0.002172208                         #This is the highest, and could be the threshold

#SCCDIRS01729
dist["SCCDIRS01729","SCCDIRSX1729"]
#[1] 0.001116273

#SCCHULS01792, SCCHULSX1792
dist["SCCHULS01792","SCCHULSX1792"]
#[1] 0.001372715 

#SCCESTS01426, SCCESTSX1426
dist["SCCESTS01426","SCCESTSX1426"]
#[1] 0.001086104

#SCCMARS01890, SCCMARSX1890
dist["SCCMARS01890","SCCMARSX1890"]  
#[1] 0.001780004

#SCCESTS01434, SCCESTSX1434
dist["SCCESTS01434","SCCESTSX1434"] 
#[1] 0.00104085

#SCCKALS01676, SCCKALSX1676
dist["SCCKALS01676","SCCKALSX1676"] 
#[1] 0.01956496

#Run in parallel the three different clustering methods
#With missing = "ignore" does not remove or replace missing data.
genclone_stats <- filter_stats(genclone,
                                    distance = bitwise.dist,
                                    stats = "All",
                                    missing = "ignore",
                                    plot = TRUE,
                                    nclone = 70,    #Based on the run pairwise similarity (Pim's script)
                                    hist = "Scott")


#The options are
#"farthest_neighbor" - in red in the figure
#(default) merges clusters based on the maximum distance between points in either cluster. This is the strictest of the three.

#"nearest_neighbor" - Green in the graph
#merges clusters based on the minimum distance between points in either cluster. This is the loosest of the three.

#"average_neighbor" - in blue
#merges clusters based on the average distance between every pair of points between clusters.

#One method described in the literature of choosing a threshold is to look for an initial, small peak in the histogram of 
#pairwise genetic distances and set the threshold to be between that peak and the larger peak `(Arnaud-Haond et al. 2007, 
#@bailleul2016rclone). This initial peak likely represents clones differentiated by a small set of random mutations. But
#in this case, the largest is at a small genetic distance.

#However, if this peak is not obvious, then another method is to look for the largest gap between all putative thresholds. 
#For this, you can use the cutoff_predictor() function with the output of filter_stats(). It should be noted that this method is 
#not a perfect solution. If we take the results from above, we can find the threshold for each algorithm:

print(farthest_thresh <- cutoff_predictor(genclone_stats$farthest$THRESHOLDS))
#[1] 0.003107464
print(average_thresh <- cutoff_predictor(genclone_stats$average$THRESHOLDS))
#[1] 0.003142848
print(nearest_thresh <- cutoff_predictor(genclone_stats$nearest$THRESHOLDS))
#[1] 0.0004525433

#Our clonality threshold seats close to the average, so I'll move forward with with that one

mlg.filter(genclone,
           missing = "ignore",
           algorithm = "average",
           distance = bitwise.dist
) <- average_thresh

genclone
# This is a genclone object
# -------------------------
#   Genotype information:
#   
#   234 contracted multilocus genotypes
# (0.003) [t], (bitwise.dist) [d], (average) [a] 
# 453 diploid individuals
# 33146 codominant loci
# 
# Population information:
#   
#   0 strata. 
# 7 populations defined - DIR, EST, HUL, KAL, MAR, SEA, SNA

#Access to the clone groups 
samples <- ind
poppr_clones <- as.data.frame(genclone@mlg[,"contracted"])

poppr_clones <- bind_cols(samples, poppr_clones)
colnames(poppr_clones)[3] <- "clonal_pairs"

#Check duplicates
unique(poppr_clones$clonal_pairs)
#234
n_occur <- data.frame(table(poppr_clones$clonal_pairs))
n_dup <- n_occur %>% group_by(Freq) %>% summarize(Total = sum(Freq))
n_occur[n_occur$Freq > 1,]

# write.table(poppr_clones, file = "sco_poppr_clones.txt", sep = "\t",
#           row.names = FALSE, col.names = TRUE)


#Adding poppr clones to the list created from the pairwise comparison
#Momentarily changing the name of column to match the two datasets
colnames(above_99.5)[1] <- "Sample"
above_99.5$clonal_pairs <- with(poppr_clones,clonal_pairs[match(above_99.5$Sample,Sample)])
colnames(above_99.5)[1] <- "ind1"
colnames(above_99.5)[11] <- "clonal_pairs1"

colnames(above_99.5)[2] <- "Sample"
above_99.5$clonal_pairs <- with(poppr_clones,clonal_pairs[match(above_99.5$Sample,Sample)])
colnames(above_99.5)[2] <- "ind2"
colnames(above_99.5)[12] <- "clonal_pairs2"

# write.table(above_99.5, file = "sco_clonality_assessment_99.5.txt", sep = "\t",
#           row.names = FALSE, col.names = TRUE)


#Wrangling to visualize it in the tree
group1 <- above_99.5[,c(1,9)]
group2 <- above_99.5[,c(2,10)]

#Substitute column names to merge
colnames(group1)[1] <- "Sample"
colnames(group2)[1] <- "Sample"

colnames(group1)[2] <- "group"
colnames(group2)[2] <- "group"

#Merge and make unique
clonal_groups <- bind_rows(group1, group2)
clonal_groups <-  clonal_groups %>% distinct(Sample, .keep_all = TRUE)
  

###############################################################################################################
#Highlight in the tree the clonal groups for which the two methods coincide
tree <- read.nexus("/Users/alejhernandez/Dropbox/Rare species/scolymia/dart/C - Filtering and QC of loci and SNPs/C4 - Detect and eliminate clones/sco_c3_80.tre")

#Modify root
tree <- midpoint.root(tree)
tree

#Add missing data
missing <- read_table2("/Users/alejhernandez/Dropbox/Rare species/scolymia/dart/C - Filtering and QC of loci and SNPs/C3 - Sample performance and NJ tree/sco_performance_c3_80.txt")

#Modifying column names
missing <- missing[,-6]
colnames(missing)[5] <- "Genotyped"


# Generate 60 distinct colors
colors <- colorRampPalette(brewer.pal(12, "Set3"))(70)

theme2 <- theme(axis.text = element_text(size=8),axis.text.x = element_blank(),
                axis.title = element_blank(),legend.position="none",strip.text = element_blank(),
                strip.background = element_blank(),axis.line.x = element_line(color="black", size = 0.3))


#Plot
tree_plot <- ggtree(tree, layout = "rectangular", branch.length = "TRUE", ladderize = T)
tree_plot

p1 <- tree_plot %<+% clonal_groups +
  geom_tiplab(aes(color=factor(group)),size=0.9, geom="text")+
  scale_colour_manual(values=colors, na.value = "black")+
  theme2

p1 

p2 <- facet_plot(p1, panel = "Genotyped",
                 data = missing, geom = geom_barh,
                 mapping = aes(x = missing$Genotyped),
                 stat='identity', width = 1)

p2

#ggsave(plot = p2, filename = "tree_clonality_coincide.png", dpi = 300, limitsize = FALSE, width = 10, height = 16)



#Highlighting inconsistencies:
#Select the 10 groups with inconsistencies (using vcf_clone_detect clonal group)
inconsistent <- c(2, 8, 11, 24, 35, 60, 65, 68, 69, 70)
inconsistent2 <- subset(above_99.5, group1 %in% inconsistent)
str(inconsistent2)

#Iterate between the levels of inconsistent groups to identify number of colors necessary
group_inconsistent <- subset(inconsistent2, group1 == 70)
str(group_inconsistent)
group_inconsistent$clonal_pairs1 <- as.factor(group_inconsistent$clonal_pairs1)
levels(group_inconsistent$clonal_pairs1)

#Repeat previous steps with group2
inconsistent3 <- subset(above_99.5, group2 %in% inconsistent)
str(inconsistent3)
group_inconsistent2 <- subset(inconsistent3, group1 == 70)
group_inconsistent2$clonal_pairs2 <- as.factor(group_inconsistent2$clonal_pairs2)
levels(group_inconsistent2$clonal_pairs2)


colors2 <- c("143"="#e9c46a", "150"="#ee9b00",
            "196"="#333d29","227"="#656d4a","433"="#a4ac86",
            "51"="#bb3e03","107"="#bb3e03","115"="#ca6702",
            "15"="#ca7df9", "318"="#e0aaff",
            "156"="#03045e","269"="#0077b6","270"="#00b4d8","281"="#90e0ef",
            "175"="#d90429","350"="#ef233c",
            "254"="#ff5d8f","345"="#ff97b7",
            "253", "#7cb518","397"="#5c8001",
            "302"= "#00f5d4","240"="#246a73",
            "407"="#99582a", "323"="#772f1a")


#Wrangling to visualize it in the tree
group1.2 <- inconsistent2[,c(1,11)]
group2.2 <- inconsistent2[,c(2,12)]

#Substitute column names to merge
colnames(group1.2)[1] <- "Sample"
colnames(group2.2)[1] <- "Sample"

colnames(group1.2)[2] <- "clone_pair"
colnames(group2.2)[2] <- "clone_pair"

#Merge and make unique
clonal_groups2 <- bind_rows(group1.2, group2.2)
clonal_groups2 <-  clonal_groups2 %>% distinct(Sample, .keep_all = TRUE)

clonal_groups2$clone_pair <- as.factor(clonal_groups2$clone_pair)
levels(clonal_groups2$clone_pair)


#Plot


p3 <- tree_plot %<+% clonal_groups2 +
  geom_tiplab(aes(color=factor(clone_pair)),size=0.9, hjust = 1,geom="text")+
  scale_colour_manual(values=colors2, na.value = "black")+
  theme2+
  scale_x_reverse()+
  ggtitle("Inconsistencies highlighted")

p3 

p4 <- p1 + p3

p4

# p5 <- facet_plot(p4, panel = "Genotyped",
#                  data = missing, geom = geom_barh,
#                  mapping = aes(x = missing$Genotyped),
#                  stat='identity', width = 1)
# 
# p5

ggsave(plot = p4, filename = "tree_clonality_inconsistencies.png", dpi = 300, limitsize = FALSE, width = 10, height = 16)


#Extra checking some samples not highlighted with poppr
# Define the search term
search <- "SCCKALS01714"

# Find rows containing the search term in any column
matching_rows <- inconsistent3[apply(inconsistent3, 1, function(row) any(grepl(search, row))), ]

# Print the matching rows
print(matching_rows)


#The inconsistency between the two methods are observed in 10 clonal groups. Clonal groups identified with 
#vcf_clone_detect.py are splitted in two or more groups with poppr. After assessing the trees along with 
#the allelic similarities, that splitting occurs when allelic similarities between samples are close to 99.5%
#(the selected threshold), the branches are not shared or slightly larger than in clonal groups where both 
#methods coincide. Based on the histogram, increasing the threshold does not seem the way to go as the natural
#break in the distribution is at 99.5%. Thus, taking a conservative approch by continuing with the clonal 
#groups identified through  vcf_clone_detect.py but splitting the clonal groups as in poppr on those 10 clonal groups.  

#Select no problematic clonal groups 
non_problematic <- subset(above_99.5, !group1 %in% inconsistent)

#Colapse the samples into a two column dataframe
#Select columns
clonal_pair1 <- non_problematic[,c(1,3,9)]
clonal_pair2 <- non_problematic[,c(2,4,10)]

#Rename columns and merge
colnames(clonal_pair1)[1] <- "Sample"
colnames(clonal_pair2)[1] <- "Sample"

colnames(clonal_pair1)[2] <- "snps"
colnames(clonal_pair2)[2] <- "snps"

colnames(clonal_pair1)[3] <- "group"
colnames(clonal_pair2)[3] <- "group"

clonal_pairs <- bind_rows(clonal_pair1, clonal_pair2)

#Remove duplicated samples
clonal_pairs <-  clonal_pairs %>% distinct(Sample, .keep_all = TRUE)


#Based on number of snps, select the individual with the highest number of snps per group
unique(clonal_pairs$group)
#60

ind_keep1 <- clonal_pairs %>%
  group_by(group) %>%
  filter(snps == max(snps)) %>%
  pull(Sample)

print(ind_keep1)
#60, OK

#Make a list of samples to remove
ind_to_remove1 <- subset(clonal_pairs, !Sample %in% ind_keep1)
#208-60 = 148, OK

##################################
#Repeat same steps for the problematic groups but taking the grouping from poppr
problematic <- subset(above_99.5, group1 %in% inconsistent)

#Select columns, rename and merge
clonal_pair3 <- problematic[,c(1,3,11)]
clonal_pair4 <- problematic[,c(2,4,12)]

colnames(clonal_pair3)[1] <- "Sample"
colnames(clonal_pair4)[1] <- "Sample"

colnames(clonal_pair3)[2] <- "snps"
colnames(clonal_pair4)[2] <- "snps"

colnames(clonal_pair3)[3] <- "group"
colnames(clonal_pair4)[3] <- "group"

clonal_pairs2 <- bind_rows(clonal_pair3, clonal_pair4)
unique(clonal_pairs2$group)
#24

#Remove duplicated samples and select the sample with the hightest number of snps
clonal_pairs2 <-  clonal_pairs2 %>% distinct(Sample, .keep_all = TRUE)

ind_keep2 <- clonal_pairs2 %>%
  group_by(group) %>%
  filter(snps == max(snps))  %>%
  pull(Sample)

print(ind_keep2)
#24, OK

#Make a list with those to remove
ind_to_remove2 <- subset(clonal_pairs2, !Sample %in% ind_keep2)
#93-24 = 69

#Merge the two lists
ind_to_remove <- bind_rows(ind_to_remove1, ind_to_remove2)
#148 + 69 = 217

#Make sure there are no repeated samples
ind_to_remove %>% distinct(Sample, .keep_all = TRUE)
#OK

#Make sure that one of the copies of each technical replicate is included
#Iterate between the technical samples
# Define the search term
search <- "SCCKALSX1676"

# Find rows containing the search term in any column
ind_to_remove[apply(ind_to_remove, 1, function(row) any(grepl(search, row))), ]

#SCCESTSX1418 - Yes
#SCCDIRSX1729 - Yes
#SCCHULSX1792 - Yes, as SCCHULS01792
#SCCESTSX1426 - Yes
#SCCMARSX1890 - Yes, as SCCMARS01890
#SCCESTSX1434 - Yes
#SCCKALSX1676 - Yes

#Visualizing it in the tree
str(ind_to_remove)
ind_to_remove$group <- as.factor(ind_to_remove$group)
levels(ind_to_remove$group)

colors3 <- hue_pal()(82)

removing <- tree_plot %<+% ind_to_remove +
  geom_tiplab(aes(color=factor(group)),size=0.9, geom="text")+
  theme2+
  scale_colour_manual(values=colors3, na.value = "black")

removing

ggsave(plot = removing, filename = "tree_clonality_removing.png", dpi = 300, limitsize = FALSE, width = 10, height = 16)


#Export the list to remove
#samples currenty in the dataset - samples to remove
453-217
#[1] 236

write.table(ind_to_remove[,1], here("clones_to_remove_c4.txt"), quote = FALSE,row.names = FALSE, col.names = FALSE, sep = "\t")

