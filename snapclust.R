#Snapclust
#Alejandra Hernandez


#############################################################
#Dependencies
library("here")
library("vcfR")
library("adegenet")
library("treeio")
library("ggtree")
library("phytools")
library("reshape2")
library("ggplot2")
library("ggstance")
library("grid")


#############################################################
# Set the main working directory:
setwd("/Users/alejhernandez/Dropbox/Rare species/scolymia/dart/D - Overall genetic structure/D1 - SNAPCLUST")



#############################################################
#Open data
#scolymia <- read.vcfR("sco_c4.vcf")
scolymia <- read.vcfR("sco_singlesnp3_c4.vcf")


# Convert to genind object for snapclust
scolymia_genind <- vcfR2genind(scolymia)


# Run snapclust for various K's
# scolymia_snap_k1 <- snapclust(scolymia_genind, 1)
scolymia_snap_k2 <- snapclust(scolymia_genind, 2)
scolymia_snap_k3 <- snapclust(scolymia_genind, 3)
scolymia_snap_k4 <- snapclust(scolymia_genind, 4)
scolymia_snap_k5 <- snapclust(scolymia_genind, 5)
scolymia_snap_k6 <- snapclust(scolymia_genind, 6)
scolymia_snap_k7 <- snapclust(scolymia_genind, 7)
# scolymia_snap_k8 <- snapclust(scolymia_genind, 8)
# scolymia_snap_k9 <- snapclust(scolymia_genind, 9)
# scolymia_snap_k10 <- snapclust(scolymia_genind, 10)
# scolymia_snap_k11 <- snapclust(scolymia_genind, 11)
# scolymia_snap_k12 <- snapclust(scolymia_genind, 12)
# scolymia_snap_k13 <- snapclust(scolymia_genind, 13)
# scolymia_snap_k14 <- snapclust(scolymia_genind, 14)




# Visualize results
# compoplot(scolymia_snap_k1)
compoplot(scolymia_snap_k2)
compoplot(scolymia_snap_k3)
compoplot(scolymia_snap_k4)
compoplot(scolymia_snap_k5)
compoplot(scolymia_snap_k6)
compoplot(scolymia_snap_k7)
# compoplot(scolymia_snap_k8)
# compoplot(scolymia_snap_k9)
# compoplot(scolymia_snap_k10)
# compoplot(scolymia_snap_k11)
# compoplot(scolymia_snap_k12)
# compoplot(scolymia_snap_k13)
# compoplot(scolymia_snap_k14)

# Decide on optimal K
#AIC
# aic <- snapclust.choose.k(max = 14, scolymia_genind)
# plot(1:14, aic, xlab = 'Number of clusters (K)', ylab = 'AIC', type = 'b', pch = 20, cex = 3) # K = 6

aic <- snapclust.choose.k(max = 7, scolymia_genind)
plot(1:7, aic, xlab = 'Number of clusters (K)', ylab = 'AIC', type = 'b', pch = 20, cex = 3) # K = 5

# BIC
# bic <- snapclust.choose.k(max = 14, scolymia_genind, IC = BIC)
# plot(1:14, bic, xlab = 'Number of clusters (K)', ylab = 'BIC', type = 'b', pch = 20, cex = 3) # K = 3

bic <- snapclust.choose.k(max = 7, scolymia_genind, IC = BIC)
plot(1:7, bic, xlab = 'Number of clusters (K)', ylab = 'BIC', type = 'b', pch = 20, cex = 3) # K = 5

# Write K results for use in composite plot
# write.table(scolymia_snap_k1$proba, here("scolymia_snap_k1.txt"), quote = FALSE, col.names = FALSE, sep = "\t")
# write.table(scolymia_snap_k2$proba, here("scolymia_snap_k2.txt"), quote = FALSE, col.names = FALSE, sep = "\t")
# write.table(scolymia_snap_k3$proba, here("scolymia_snap_k3.txt"), quote = FALSE, col.names = FALSE, sep = "\t")
# write.table(scolymia_snap_k4$proba, here("scolymia_snap_k4.txt"), quote = FALSE, col.names = FALSE, sep = "\t")
# write.table(scolymia_snap_k5$proba, here("scolymia_snap_k5.txt"), quote = FALSE, col.names = FALSE, sep = "\t")
# write.table(scolymia_snap_k6$proba, here("scolymia_snap_k6.txt"), quote = FALSE, col.names = FALSE, sep = "\t")
# write.table(scolymia_snap_k7$proba, here("scolymia_snap_k7.txt"), quote = FALSE, col.names = FALSE, sep = "\t")
# write.table(scolymia_snap_k8$proba, here("scolymia_snap_k8.txt"), quote = FALSE, col.names = FALSE, sep = "\t")
# write.table(scolymia_snap_k9$proba, here("scolymia_snap_k9.txt"), quote = FALSE, col.names = FALSE, sep = "\t")
# write.table(scolymia_snap_k10$proba, here("scolymia_snap_k10.txt"), quote = FALSE, col.names = FALSE, sep = "\t")
# write.table(scolymia_snap_k11$proba, here("scolymia_snap_k11.txt"), quote = FALSE, col.names = FALSE, sep = "\t")
# write.table(scolymia_snap_k12$proba, here("scolymia_snap_k12.txt"), quote = FALSE, col.names = FALSE, sep = "\t")
# write.table(scolymia_snap_k13$proba, here("scolymia_snap_k13.txt"), quote = FALSE, col.names = FALSE, sep = "\t")
# write.table(scolymia_snap_k14$proba, here("scolymia_snap_k14.txt"), quote = FALSE, col.names = FALSE, sep = "\t")

write.table(scolymia_snap_k2$proba, here("scolymia_snp_snap_k2.txt"), quote = FALSE, col.names = FALSE, sep = "\t")
write.table(scolymia_snap_k3$proba, here("scolymia_snp_snap_k3.txt"), quote = FALSE, col.names = FALSE, sep = "\t")
write.table(scolymia_snap_k4$proba, here("scolymia_snp_snap_k4.txt"), quote = FALSE, col.names = FALSE, sep = "\t")
write.table(scolymia_snap_k5$proba, here("scolymia_snp_snap_k5.txt"), quote = FALSE, col.names = FALSE, sep = "\t")
write.table(scolymia_snap_k6$proba, here("scolymia_snp_snap_k6.txt"), quote = FALSE, col.names = FALSE, sep = "\t")
write.table(scolymia_snap_k7$proba, here("scolymia_snp_snap_k7.txt"), quote = FALSE, col.names = FALSE, sep = "\t")


###########################################################################################################


#load in tree
tree <- read.nexus("/Users/alejhernandez/Dropbox/Rare species/scolymia/dart/C - Filtering and QC of loci and SNPs/C4 - Detect and eliminate clones/sco_c4.tre")
tree <- midpoint.root(tree)

#add snapclust data
k2 <- as.data.frame(scolymia_snap_k2$proba)
k3 <- as.data.frame(scolymia_snap_k3$proba)
k4 <- as.data.frame(scolymia_snap_k4$proba)
k5 <- as.data.frame(scolymia_snap_k5$proba)
k6 <- as.data.frame(scolymia_snap_k6$proba)
k7 <- as.data.frame(scolymia_snap_k7$proba)
# k8 <- as.data.frame(scolymia_snap_k8$proba)
# k9 <- as.data.frame(scolymia_snap_k9$proba)
# k10 <- as.data.frame(scolymia_snap_k10$proba)
# k11 <- as.data.frame(scolymia_snap_k11$proba)
# k12 <- as.data.frame(scolymia_snap_k12$proba)

k2$INDV <- row.names(k2)
k3$INDV <- row.names(k3)
k4$INDV <- row.names(k4)
k5$INDV <- row.names(k5)
k6$INDV <- row.names(k6)
k7$INDV <- row.names(k7)
# k8$INDV <- row.names(k8)
# k9$INDV <- row.names(k9)
# k10$INDV <- row.names(k10)
# k11$INDV <- row.names(k11)
# k12$INDV <- row.names(k12)

k2 <- k2[,c(3,1,2)]
k3 <- k3[,c(4,1:3)]
k4 <- k4[,c(5,1:4)]
k5 <- k5[,c(6,1:5)]
k6 <- k6[,c(7,1:6)]
k7 <- k7[,c(8,1:7)]
# k8 <- k8[,c(9,1:8)]
# k9 <- k9[,c(10,1:9)]
# k10 <- k10[,c(11,1:10)]
# k11 <- k11[,c(12,1:11)]
# k12 <- k12[,c(14,1:13)]

colnames(k2) <- c("INDV", "pop1", "pop2")
colnames(k3) <- c("INDV", "pop1", "pop2","pop3")
colnames(k4) <- c("INDV", "pop1", "pop2","pop3","pop4")
colnames(k5) <- c("INDV", "pop1", "pop2","pop3","pop4","pop5")
colnames(k6) <- c("INDV", "pop1", "pop2","pop3","pop4","pop5","pop6")
colnames(k7) <- c("INDV", "pop1", "pop2","pop3","pop4","pop5","pop6","pop7")
# colnames(k8) <- c("INDV", "pop1", "pop2","pop3","pop4","pop5","pop6","pop7","pop8")
# colnames(k9) <- c("INDV", "pop1", "pop2","pop3","pop4","pop5","pop6","pop7","pop8","pop9")
# colnames(k10) <- c("INDV", "pop1", "pop2","pop3","pop4","pop5","pop6","pop7","pop8","pop9","pop10")
# colnames(k11) <- c("INDV", "pop1", "pop2","pop3","pop4","pop5","pop6","pop7","pop8","pop9","pop10","pop11")
# colnames(k12) <- c("INDV", "pop1", "pop2","pop3","pop4","pop5","pop6","pop7","pop8","pop9","pop10","pop11","pop12","pop13")

d_k2 <- data.frame(INDV=tree$tip.label)
d_k3 <- data.frame(INDV=tree$tip.label)
d_k4 <- data.frame(INDV=tree$tip.label)
d_k5 <- data.frame(INDV=tree$tip.label)
d_k6 <- data.frame(INDV=tree$tip.label)
d_k7 <- data.frame(INDV=tree$tip.label)
# d_k8 <- data.frame(INDV=tree$tip.label)
# d_k9 <- data.frame(INDV=tree$tip.label)
# d_k10 <- data.frame(INDV=tree$tip.label)
# d_k11 <- data.frame(INDV=tree$tip.label)
# d_k12 <- data.frame(INDV=tree$tip.label)

d_k2$p1 <- with(k2, pop1 [match(d_k2$INDV, INDV)])
d_k2$p2 <- with(k2, pop2 [match(d_k2$INDV, INDV)])

d_k3$p1 <- with(k3, pop1 [match(d_k3$INDV, INDV)])
d_k3$p2 <- with(k3, pop2 [match(d_k3$INDV, INDV)])
d_k3$p3 <- with(k3, pop3 [match(d_k3$INDV, INDV)])

d_k4$p1 <- with(k4, pop1 [match(d_k4$INDV, INDV)])
d_k4$p2 <- with(k4, pop2 [match(d_k4$INDV, INDV)])
d_k4$p3 <- with(k4, pop3 [match(d_k4$INDV, INDV)])
d_k4$p4 <- with(k4, pop4 [match(d_k4$INDV, INDV)])

d_k5$p1 <- with(k5, pop1 [match(d_k5$INDV, INDV)])
d_k5$p2 <- with(k5, pop2 [match(d_k5$INDV, INDV)])
d_k5$p3 <- with(k5, pop3 [match(d_k5$INDV, INDV)])
d_k5$p4 <- with(k5, pop4 [match(d_k5$INDV, INDV)])
d_k5$p5 <- with(k5, pop5 [match(d_k5$INDV, INDV)])

d_k6$p1 <- with(k6, pop1 [match(d_k6$INDV, INDV)])
d_k6$p2 <- with(k6, pop2 [match(d_k6$INDV, INDV)])
d_k6$p3 <- with(k6, pop3 [match(d_k6$INDV, INDV)])
d_k6$p4 <- with(k6, pop4 [match(d_k6$INDV, INDV)])
d_k6$p5 <- with(k6, pop5 [match(d_k6$INDV, INDV)])
d_k6$p6 <- with(k6, pop6 [match(d_k6$INDV, INDV)])

d_k7$p1 <- with(k7, pop1 [match(d_k7$INDV, INDV)])
d_k7$p2 <- with(k7, pop2 [match(d_k7$INDV, INDV)])
d_k7$p3 <- with(k7, pop3 [match(d_k7$INDV, INDV)])
d_k7$p4 <- with(k7, pop4 [match(d_k7$INDV, INDV)])
d_k7$p5 <- with(k7, pop5 [match(d_k7$INDV, INDV)])
d_k7$p6 <- with(k7, pop6 [match(d_k7$INDV, INDV)])
d_k7$p7 <- with(k7, pop7 [match(d_k7$INDV, INDV)])
# 
# d_k8$p1 <- with(k8, pop1 [match(d_k8$INDV, INDV)])
# d_k8$p2 <- with(k8, pop2 [match(d_k8$INDV, INDV)])
# d_k8$p3 <- with(k8, pop3 [match(d_k8$INDV, INDV)])
# d_k8$p4 <- with(k8, pop4 [match(d_k8$INDV, INDV)])
# d_k8$p5 <- with(k8, pop5 [match(d_k8$INDV, INDV)])
# d_k8$p6 <- with(k8, pop6 [match(d_k8$INDV, INDV)])
# d_k8$p7 <- with(k8, pop7 [match(d_k8$INDV, INDV)])
# d_k8$p8 <- with(k8, pop8 [match(d_k8$INDV, INDV)])
# 
# d_k9$p1 <- with(k9, pop1 [match(d_k9$INDV, INDV)])
# d_k9$p2 <- with(k9, pop2 [match(d_k9$INDV, INDV)])
# d_k9$p3 <- with(k9, pop3 [match(d_k9$INDV, INDV)])
# d_k9$p4 <- with(k9, pop4 [match(d_k9$INDV, INDV)])
# d_k9$p5 <- with(k9, pop5 [match(d_k9$INDV, INDV)])
# d_k9$p6 <- with(k9, pop6 [match(d_k9$INDV, INDV)])
# d_k9$p7 <- with(k9, pop7 [match(d_k9$INDV, INDV)])
# d_k9$p8 <- with(k9, pop8 [match(d_k9$INDV, INDV)])
# d_k9$p9 <- with(k9, pop9 [match(d_k9$INDV, INDV)])
# 
# d_k10$p1 <- with(k10, pop1 [match(d_k10$INDV, INDV)])
# d_k10$p2 <- with(k10, pop2 [match(d_k10$INDV, INDV)])
# d_k10$p3 <- with(k10, pop3 [match(d_k10$INDV, INDV)])
# d_k10$p4 <- with(k10, pop4 [match(d_k10$INDV, INDV)])
# d_k10$p5 <- with(k10, pop5 [match(d_k10$INDV, INDV)])
# d_k10$p6 <- with(k10, pop6 [match(d_k10$INDV, INDV)])
# d_k10$p7 <- with(k10, pop7 [match(d_k10$INDV, INDV)])
# d_k10$p8 <- with(k10, pop8 [match(d_k10$INDV, INDV)])
# d_k10$p9 <- with(k10, pop9 [match(d_k10$INDV, INDV)])
# d_k10$p10 <- with(k10, pop10 [match(d_k10$INDV, INDV)])
# 
# d_k11$p1 <- with(k11, pop1 [match(d_k11$INDV, INDV)])
# d_k11$p2 <- with(k11, pop2 [match(d_k11$INDV, INDV)])
# d_k11$p3 <- with(k11, pop3 [match(d_k11$INDV, INDV)])
# d_k11$p4 <- with(k11, pop4 [match(d_k11$INDV, INDV)])
# d_k11$p5 <- with(k11, pop5 [match(d_k11$INDV, INDV)])
# d_k11$p6 <- with(k11, pop6 [match(d_k11$INDV, INDV)])
# d_k11$p7 <- with(k11, pop7 [match(d_k11$INDV, INDV)])
# d_k11$p8 <- with(k11, pop8 [match(d_k11$INDV, INDV)])
# d_k11$p9 <- with(k11, pop9 [match(d_k11$INDV, INDV)])
# d_k11$p10 <- with(k11, pop10 [match(d_k11$INDV, INDV)])
# d_k11$p11 <- with(k11, pop11 [match(d_k11$INDV, INDV)])
# 
# d_k12$p1 <- with(k12, pop1 [match(d_k12$INDV, INDV)])
# d_k12$p2 <- with(k12, pop2 [match(d_k12$INDV, INDV)])
# d_k12$p3 <- with(k12, pop3 [match(d_k12$INDV, INDV)])
# d_k12$p4 <- with(k12, pop4 [match(d_k12$INDV, INDV)])
# d_k12$p5 <- with(k12, pop5 [match(d_k12$INDV, INDV)])
# d_k12$p6 <- with(k12, pop6 [match(d_k12$INDV, INDV)])
# d_k12$p7 <- with(k12, pop7 [match(d_k12$INDV, INDV)])
# d_k12$p8 <- with(k12, pop8 [match(d_k12$INDV, INDV)])
# d_k12$p9 <- with(k12, pop9 [match(d_k12$INDV, INDV)])
# d_k12$p10 <- with(k12, pop10 [match(d_k12$INDV, INDV)])
# d_k12$p11 <- with(k12, pop11 [match(d_k12$INDV, INDV)])
# d_k12$p12 <- with(k12, pop12 [match(d_k12$INDV, INDV)])
# d_k12$p13 <- with(k12, pop13 [match(d_k12$INDV, INDV)])



#Preparing data for the graph
k2_melt <- melt(d_k2, id.vars = "INDV")
k3_melt <- melt(d_k3, id.vars = "INDV")
k4_melt <- melt(d_k4, id.vars = "INDV")
k5_melt <- melt(d_k5, id.vars = "INDV")
k6_melt <- melt(d_k6, id.vars = "INDV")
k7_melt <- melt(d_k7, id.vars = "INDV")
# k8_melt <- melt(d_k8, id.vars = "INDV")
# k9_melt <- melt(d_k9, id.vars = "INDV")
# k10_melt <- melt(d_k10, id.vars = "INDV")
# k11_melt <- melt(d_k11, id.vars = "INDV")
# k12_melt <- melt(d_k12, id.vars = "INDV")





#Plot
colors <- c("#001219", "#005f73","#0a9396", "#94d2bd", "#e9d8a6","#ee9b00", "#ca6702", "#bb3e03",
            "#ae2012", "#9b2226", "#6a040f", "#03071e","#e3d5ca")



tree_plot <- ggtree(tree, layout = "rectangular", size=0.15, open.angle = 5,tree_width = 3) + #geom_tippoint(size = 0.1) + 
  geom_tiplab(align = TRUE, size=0.5,linetype = "dotted", linesize = 0.1) +
  theme_tree2()+
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))
tree_plot


p1 <- tree_plot + geom_facet(panel = "k=2", 
                             data = k2_melt, geom = geom_barh, 
                             mapping = aes(x = value, fill = variable), 
                             stat='identity', width = 1) + scale_fill_manual(values = colors)


p2 <- p1 + geom_facet(panel = "k=3", 
                      data = k3_melt, geom = geom_barh, 
                      mapping = aes(x = value, fill = variable), 
                      stat='identity', width = 1) 

p3 <- p2 + geom_facet(panel = "k=4", 
                      data = k4_melt, geom = geom_barh, 
                      mapping = aes(x = value, fill = variable), 
                      stat='identity', width = 1) 


p4 <- p3 + geom_facet(panel = "k=5", 
                      data = k5_melt, geom = geom_barh, 
                      mapping = aes(x = value, fill = variable), 
                      stat='identity', width = 1) 

p5 <- p4 + geom_facet(panel = "k=6", 
                      data = k6_melt, geom = geom_barh, 
                      mapping = aes(x = value, fill = variable), 
                      stat='identity', width =1) 

p6 <- p5 + geom_facet(panel = "k=7",
                      data = k7_melt, geom = geom_barh,
                      mapping = aes(x = value, fill = variable),
                      stat='identity', width = 1)+ theme(legend.position = "none")
# 
# p7 <- p6 + geom_facet(panel = "k=8", 
#                       data = k8_melt, geom = geom_barh, 
#                       mapping = aes(x = value, fill = variable), 
#                       stat='identity', width = 1) + theme(legend.position = "none")
# 
# p8 <- p7 + geom_facet(panel = "k=9", 
#                       data = k9_melt, geom = geom_barh, 
#                       mapping = aes(x = value, fill = variable), 
#                       stat='identity', width = 1) 
# 
# p9 <- p8 + geom_facet(panel = "k=10", 
#                       data = k10_melt, geom = geom_barh, 
#                       mapping = aes(x = value, fill = variable), 
#                       stat='identity', width = 1) 
# 
# 
# p10 <- p9 + geom_facet(panel = "k=11", 
#                       data = k11_melt, geom = geom_barh, 
#                       mapping = aes(x = value, fill = variable), 
#                       stat='identity', width = 1)
# 
# p11 <- p10 + geom_facet(panel = "k=12", 
#                        data = k12_melt, geom = geom_barh, 
#                        mapping = aes(x = value, fill = variable), 
#                        stat='identity', width = 1) + theme(legend.position = "none")
# 
# p11

p6

ggsave(p6, file="sco_SNAPCLUST_singlesnp.png", height = 8, width = 14, dpi = 300)

