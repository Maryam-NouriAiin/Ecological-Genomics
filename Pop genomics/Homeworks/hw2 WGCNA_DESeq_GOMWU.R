## Set your working directory
setwd("/Users/maryamnouri-aiin/Desktop/Ecological-Genomics/Pop genomics/results/Transcriptome/")

# Load the package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GO.db")

BiocManager::install("impute")
###WGCNA weighted gene co-expression network analysis
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

library(DESeq2)
library(ggplot2)
library(lattice)
library(tidyverse)
library(ggpubr)
library(ggrepel)# install.packages("remotes")
# remotes::install_github("kevinblighe/CorLevelPlot")
remotes::install_github("kevinblighe/CorLevelPlot")
library(CorLevelPlot) 
library(gridExtra)

library(Rmisc) 


###################################################################################### 

# 1. Import the counts matrix and metadata and filter using DESeq2,,,, unfliter was too much unmanageable data set


countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample description table
# conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
# head(conds)

sample_metadata = read.table(file = "Ahud_trait_data.txt",header=T, row.names = 1)
### design as 1 to keep WGCNA agnostic
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=sample_metadata, 
                              design= ~ 1)

dim(dds)

# Filter out genes with too few reads - remove all genes with counts < 15 in more than 75% of samples, so ~28)
## suggested by WGCNA on RNAseq FAQ

dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 
# [1] 25260, that have at least 15 reads (a.k.a counts) in 75% of the samples

# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)

###############################################################################
# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(countsTable))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)
### if in table we have any false we need to filter them (The number here is 38 which the number of our samples, it means all samples are there and non of them are ID as outliyers even if two of them look a little off, those two are the one with the lowest counts)

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(countsTable)), method = "average")
plot(htree) 



# pca - method 2

pca <- prcomp(t(countsTable))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

###############################################################################

# 3. Normalization using DESeq----------------------------------------------------------------------

colData <- row.names(sample_metadata)

# making the rownames and column names identical
# make sure everything is in the order
all(rownames(colData) %in% colnames(countsTableRound)) # to see if all samples are present in both
all(rownames(colData) == colnames(countsTableRound))  # to see if all samples are in the same order



# perform variance stabilization
dds_norm <- vst(dds)
# dds_norm <- vst(normalized_counts)

# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()

#################################################################################


# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers,, generating set of things to run which are called power
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

power6 <- c(c(1:10), seq(from = 12, to = 50, by = 6))

power12 <- c(c(1:10), seq(from = 12, to = 50, by = 12))

# Call the network topology analysis function; this step takes a couple minutes
#signed means if it is postively regulated or not, if not signed is just an absolute value or genes that uncorrelated
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft6 <- pickSoftThreshold(norm.counts,
                         powerVector = power6,
                         networkType = "signed",
                         verbose = 5)

sft12 <- pickSoftThreshold(norm.counts,
                          powerVector = power12,
                          networkType = "signed",
                          verbose = 5)

sft.data <- sft$fitIndices
sft.data6 <- sft6$fitIndices
sft.data12 <- sft12$fitIndices
######################### Extract the necessary data from sft.data



# visualization to pick power

# a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
#   geom_point() +
#   geom_text(nudge_y = 0.1) +
#   geom_hline(yintercept = 0.8, color = 'red') +
#   labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
#   theme_classic()

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point(size = 2, shape = 16, color = "midnightblue") +
  geom_text(nudge_y = 0.1, color = "black") +
  geom_hline(yintercept = 0.8, color = "orangered4", size = 1.2, linetype = "dashed") +
  labs(x = "Power", y = "Scale free topology model fit, signed R^2", title = "") +
  theme_minimal()

# Display the plot
a1

a2 <- ggplot(sft.data6, aes(Power, SFT.R.sq, label = Power)) +
  geom_point(size = 2, shape = 16, color = "midnightblue") +
  geom_text(nudge_y = 0.1, color = "black") +
  geom_hline(yintercept = 0.8, color = "orangered4", size = 1.2, linetype = "dashed") +
  labs(x = "Power", y = "Scale free topology model fit, signed R^2", title = "") +
  theme_minimal()

a3 <- ggplot(sft.data12, aes(Power, SFT.R.sq, label = Power)) +
  geom_point(size = 2, shape = 16, color = "midnightblue") +
  geom_text(nudge_y = 0.1, color = "black") +
  geom_hline(yintercept = 0.8, color = "orangered4", size = 1.2, linetype = "dashed") +
  labs(x = "Power", y = "Scale free topology model fit, signed R^2", title = "") +
  theme_minimal()
# a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
#   geom_point(size = 3, shape = 16, color = "midnightblue", fill = "white", stroke = 1) +
#   geom_text_repel(
#     nudge_y = 0.1,
#     color = "black",
#     segment.color = "gray",
#     segment.size = 0.5,
#     box.padding = 0.4,
#     force = 10,
#     min.segment.length = 0.1,
#     max.segments = 20
#   ) +
#   geom_hline(yintercept = 0.8, color = "orangered4", size = 1.2, linetype = "dashed") +
#   labs(x = "Power", y = "Scale free topology model fit, signed R^2", title = "Power vs. Scale Free Topology") +
#   theme_minimal() +
#   theme(
#     text = element_text(size = 10),
#     axis.text = element_text(size = 10),
#     axis.title = element_text(size = 12, face = "bold"),
#     plot.title = element_text(size = 14, face = "bold")
#   )
# 
# a1


# 
a4 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point(size = 2, shape = 16, color = "midnightblue")  +
  geom_text_repel(nudge_y = 0.1, color = "black")  +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_minimal()

a5 <- ggplot(sft.data6, aes(Power, mean.k., label = Power)) +
  geom_point(size = 2, shape = 16, color = "midnightblue")  +
  geom_text_repel(nudge_y = 0.1, color = "black")  +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_minimal()

a6 <- ggplot(sft.data12, aes(Power, mean.k., label = Power)) +
  geom_point(size = 2, shape = 16, color = "midnightblue")  +
  geom_text_repel(nudge_y = 0.1, color = "black")  +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_minimal()

# a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
#   geom_point(size = 3, shape = 16, color = "midnightblue", fill = "white", stroke = 1) +
#   geom_text_repel(
#     nudge_y = 0.1,
#     color = "black",
#     segment.color = "gray",
#     segment.size = 0.5,
#     box.padding = 0.4,
#     force = 10,
#     min.segment.length = 0.1,
#     max.segments = 20
#   ) +
#   labs(x = 'Power', y = 'Mean Connectivity') +
#   theme_minimal() +
#   theme(
#     text = element_text(size = 10),
#     axis.text = element_text(size = 10),
#     axis.title = element_text(size = 12, face = "bold")
#   )

a2


grid.arrange(a1, a2, nrow = 1)
# based on this plot, choose a soft power to maximize R^2 (above 0.8) and minimize connectivity
# for these ahud data: 6-8; Higher R2 should yield more modules.
tiff('Power', units="in", width=18, height=11, res=300)
ggarrange(
  a1,
  a4,
  a2,
  a5,
  a3,
  a6,
  nrow = 3,
  ncol = 2,
  labels = c("A", "B", "C","D","E", "F"
  )
)

dev.off()


tiff('Power6', units="in", width=14, height=8, res=300)

ggarrange(
  a3,
  a4,
  nrow = 2,
  ncol = 1,
  labels = c("C", "D"
  )
)

dev.off()



# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power2<- 2
soft_power6<- 6
soft_power12<- 12 #should try running with 7 based on R^2
temp_cor <- cor
cor <- WGCNA::cor # use the 'cor' function from the WGCNA package


# this step also takes a few minutes; ideally your maxBlockSize is larger than your number of genes to run the memory-intensive network construction all at once.
bwnet2<- blockwiseModules(norm.counts,
                          maxBlockSize = 26000,
                          minModuleSize = 30, 
                          reassignThreshold=0,
                          TOMType = "signed",
                          power = soft_power2,
                          mergeCutHeight = 0.25,
                          numericLabels = F,
                          randomSeed = 1234,
                          verbose = 3)

bwnet6<- blockwiseModules(norm.counts,
                          maxBlockSize = 26000,
                          minModuleSize = 30, 
                          reassignThreshold=0,
                          TOMType = "signed",
                          power = soft_power6,
                          mergeCutHeight = 0.25,
                          numericLabels = F,
                          randomSeed = 1234,
                          verbose = 3)

bwnet12<- blockwiseModules(norm.counts,
                          maxBlockSize = 26000,
                          minModuleSize = 30, 
                          reassignThreshold=0,
                          TOMType = "signed",
                          power = soft_power12,
                          mergeCutHeight = 0.25,
                          numericLabels = F,
                          randomSeed = 1234,
                          verbose = 3)

# TOMtype (Topological Overlap Matrix type) parameter - unsigned - doesn't consider positive/negative co-expression
# signed - when you want to consider the direction of co-expression interaction, e.g., activating or inhibiting
# WGCNA often uses a dendrogram-based approach to identify modules. The choice of the 
# height cut in the dendrogram can determine the number of modules. Selecting a higher
# cut height results in fewer, larger modules, while a lower cut height leads to more, 
# smaller modules.
# What I did to save it:
saveRDS(bwnet2, file = "bwnet2.rds")

# To load the object
bwnet2 <- readRDS("bwnet2.rds")
cor <- temp_cor

###############################################################################################
# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes2 <- bwnet2$MEs

head(module_eigengenes2)


# get number of genes for each module### higher number higher expression
table(bwnet2$colors)
#Power 6= 11 modules:     black      blue     brown     green      grey   magenta      pink    purple       red turquoise    yellow 
#233      7919      3614       367      1588       144       164        91       331      9806      1003 

# What I did to save it:
saveRDS(bwnet6, file = "bwnet6.rds")

# To load the object
bwnet6 <- readRDS("bwnet6.rds")
cor <- temp_cor

###############################################################################################
# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes6 <- bwnet6$MEs

head(module_eigengenes6)


# get number of genes for each module### higher number higher expression
table(bwnet6$colors)

# What I did to save it:
saveRDS(bwnet12, file = "bwnet12.rds")

# To load the object
bwnet12 <- readRDS("bwnet12.rds")
cor <- temp_cor

###############################################################################################
# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes12 <- bwnet12$MEs

head(module_eigengenes12)


# get number of genes for each module### higher number higher expression
table(bwnet12$colors)

# Plot the dendrogram and the module colors before and after merging underneath
# high eigengenes means more expression,,,, the colores are the modules
dendro2 <- plotDendroAndColors(bwnet2$dendrograms[[1]], cbind(bwnet2$unmergedColors, bwnet2$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05,
                    main="")

dendro6 <- plotDendroAndColors(bwnet6$dendrograms[[1]], cbind(bwnet6$unmergedColors, bwnet6$colors),
                               c("unmerged", "merged"),
                               dendroLabels = FALSE,
                               addGuide = TRUE,
                               hang= 0.03,
                               guideHang = 0.05,
                               main="")
tiff('Module yellow 12', units="in", width=14, height=8, res=300)
plotDendroAndColors(bwnet12$dendrograms[[1]], cbind(bwnet12$unmergedColors, bwnet12$colors),
                               c("unmerged", "merged"),
                               dendroLabels = FALSE,
                               addGuide = TRUE,
                               hang= 0.03,
                               guideHang = 0.05,
                               main="")
dev.off()




# grey module = all genes that doesn't fall into other modules were assigned to the grey module
# with higher soft power, more genes fall into the grey module

#####################################################################
# 6A. Relate modules to traits --------------------------------------------------
# module trait associations
#### c(5,8,11,14,17)] cols that we want fitness, EPR egg production rate, HS hatching success, survival, Devlopment time
traits <- sample_metadata[, c(5,8,11,14,17)]
#view(sample_metadata)

# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

### P for pierson correlation,,, and cor pvalue
module.trait.corr2 <- cor(module_eigengenes2, traits, use = 'p')
module.trait.corr.pvals2 <- corPvalueStudent(module.trait.corr2, nSamples)

module.trait.corr6 <- cor(module_eigengenes6, traits, use = 'p')
module.trait.corr.pvals6 <- corPvalueStudent(module.trait.corr6, nSamples)

module.trait.corr12 <- cor(module_eigengenes12, traits, use = 'p')
module.trait.corr.pvals12 <- corPvalueStudent(module.trait.corr12, nSamples)

# visualize module-trait association as a heatmap

heatmap.data2 <- merge(module_eigengenes2, traits, by = 'row.names')

head(heatmap.data2)

heatmap.data2 <- heatmap.data2 %>% 
  column_to_rownames(var = 'Row.names')


names(heatmap.data2)

###

heatmap.data6 <- merge(module_eigengenes6, traits, by = 'row.names')

head(heatmap.data6)

heatmap.data6 <- heatmap.data6 %>% 
  column_to_rownames(var = 'Row.names')


names(heatmap.data6)

###

heatmap.data12 <- merge(module_eigengenes12, traits, by = 'row.names')

head(heatmap.data12)

heatmap.data12 <- heatmap.data12 %>% 
  column_to_rownames(var = 'Row.names')


names(heatmap.data12)

### red positive correlation, blue no correlation
#Mean fitness is correlated with grey module, no sig for hatching success, 2 for survival, one for EPR, Dev time with red and blue module
tiff('Heatmap power12', units="in", width=7, height=6, res=300)

CorLevelPlot(heatmap.data12,
             x = names(heatmap.data2)[11:15],
             y = names(heatmap.data2)[1:10],
             col = c("white", "lightyellow", "lightblue", "coral", "darkred"))

dev.off()

module.gene.mapping12 <- as.data.frame(bwnet12$colors) # assigns module membership to each gene
module.gene.mapping12 %>% 
  filter(`bwnet12$colors` == 'yellow') %>% 
  rownames()
### pulling the number of generation
groups <- sample_metadata[,c(3,1)]
module_eigengene.metadata12 <- merge(groups, heatmap.data12, by = 'row.names')

#Create a summary data frame of a particular module eigengene information
MEyellow_summary12 <- summarySE(module_eigengene.metadata12, measurevar="MEyellow", groupvars=c("Generation","treatment"))

#Plot a line interaction plot of a particular module eigengene


#tiff('Module yellow', units="in", width=7, height=6, res=300)

MY1 <- ggplot(MEyellow_summary12, aes(x=as.factor(Generation), y=MEyellow, color=treatment, fill = treatment, shape = treatment)) +
  geom_point(size=5, stroke = 1.5 ) +
  geom_errorbar(aes(ymin=MEyellow-se, ymax=MEyellow+se), width=.15) +
  geom_line(aes(color = treatment, group = treatment, linetype = treatment)) +
  scale_fill_manual(values = c('#440154FF', '#404688FF', '#27AD81FF', 'darkred'), labels = c("Ambient", "Acidification", "Warming", "OWA")) +
  scale_color_manual(values = c('#440154FF', '#404688FF', '#27AD81FF', 'darkred')) +
  xlab("Generation") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "lightgrey", fill = NA, size = 1),
    text = element_text(size = 12),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(0, 6, 0, 6)
  ) +
  guides(fill = "none", linetype = "none")


#dev.off()
###y axis is eigen, from pos to neg decrease in expression over generation
####################################################################################
# 6B. Intramodular analysis: Identifying driver genes ---------------

# Get top hub genes (genes with highest connectivity in the network)
# pull the genes that have the highest connectivity
hubs  <-  chooseTopHubInEachModule(norm.counts, bwnet12$colors, type = "signed", omitColors = "")
hubs
# black 
# "TRINITY_DN24998_c0_g1::TRINITY_DN24998_c0_g1_i3::g.49112::m.49112" 
# blue 
# "TRINITY_DN6125_c0_g1::TRINITY_DN6125_c0_g1_i1::g.23836::m.23836" 
# brown 
# "TRINITY_DN1784_c0_g1::TRINITY_DN1784_c0_g1_i3::g.9203::m.9203" 
# green 
# "TRINITY_DN239_c0_g1::TRINITY_DN239_c0_g1_i12::g.1858::m.1858" 
# grey 
# "TRINITY_DN13235_c0_g2::TRINITY_DN13235_c0_g2_i3::g.38482::m.38482" 
# magenta 
# "TRINITY_DN4230_c0_g2::TRINITY_DN4230_c0_g2_i1::g.18061::m.18061" 
# pink 
# "TRINITY_DN865_c0_g1::TRINITY_DN865_c0_g1_i27::g.5072::m.5072" 
# purple 
# "TRINITY_DN1336_c0_g1::TRINITY_DN1336_c0_g1_i4::g.7523::m.7523" 
# red 
# "TRINITY_DN3251_c0_g1::TRINITY_DN3251_c0_g1_i25::g.14716::m.14716" 
# turquoise 
# "TRINITY_DN2215_c0_g1::TRINITY_DN2215_c0_g1_i1::g.11114::m.11114" 
# yellow 
# "TRINITY_DN11845_c0_g1::TRINITY_DN11845_c0_g1_i9::g.36434::m.36434" 
### Plot Individual genes  to check! ### 
## pick the genes from previouse code: TRINITY_DN11845_c0_g1::TRINITY_DN11845_c0_g1_i9::g.36434::m.36434
d <-plotCounts(dds, gene="TRINITY_DN63380_c0_g1::TRINITY_DN63380_c0_g1_i1::g.63074::m.63074", intgroup = (c("treatment","Generation")), returnData=TRUE)
d_summary <- summarySE(d, measurevar = "count", groupvars=c("Generation","treatment"))
#tiff('Module yellow_gene', units="in", width=10, height=5, res=300)

MY2 <- ggplot(d_summary, aes(x = as.factor(Generation), y = count, color = treatment, fill = treatment, shape = treatment)) +
  geom_point(size = 5, stroke = 1.5) +
  geom_errorbar(aes(ymin = count - se, ymax = count + se), width = 0.15) +
  geom_line(aes(color = treatment, group = treatment, linetype = treatment)) +
  scale_fill_manual(values = c('#440154FF', '#404688FF', '#27AD81FF', 'darkred'), labels = c("Ambient", "Acidification", "Warming", "OWA")) +
  scale_color_manual(values = c('#440154FF', '#404688FF', '#27AD81FF', 'darkred')) +
  xlab("Generation") +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.border = element_rect(color = "lightgrey", fill = NA, size = 1),
    text = element_text(size = 12),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(0, 6, 0, 6)
  ) +
  guides(fill = "none", linetype = "none")

tiff('Module yellow 6', units="in", width=14, height=8, res=300)

ggarrange(
MY1,
MY2,
  nrow = 1,
  ncol = 2,
  labels = c("C", "D")
)

dev.off()
#dev.off()

#### extra codes for the plot manipulation
# scale_color_manual(values = c('#6699CC',"#F2AD00","#00A08A", "#CC3333")) +
#   scale_shape_manual(values=c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "OWA"))+
#   scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"), labels = c("Ambient", "Acidification","Warming", "OWA"))+
### values of the eigen genes are pos for high expressiona and neg for los expression


# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)


module.membership.measure.pvals[1:10,1:10]

###########################################################################################
# Make a heat map of gene expressions within modules.
# Use the norm.counts matrix, subset based on module membership

library(pheatmap)
t_norm.counts <- norm.counts %>% t() %>% as.data.frame()

# purple module
purple_transcripts <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'purple') %>% 
  rownames()

t_norm.counts_purple <- t_norm.counts %>% 
  filter(row.names(t_norm.counts) %in% purple_transcripts)

t_norm.counts_purple <- t_norm.counts_purple - rowMeans(t_norm.counts_purple)
df <- as.data.frame(colData(dds)[,c("generation","treatment")])

#blue to yellow color scheme
paletteLength <- 50
myColor <- colorRampPalette(c("dodgerblue", "white", "yellow"))(paletteLength)
myBreaks <- c(seq(min(t_norm.counts_purple), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(t_norm.counts_purple)/paletteLength, max(t_norm.counts_purple), length.out=floor(paletteLength/2)))
pheatmap(t_norm.counts_purple, color = myColor, breaks = myBreaks,
         show_colnames = FALSE, show_rownames = FALSE, annotation_col = df, main = "Yellow")

#### heatmap shows 91 genes in purple module,,, could show low expression in this module, beside the top right corner that shows a little higher expression

# Make a heat map of gene expressions within modules.
# Use the norm.counts matrix, subset based on module membership
t_norm.counts <- norm.counts %>% t() %>% as.data.frame()

# Yellow module
yellow_transcripts <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'yellow') %>% 
  rownames()

t_norm.counts_yellow <- t_norm.counts %>% 
  filter(row.names(t_norm.counts) %in% yellow_transcripts)

t_norm.counts_yellow <- t_norm.counts_yellow - rowMeans(t_norm.counts_yellow)
df <- as.data.frame(colData(dds)[,c("eneration","treatment")])

#blue to yellow color scheme
paletteLength <- 50
myColor <- colorRampPalette(c("dodgerblue", "black", "yellow"))(paletteLength)
myBreaks <- c(seq(min(t_norm.counts_yellow), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(t_norm.counts_yellow)/paletteLength, max(t_norm.counts_yellow), length.out=floor(paletteLength/2)))
pheatmap(t_norm.counts_yellow, color = myColor, breaks = myBreaks,
         show_colnames = FALSE, show_rownames = FALSE, annotation_col = df, main = "Yellow")

##########################
##########################
##########################
#####GO
##########################
##########################
##########################

#FET -cutoff Fisher exact method
## Set your working directory
setwd("/Users/maryamnouri-aiin/Desktop/Ecological-Genomics/Pop genomics/results/Transcriptome/")

## Import the libraries that we're likely to need in this session
library(DESeq2)

# Try with new counts table from filtered transcriptome assembly
countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)


head(countsTable)
dim(countsTable)
#[1] 130580     38 - genes
# [1] 349516     38 - isoforms
# [1] 67916    38 - filtered assembly

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample discription table
conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)

#################### MODEL NUMBER 2 - subset to focus on effect of treatment for each generation

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ treatment)

dim(dds)
# [1] 130580     38

# Filter ----What we set up 30 reads, 28 of the sample which is 75% of them filtering the reads that do not have depth
dds <- dds[rowSums(counts(dds) >= 30) >= 28,]
nrow(dds) 

# Subset the DESeqDataSet to the specific level of the "generation" factor
dds_F0 <- subset(dds, select = generation == 'F0')
dim(dds_F0)

# Perform DESeq2 analysis on the subset
dds_F0 <- DESeq(dds_F0)

resultsNames(dds_F0)
# [1] "Intercept"           "treatment_OA_vs_AM"  "treatment_OW_vs_AM"  "treatment_OWA_vs_AM"

###pull out the treatment cutoff point at 0.05
res_F0_OWvAM <- results(dds_F0, name="treatment_OW_vs_AM", alpha=0.05)

res_F0_OWvAM <- res_F0_OWvAM[order(res_F0_OWvAM$padj),]
head(res_F0_OWvAM)

summary(res_F0_OWvAM)


res_F0_OWAvAM <- results(dds_F0, name="treatment_OWA_vs_AM", alpha=0.05)

res_F0_OWAvAM <- res_F0_OWAvAM[order(res_F0_OWAvAM$padj),]
head(res_F0_OWAvAM)

summary(res_F0_OWAvAM)


res_F0_OAvAM <- results(dds_F0, name="treatment_OA_vs_AM", alpha=0.05)

res_F0_OAvAM <- res_F0_OAvAM[order(res_F0_OAvAM$padj),]
head(res_F0_OAvAM)

summary(res_F0_OAvAM)


################## Save all the results as csv to go into GOMWU

library(tidyr)

# Make the rownames a separate column called transcriptID and make it all a dataframe
res_F0_OWvAM_df <- data.frame(transcriptID = rownames(res_F0_OWvAM), res_F0_OWvAM)

# Split the "transcriptID" column by double colons and create new columns of the parts----- getting rid of the extra part from the transcrptor program added
res_F0_OWvAM_df <- separate(res_F0_OWvAM_df, transcriptID, into = c("part1", "part2", "part3", "rest"), sep = "::", remove = FALSE) 

# Create a new column by concatenating "part1" and "part2" with double colons in between
res_F0_OWvAM_df$transcriptID_trim <- paste(res_F0_OWvAM_df$part1, res_F0_OWvAM_df$part2, sep = "::")

# Optional: Remove the "part1" and "part2" columns from the dataframe
res_F0_OWvAM_df <- res_F0_OWvAM_df[, !(names(res_F0_OWvAM_df) %in% c("part1", "part2", "part3", "rest"))]

write.table(res_F0_OWvAM_df, file = "res_F0_OWvAM.txt", sep = "\t", row.names = F)   # saves the full original for the records

# Select the two columns we want to save for the GOMWU analysis
selected_columns_OW <- res_F0_OWvAM_df[c("transcriptID_trim", "log2FoldChange")]

# Save the selected columns as a CSV file,,, quote false is important otherwise it won't read it
write.csv(selected_columns_OW, file = "res_F0_OWvAM_LFC.csv", quote = FALSE, row.names = F) # saves the selected columns for GOMWU


############ Now for OWA

# Make the rownames a separate column called transcriptID and make it all a dataframe
res_F0_OWAvAM_df <- data.frame(transcriptID = rownames(res_F0_OWAvAM), res_F0_OWAvAM)

# Split the "transcriptID" column by double colons and create new columns of the parts
res_F0_OWAvAM_df <- separate(res_F0_OWAvAM_df, transcriptID, into = c("part1", "part2", "part3", "rest"), sep = "::", remove = FALSE) 

# Create a new column by concatenating "part1" and "part2" with double colons in between
res_F0_OWAvAM_df$transcriptID_trim <- paste(res_F0_OWAvAM_df$part1, res_F0_OWAvAM_df$part2, sep = "::")

# Optional: Remove the "part1" and "part2" columns from the dataframe
res_F0_OWAvAM_df <- res_F0_OWAvAM_df[, !(names(res_F0_OWAvAM_df) %in% c("part1", "part2", "part3", "rest"))]
write.table(res_F0_OWAvAM_df, file = "res_F0_OWAvAM.txt", sep = "\t", row.names = F)   # saves the full original for the records

# Select the two columns we want to save for the GOMWU analysis
selected_columns_OWA <- res_F0_OWAvAM_df[c("transcriptID_trim", "log2FoldChange")]

# Save the selected columns as a CSV file
write.csv(selected_columns_OWA, file = "res_F0_OWAvAM_LFC.csv", quote = FALSE, row.names = F) # saves the selected columns for GOMWU



############ Now for OA

# Make the rownames a separate column called transcriptID and make it all a dataframe
res_F0_OAvAM_df <- data.frame(transcriptID = rownames(res_F0_OAvAM), res_F0_OAvAM)

# Split the "transcriptID" column by double colons and create new columns of the parts
res_F0_OAvAM_df <- separate(res_F0_OAvAM_df, transcriptID, into = c("part1", "part2", "part3", "rest"), sep = "::", remove = FALSE) 

# Create a new column by concatenating "part1" and "part2" with double colons in between
res_F0_OAvAM_df$transcriptID_trim <- paste(res_F0_OAvAM_df$part1, res_F0_OAvAM_df$part2, sep = "::")

# Optional: Remove the "part1" and "part2" columns from the dataframe
res_F0_OAvAM_df <- res_F0_OAvAM_df[, !(names(res_F0_OAvAM_df) %in% c("part1", "part2", "part3", "rest"))]
write.table(res_F0_OAvAM_df, file = "res_F0_OAvAM.txt", sep = "\t", row.names = F)   # saves the full original for the records

# Select the two columns we want to save for the GOMWU analysis
selected_columns_OA <- res_F0_OAvAM_df[c("transcriptID_trim", "log2FoldChange")]

# Save the selected columns as a CSV file
write.csv(selected_columns_OA, file = "res_F0_OAvAM_LFC.csv", quote = FALSE, row.names = F) # saves the selected columns for GOMWU

