## Set your working directory
setwd("/Users/maryamnouri-aiin/Desktop/Ecological-Genomics/Pop genomics/results/Transcriptome/")

## Import the libraries that we're likely to need in this session

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn) 
############A. Import Counts Matrix and Sample ID tables into R and DESeq2

# Import the counts matrix
countsTable <- read.table("salmon.xisoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample discription table
conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)

#####B. Explore the counts data with simple calculations and plots

# Let's see how many reads we have from each sample, total number of reads per sample
colSums(countsTableRound)
mean(colSums(countsTableRound)) 

barplot(colSums(countsTableRound), names.arg=colnames(countsTableRound),cex.names=0.5, las=3,ylim=c(0,20000000))
abline(h=mean(colSums(countsTableRound)), col="blue", lwd=2)

# the average number of counts per gene
rowSums(countsTableRound)
mean(rowSums(countsTableRound)) # [1] #8217.81 - hudsonica genes, 2269 - hudsonica isoform, 8218 - hudsonica filtered
median(rowSums(countsTableRound)) # [1] 377 - hudsonica, 109 - hudsonica isoforms, 377 - hudsonica filtered

apply(countsTableRound,2,mean) # 2 in the apply function does the action across columns
apply(countsTableRound,1,mean) # 1 in the apply function does the action across rows
hist(apply(countsTableRound,1,mean),xlim=c(0,10000), ylim=c(0,40000),breaks=10000)

############## Satrt working with DESeq ################

#### Create a DESeq object and define the experimental design here with the tilda

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ generation + treatment)

dim(dds)

# Filter out genes with too few reads - remove all genes with counts < 15 in more than 75% of samples, so ~28)
## suggested by WGCNA on RNAseq FAQ

dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 
# [1] 25260 that have at least 15 reads (a.k.a reads) in 75% of the samples

# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)
# produce a data matrix by resultsNames(dds)
#resultsNames(dds)
#[1] "Intercept"            "generation_F11_vs_F0" "generation_F2_vs_F0" 
#[4] "generation_F4_vs_F0"  "treatment_OA_vs_AM"   "treatment_OW_vs_AM"  
#[7] "treatment_OWA_vs_AM"

#####################  Check the quality of the data by sample clustering and visualization
########### D. Check the quality of the data by sample clustering and visualization
###############################################################

# Check the quality of the data by sample clustering and visualization
# The goal of transformation "is to remove the dependence of the variance on the mean, particularly the high variance of the logarithm of count data when the mean is low."


dds$treatment <- relevel(dds$treatment, ref = "AM")
dds$generation <- relevel(dds$generation, ref = "F0")


dds$group <- factor(paste0(dds$treatment, dds$generation))
design(dds) <- ~ group
dds <- DESeq(dds)
resultsNames(dds)
# [1] "Intercept"            "group_AMF11_vs_AMF0"  "group_AMF2_vs_AMF0"   "group_AMF4_vs_AMF0"  
# [5] "group_OAF0_vs_AMF0"   "group_OAF2_vs_AMF0"   "group_OAF4_vs_AMF0"   "group_OWAF0_vs_AMF0" 
# [9] "group_OWAF11_vs_AMF0" "group_OWAF4_vs_AMF0"  "group_OWF0_vs_AMF0"   "group_OWF2_vs_AMF0"  
# [13] "group_OWF4_vs_AMF0" 

res_F0_AMvOW <- results(dds, contrast=c("group","OWF0","AMF0"), alpha = 0.05)

res_F0_AMvOW <- res_F0_AMvOW[order(res_F0_AMvOW$padj),]
head(res_F0_AMvOW) 
# log2 fold change (MLE): group OWF0 vs AMF0 
# Wald test p-value: group OWF0 vs AMF0 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   TRINITY_DN83225_c0_g1_i29   209.988      -24.26072  2.415427 -10.04407 9.75617e-24 4.15320e-19
# TRINITY_DN2142_c0_g1_i6     355.201      -11.41509  1.563258  -7.30211 2.83285e-13 6.02971e-09
# TRINITY_DN37167_c0_g1_i7    837.120        2.27357  0.330139   6.88672 5.70947e-12 8.10174e-08
# TRINITY_DN48447_c0_g1_i4    461.021      -10.89509  1.628017  -6.69224 2.19775e-11 2.33896e-07
# TRINITY_DN10274_c0_g1_i6    135.360        2.70207  0.410322   6.58525 4.54134e-11 3.17091e-07
# TRINITY_DN4949_c0_g1_i1    1213.253        3.81937  0.581805   6.56469 5.21409e-11 3.17091e-07

summary(res_F0_AMvOW)
library("pheatmap")
library("vsn")

# this gives log2(n + 1) # Normalization 
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))

# Variance stabilizing transformation # Reduces the low reads transcriptomes
vsd <- vst(dds, blind=FALSE)
meanSdPlot(assay(vsd))


sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$treatment, vsd$generation, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# AM_F11_Rep3 and OA_F2_Rep2 could be outliers
#OW ocean warming, AM ambient, OA ocean acidification, OWA Ocean Warming Acidification

############## Visualize the global gene expression patterns using PCA
###############################################################

# PCA to visualize global gene expression patterns

# first transform the data for plotting using variance stabilization
vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("treatment","generation"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData,"percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=treatment, shape=generation)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
 

#### cluster more with generation rather than treatment
### stong gene experssion response, second gene establish, and again F11

###Make a more advanced PCA plot and plot by generation in four panels

############################################################### 

# Let's plot the PCA by generation in four panels

data <- plotPCA(vsd, intgroup=c("treatment","generation"), returnData=TRUE)
percentVar <- round(100 * attr(data,"percentVar"))

###########  

dataF0 <- subset(data, generation == 'F0')

F0 <- ggplot(dataF0, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-10, 25) + xlim(-40, 10)+ # zoom for F0 with new assembly
  #ylim(-40, 25) + xlim(-50, 50)+ # new assembly limits
  #ylim(-40, 20) + xlim(-50, 30)+
  scale_shape_manual(values=c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  ##theme(legend.position = c(0.83,0.85), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) +
  #guides(shape = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  #guides(fill = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  #guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank())

F0


#png("PCA_F0.png", res=300, height=5, width=5, units="in")

#ggarrange(F0, nrow = 1, ncol=1)

#dev.off()

################# F2

dataF2 <- subset(data, generation == 'F2')

F2 <- ggplot(dataF2, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-40, 25) + xlim(-50, 55)+
  #ylim(-40, 20) + xlim(-50, 30)+
  scale_shape_manual(values=c(21,22,23), labels = c("Ambient", "Acidification","Warming"))+
  # scale_color_manual(values = c('#6699CC',"#F2AD00","#00A08A", "#CC3333")) + 
  #scale_color_manual(values=c('black')) +
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A"), labels = c("Ambient", "Acidification","Warming"))+
  theme(legend.position = c(0.83,0.85), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) +
  #scale_size(guide="none") +
  guides(shape = guide_legend(override.aes = list(shape = c( 21,22, 23))))+
  guides(fill = guide_legend(override.aes = list(shape = c( 21,22, 23))))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank())
F2



###### one ambient sample is missing  in the data set

# png("PCA_F2.png", res=300, height=5, width=5, units="in")
# 
# ggarrange(F2, nrow = 1, ncol=1)
# 
# dev.off()

# Yes - F2 is missing one ambient replicate

################################ F4

dataF4 <- subset(data, generation == 'F4')

F4 <- ggplot(dataF4, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-40, 25) + xlim(-50, 55)+ # limits with filtered assembly
  #ylim(-20, 10) + xlim(-40, 25)+  # zoom with filtered assembly
  #ylim(-40, 20) + xlim(-50, 30)+
  scale_shape_manual(values=c(21,22,23,24), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  # scale_color_manual(values = c('#6699CC',"#F2AD00","#00A08A", "#CC3333")) + 
  #scale_color_manual(values=c('black')) +
  scale_fill_manual(values=c('#6699CC',"#F2AD00","#00A08A", "#CC3333"), labels = c("Ambient", "Acidification","Warming", "OWA"))+
  #theme(legend.position = c(0.83,0.85), legend.background = element_blank(), legend.box.background = element_rect(colour = "black")) +
  #scale_size(guide="none") +
  guides(shape = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  guides(fill = guide_legend(override.aes = list(shape = c( 21,22, 23, 24))))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank())
F4


# png("PCA_F4.png", res=300, height=5, width=5, units="in")
# 
# ggarrange(F4, nrow = 1, ncol=1)
# 
# dev.off()


################# F11

dataF11 <- subset(data, generation == 'F11')

F11 <- ggplot(dataF11, aes(PC1, PC2)) +
  geom_point(size=10, stroke = 1.5, aes(fill=treatment, shape=treatment)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylim(-50, 25) + xlim(-50, 55)+
  #ylim(-100, 20) + xlim(-100, 30)+
  scale_shape_manual(values=c(21,24), labels = c("Ambient", "OWA"))+
  scale_fill_manual(values=c('#6699CC', "#CC3333"), labels = c("Ambient", "OWA"))+
  guides(shape = guide_legend(override.aes = list(shape = c( 21, 24))))+
  guides(fill = guide_legend(override.aes = list(shape = c( 21, 24))))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 4))+
  theme(text = element_text(size = 20)) +
  theme(legend.title = element_blank())
F11

############## fix all axes size before the final plot
# png("PCA_F11.png", res=300, height=5, width=5, units="in")
# 
# ggarrange(F11, nrow = 1, ncol=1)
# 
# dev.off()

ggarrange(F0, F2, F4, F11, nrow = 2, ncol=2)


##########E. Explore differential expression with full model, GE ~ generation + treatment
###########Now that we have a sense of our data, let’s explore models and contrasts for testing for differential expression

## Check on the DE results from the DESeq command way above ##

#resultsNames(dds)
#[1] "Intercept"            "generation_F11_vs_F0" "generation_F2_vs_F0" 
#[4] "generation_F4_vs_F0"  "treatment_OA_vs_AM"   "treatment_OW_vs_AM"  
#[7] "treatment_OWA_vs_AM"
resAM_OWA <- results(dds, name="treatment_OWA_vs_AM", alpha=0.05)

resAM_OWA <- resAM_OWA[order(resAM_OWA$padj),]
head(resAM_OWA)  

# log2 fold change (MLE): treatment OWA vs AM 
# Wald test p-value: treatment OWA vs AM 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue
# <numeric>      <numeric> <numeric> <numeric>   <numeric>
#   TRINITY_DN22748_c0_g1::TRINITY_DN22748_c0_g1_i4::g.47585::m.47585   209.899      -1.799165  0.282811  -6.36172 1.99509e-10
# TRINITY_DN3600_c0_g1::TRINITY_DN3600_c0_g1_i2::g.16079::m.16079      91.998      -2.256764  0.355515  -6.34788 2.18302e-10
# TRINITY_DN505_c0_g1::TRINITY_DN505_c0_g1_i6::g.3715::m.3715         570.049      -0.698475  0.112387  -6.21493 5.13460e-10
# TRINITY_DN14069_c0_g1::TRINITY_DN14069_c0_g1_i2::g.39418::m.39418   193.197       3.555491  0.625276   5.68627 1.29844e-08
# TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i14::g.733::m.733            699.038       2.606874  0.457072   5.70342 1.17424e-08
# TRINITY_DN8532_c1_g1::TRINITY_DN8532_c1_g1_i3::g.30152::m.30152     184.374       0.663567  0.119300   5.56217 2.66438e-08
summary(resAM_OWA)
# out of 25260 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 99, 0.39% up regulated
# LFC < 0 (down)     : 75, 0.3% Down regulated
# outliers [1]       : 4, 0.016%
# low counts [2]     : 0, 0%
# (mean count < 17)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# ocean warming vs ambient
resAM_OW <- results(dds, name="treatment_OW_vs_AM", alpha=0.05)

resAM_OW <- resAM_OW[order(resAM_OW$padj),]
head(resAM_OW)  

# log2 fold change (MLE): treatment OW vs AM 
# Wald test p-value: treatment OW vs AM 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i3::g.747::m.747             471.0065        1.75150  0.307726   5.69173 1.25756e-08 0.000286812
# TRINITY_DN17372_c0_g1::TRINITY_DN17372_c0_g1_i7::g.43175::m.43175    51.6168        1.58846  0.286984   5.53500 3.11221e-08 0.000354901
# TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i3::g.744::m.744            3876.7153        1.42685  0.262666   5.43219 5.56671e-08 0.000423200
# TRINITY_DN1275_c0_g1::TRINITY_DN1275_c0_g1_i9::g.6875::m.6875       168.9233       -1.64422  0.341342  -4.81694 1.45779e-06 0.008311970
# TRINITY_DN17372_c0_g1::TRINITY_DN17372_c0_g1_i11::g.43161::m.43161   87.3986        2.49827  0.525888   4.75057 2.02842e-06 0.009252455
# TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i14::g.733::m.733            699.0377        2.20090  0.471974   4.66318 3.11356e-06 0.011442688
summary(resAM_OW)

# out of 25260 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 9, 0.036%
# LFC < 0 (down)     : 8, 0.032%
# outliers [1]       : 4, 0.016%
# low counts [2]     : 2449, 9.7%
# (mean count < 31)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


### Plot Individual genes ### 

# Counts of specific top interaction gene! (important validatition that the normalization, model is working)

## The gen can be changed from the code above "TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i3::g.747::m.747"
d <-plotCounts(dds, gene="TRINITY_DN29_c1_g2::TRINITY_DN29_c1_g2_i3::g.747::m.747", intgroup = (c("treatment","generation")), returnData=TRUE)
d

p <-ggplot(d, aes(x=treatment, y=count, color=treatment, shape=generation)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p


################3F. Run another model that focuses on treatment effects for each generation

#################### MODEL NUMBER 2 - subset to focus on effect of treatment for each generation

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ treatment)

dim(dds)
# [1] 130580     38

# Filter 
dds <- dds[rowSums(counts(dds) >= 30) >= 28,]
nrow(dds) 

# Subset the DESeqDataSet to the specific level of the "generation" factor
dds_sub <- subset(dds, select = generation == 'F0')
dim(dds_sub)

# Perform DESeq2 analysis on the subset
dds_sub <- DESeq(dds_sub)

resultsNames(dds_sub)

res_F0_OWvAM <- results(dds_sub, name="treatment_OW_vs_AM", alpha=0.05)

res_F0_OWvAM <- res_F0_OWvAM[order(res_F0_OWvAM$padj),]
head(res_F0_OWvAM)

summary(res_F0_OWvAM)


### Plot Individual genes ### 

# Counts of specific top interaction gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds_sub, gene="TRINITY_DN30_c0_g2::TRINITY_DN30_c0_g2_i1::g.130::m.130", intgroup = (c("treatment","generation")), returnData=TRUE)
d

p <-ggplot(d, aes(x=treatment, y=count, color=treatment, shape=generation)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p

# higher expression in Ocean Warming
#Lower in Ambient and Ocean acidification
# strong clustering in ambient in all generation
#IN OWA F4 low experssion 


############ F. Run another model that focuses on treatment effects for each generation
#################### MODEL NUMBER 2 - subset to focus on effect of treatment for each generation
# remake dds by treatments
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ treatment)

dim(dds)
# [1]  67916    38

# Filter 
dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 
# [1] 25260    12

# Subset the DESeqDataSet to the specific level of the "generation" factor
dds_sub <- subset(dds, select = generation == 'F0')
dim(dds_sub)
# [1] 25260    12
# Perform DESeq2 analysis on the subset
dds_sub <- DESeq(dds_sub)

resultsNames(dds_sub)
# [1] "Intercept"           "treatment_OA_vs_AM"  "treatment_OW_vs_AM"  "treatment_OWA_vs_AM"


## look at OW and AM, they were far apart on PCA
res_F0_OWvAM <- results(dds_sub, name="treatment_OW_vs_AM", alpha=0.05)

res_F0_OWvAM <- res_F0_OWvAM[order(res_F0_OWvAM$padj),]
head(res_F0_OWvAM)

# log2 fold change (MLE): treatment OW vs AM 
# Wald test p-value: treatment OW vs AM 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE
# <numeric>      <numeric> <numeric>
#   TRINITY_DN30_c0_g2::TRINITY_DN30_c0_g2_i1::g.130::m.130           32956.030        3.74319  0.292239
# TRINITY_DN33_c0_g1::TRINITY_DN33_c0_g1_i18::g.235::m.235           4622.325        3.72565  0.296419
# TRINITY_DN9365_c0_g1::TRINITY_DN9365_c0_g1_i1::g.32056::m.32056    3814.681        3.62002  0.292971
# TRINITY_DN258_c0_g1::TRINITY_DN258_c0_g1_i18::g.1704::m.1704       4670.245        3.72840  0.317074
# TRINITY_DN258_c0_g1::TRINITY_DN258_c0_g1_i13::g.1700::m.1700       5937.496        3.30257  0.280896
# TRINITY_DN10209_c0_g2::TRINITY_DN10209_c0_g2_i3::g.33784::m.33784   279.313        3.00313  0.264004
# stat      pvalue        padj
# <numeric>   <numeric>   <numeric>
#   TRINITY_DN30_c0_g2::TRINITY_DN30_c0_g2_i1::g.130::m.130             12.8086 1.46698e-37 3.63048e-33
# TRINITY_DN33_c0_g1::TRINITY_DN33_c0_g1_i18::g.235::m.235            12.5689 3.13154e-36 3.87496e-32
# TRINITY_DN9365_c0_g1::TRINITY_DN9365_c0_g1_i1::g.32056::m.32056     12.3562 4.50893e-35 3.71957e-31
# TRINITY_DN258_c0_g1::TRINITY_DN258_c0_g1_i18::g.1704::m.1704        11.7588 6.36643e-32 3.20614e-28
# TRINITY_DN258_c0_g1::TRINITY_DN258_c0_g1_i13::g.1700::m.1700        11.7573 6.47758e-32 3.20614e-28
# TRINITY_DN10209_c0_g2::TRINITY_DN10209_c0_g2_i3::g.33784::m.33784   11.3753 5.55141e-30 2.28977e-26

summary(res_F0_OWvAM)

# out of 25260 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3468, 14%
# LFC < 0 (down)     : 2049, 8.1%
# outliers [1]       : 22, 0.087%
# low counts [2]     : 490, 1.9%
# (mean count < 20)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

### Plot Individual genes ### 

# Counts of specific top interaction gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds_sub, gene="TRINITY_DN30_c0_g2::TRINITY_DN30_c0_g2_i1::g.130::m.130", intgroup = (c("treatment","generation")), returnData=TRUE)
d

p <-ggplot(d, aes(x=treatment, y=count, color=treatment, shape=generation)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p
###### upregulated in response to warming, major diffrentially expressed gene

### try another gene
d <-plotCounts(dds_sub, gene="TRINITY_DN10209_c0_g2::TRINITY_DN10209_c0_g2_i3::g.33784::m.33784", intgroup = (c("treatment","generation")), returnData=TRUE)
d

p <-ggplot(d, aes(x=treatment, y=count, color=treatment, shape=generation)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p


######################### We can make an MA (or sideways volcano) plot now that we’re contrasting just two groups

# We can make an MA plot
plotMA(res_F0_OWvAM, ylim=c(-5,5))
#triangle when they are higher
#We can also make a heatmap of the DEGs
########################################

# Heatmap of top 20 genes sorted by pvalue

library(pheatmap)

# By environment
vsd <- vst(dds_sub, blind=FALSE)

###We can change the op number 
topgenes <- head(rownames(res_F0_OWvAM),100)

topgenes <- head(rownames(res_F0_OWvAM),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
## normalizing by substracting means to fit the dat in a heat map
df <- as.data.frame(colData(dds_sub)[,c("generation","treatment")])
pheatmap(mat, annotation_col=df)
pheatmap(mat, annotation_col=df, cluster_cols = F)
### dendograms are clustering gene together

#Let’s make a Venn (or Euler) Diagram
#Last but not least, since we have three contrasts, let’s make a Venn (or Euler) Diagram to see how similar or different the DGE is from Ambient for OA, OW, and OWA at F0!
  
#################################################################

#### PLOT OVERLAPPING DEGS IN VENN EULER DIAGRAM

#################################################################

# For OW vs AM
res_F0_OWvAM <- results(dds_sub, name="treatment_OW_vs_AM", alpha=0.05)
res_F0_OWvAM <- res_F0_OWvAM[order(res_F0_OWvAM$padj),]
head(res_F0_OWvAM)

summary(res_F0_OWvAM)
# out of 25260 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 3468, 14%
# LFC < 0 (down)     : 2049, 8.1%
# outliers [1]       : 22, 0.087%
# low counts [2]     : 490, 1.9%
# (mean count < 20)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
res_F0_OWvAM <- res_F0_OWvAM[!is.na(res_F0_OWvAM$padj),]
degs_F0_OWvAM <- row.names(res_F0_OWvAM[res_F0_OWvAM$padj < 0.05,])
#degs are diffrential expressed genes
# For OA vs AM
res_F0_OAvAM <- results(dds_sub, name="treatment_OA_vs_AM", alpha=0.05)
res_F0_OAvAM <- res_F0_OAvAM[order(res_F0_OAvAM$padj),]
head(res_F0_OAvAM)

summary(res_F0_OAvAM)
res_F0_OAvAM <- res_F0_OAvAM[!is.na(res_F0_OAvAM$padj),]
degs_F0_OAvAM <- row.names(res_F0_OAvAM[res_F0_OAvAM$padj < 0.05,])


########3 we could separately do the diagram for unregulated  or downregulated by filtering the results by log fold or statstics, here is a mix
# For OWA vs AM
res_F0_OWAvAM <- results(dds_sub, name="treatment_OWA_vs_AM", alpha=0.05)
res_F0_OWAvAM <- res_F0_OWAvAM[order(res_F0_OWAvAM$padj),]
head(res_F0_OWAvAM)

summary(res_F0_OWAvAM)
res_F0_OWAvAM <- res_F0_OWAvAM[!is.na(res_F0_OWAvAM$padj),]
degs_F0_OWAvAM <- row.names(res_F0_OWAvAM[res_F0_OWAvAM$padj < 0.05,])

library(eulerr)

# Total
length(degs_F0_OAvAM)  # 602
length(degs_F0_OWvAM)  # 5517 
length(degs_F0_OWAvAM)  # 3918

# Intersections
length(intersect(degs_F0_OAvAM,degs_F0_OWvAM))  # 387, 444
length(intersect(degs_F0_OAvAM,degs_F0_OWAvAM))  # 340, 380
length(intersect(degs_F0_OWAvAM,degs_F0_OWvAM))  # 2585, 2743

intWA <- intersect(degs_F0_OAvAM,degs_F0_OWvAM)
length(intersect(degs_F0_OWAvAM,intWA)) # 308, 338

# Number unique

602-444-380+338 # 101, 116 OA
5517-444-2743+338 # 2177, 2668 OW 
3918-380-2743+338 # 1125, 1133 OWA

444-338 # 79, 106 OA & OW
380-338 # 32, 42 OA & OWA
2743-338 # 2277, 2405 OWA & OW


# Note that the names are important and have to be specific to line up the diagram
fit1 <- euler(c("OA" = 116, "OW" = 2668, "OWA" = 1133, "OA&OW" = 106, "OA&OWA" = 42, "OW&OWA" = 2405, "OA&OW&OWA" = 338))


plot(fit1,  lty = 1:3, quantities = TRUE)
# lty changes the lines

plot(fit1, quantities = TRUE, fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))


#cross check
2668+2405+338+106 # 4841, 5517 total OW should be the same number of total DFE s in this code # length(degs_F0_OWvAM)  # 5517 
1133+2405+338+42 # 3742, 3918 total OWA
116+42+106+338    # 520, 602 total OA
