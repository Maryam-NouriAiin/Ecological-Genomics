library(ggplot2) # plotting
library(ggpubr) # plotting

setwd("/Users/maryamnouri-aiin/Desktop/Ecological-Genomics/Pop genomics/homeworks/") # set the path to where you saved the pcANGSD results on your laptop

## First, let's work on the genetic PCA:

COV <- as.matrix(read.table("allRS_poly.cov")) # read in the genetic covariance matrix

PCA <- eigen(COV) # extract the principal components from the COV matrix
pca <- prcomp(COV)
eigenvalues <- pca$sdev^2
plot(eigenvalues, type= "b", xlab="No. of PCs", ylab="Variance Explained")
elbow <- 0
for(i in 2:length(eigenvalues)) {gradient <- eigenvalues[i-1]-
  eigenvalues[i]
if (gradient<0.1*eigenvalues[i-1]){
  elbow<-i-1
  break
}
}
cat("Number of retained PCs:", elbow,"\n")

## How much variance is explained by the first few PCs?

var <- round(PCA$values/sum(PCA$values),3)

var[1:3]

# A "screeplot" of the eigenvalues of the PCA:

barplot(var, 
        xlab="Eigenvalues of the PCA", 
        ylab="Proportion of variance explained")
## notes on screeplot: xlab and ylab: SNP1 and 2
### in the plot, 5% of variations (and the most)PCA1 explained, next bar is PCA2, third is PCA3

## Bring in the bam.list file and extract the sample info:

names <- read.table("allRS_bam.list")
names <- unlist(strsplit(basename(as.character(names[,1])), split = ".sorted.rmdup.bam"))
split = strsplit(names, "_")
pops <- data.frame(names[1:95], do.call(rbind, split[1:95]))
names(pops) = c("Ind", "Pop", "Row", "Col")

## A quick and humble PCA plot:

plot(PCA$vectors[,1:2],
     col=as.factor(pops[,2]),
     xlab="PC1",ylab="PC2", 
     main="Genetic PCA")
### same as beautiful plot next one

## A more beautiful PCA plot using ggplot :)

data=as.data.frame(PCA$vectors)
data=data[,c(1:3)]
data= cbind(data, pops)

cols=c("#377eB8","#EE9B00","#0A9396","#94D2BD","#FFCB69","#005f73","#E26D5C","#AE2012", "#6d597a", "#7EA16B","#d4e09b", "gray70")

ggscatter(data, x = "V1", y = "V2",
          color = "Pop",
          mean.point = TRUE,
          star.plot = TRUE) +
  theme_bw(base_size = 13, base_family = "Times") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text=element_text(size=rel(.7)), 
        axis.text = element_text(size=13), 
        legend.position = "bottom") +
  labs(x = paste0("PC1: (",var[1]*100,"%)"), y = paste0("PC2: (",var[2]*100,"%)")) +
  scale_color_manual(values=c(cols), name="Source population") +
  guides(colour = guide_legend(nrow = 2))

#### Results: on PCA 2 2021 is seprate from the other


## Next, we can look at the admixture clustering:

# import the ancestry scores (these are the .Q files)

# plot With k=3

# Set up the layout for the multi-panel figure

tiff('K1', units="in", width=8, height=10, res=300)
par(mfrow = c(3, 1))
q2 <- read.table("allRS_poly.admix.2.Q", sep=" ", header=F)

K2=dim(q)[2] #Find the level of K modeled

## order according to population code
ord<-order(pops[,2])

# make the plot:
barplot(t(q2)[,ord],
        col=cols[1:K2],
        space=0,border=NA,
        xlab="",ylab="Admixture proportions",
        main=paste0("K=",K))
text(tapply(1:nrow(pops),pops[ord,2],mean),-0.05,unique(pops[ord,2]),xpd=T)
abline(v=cumsum(sapply(unique(pops[ord,2]),function(x){sum(pops[ord,2]==x)})),col=1,lwd=1.2)



q3 <- read.table("allRS_poly.admix.3.Q", sep=" ", header=F)

K3=dim(q3)[2] #Find the level of K modeled

## order according to population code
ord<-order(pops[,2])
# make the plot:
K02 <- barplot(t(q3)[,ord],
        col=cols[1:K3],
        space=0,border=NA,
        xlab="",ylab="Admixture proportions",
        main=paste0("K=",K3))
text(tapply(1:nrow(pops),pops[ord,2],mean),-0.05,unique(pops[ord,2]),xpd=T)
abline(v=cumsum(sapply(unique(pops[ord,2]),function(x){sum(pops[ord,2]==x)})),col=1,lwd=1.2)


## 2025 for example is a mix, could be because they are closer to black spruce in Canada
## this is plot has the exact data in PCA plot
# Each thin bar is an individual
# Admixture proportions: gene flow, back crossing 
# o x b
# o x F1 (prop 0/0.5)
# o x bc (0.75)
# BC2 (0.875)
# Comparing Fst with this plot: Fst reduces as admixture increases e.g. in 2025 and 2103....
# The lower Fst is the less genetic variation, coulde be the reason of isolation, low genetic diversity
# less diversity in southern states, hybridization is happening in northern states (hybridization is quicker in introducing variation than mutation)


# plot With k=3
# need to first change k or e number in bash code in vim
q4 <- read.table("allRS_poly.admix.4.Q", sep=" ", header=F)

K4=dim(q4)[2] #Find the level of K modeled

## order according to population code
ord<-order(pops[,2])
# make the plot:
K03 <- barplot(t(q4)[,ord],
        col=cols[1:K4],
        space=0,border=NA,
        xlab="Populations",ylab="Admixture proportions",
        main=paste0("K=",K4))
text(tapply(1:nrow(pops),pops[ord,2],mean),-0.05,unique(pops[ord,2]),xpd=T)
abline(v=cumsum(sapply(unique(pops[ord,2]),function(x){sum(pops[ord,2]==x)})),col=1,lwd=1.2)
dev.off()

# plot With k=4
# need to first change k or e number in bash code in vim

