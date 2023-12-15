library(raster)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(data.table)
library(gear)

library(RcppCNPy) # for reading python numpy (.npy) files

setwd("/Users/maryamnouri-aiin/Desktop/Ecological-Genomics/Pop genomics/results/")

list.files()

### read in selection statistics (these are chi^2 distributed)

s<-npyLoad("allRS_poly.selection.npy")

# convert test statistic to p-value
pval <- as.data.frame(1-pchisq(s,1))
names(pval) = c("p_PC1","p_PC2")

## read positions
p <- read.table("allRS_poly_mafs.sites",sep="\t",header=T, stringsAsFactors=T)
dim(p)

p_filtered = p[which(p$kept_sites==1),]
dim(p_filtered)

# get all the outliers with p-values below some cutoff
cutoff=1e-3   

outliers_PC1 <- p_filtered[which(pval$p_PC1<cutoff),c("chromo","position")]

# how many outlier loci < the cutoff?
dim(outliers_PC1)[1]

# write them out to a file
write.table(outliers_PC1,
            "allRS_poly_outliers_PC1.txt", 
            sep=":",
            quote=F,
            row.names=F,
            col.names=F)

COV <- as.matrix(read.table("allRS_poly.cov"))

PCA <- eigen(COV)

data=as.data.frame(PCA$vectors)
data=data[,c(1:2)] # the second number here is the number of PC axes you want to keep

write.table(data,
            "allRS_poly_genPC1_2.txt",
            sep="\t",
            quote=F,
            row.names=F,
            col.names=F)
# You made need to use "install.packages" if you don't have some of the below libraries already

library(raster)
library(FactoMineR)
library(factoextra)
library(corrplot)

setwd("/Users/maryamnouri-aiin/Desktop/Ecological-Genomics/Pop genomics/results/")

bio <- getData("worldclim",var="bio",res=10)
###########
bio <- getData("worldclim",var="bio",res=10)
bio10 <- bio[["bio10"]]
bio12 <- bio[["bio12"]]

coords <- read.csv("https://www.uvm.edu/~kellrlab/forClass/colebrookSampleMetaData.csv", header=T)

#The chunk below refers to your bamlist file that you transferred during last week's PCA/admixture analysis.  It should be the same one you want to use here -- if your sample list for analysis changes in the future, you'll need a different bamlist!
  
names <- read.table("allRS_bam.list")
names <- unlist(strsplit(basename(as.character(names[,1])), split = ".sorted.rmdup.bam"))
split = strsplit(names, "_")
pops <- data.frame(names[1:95], do.call(rbind, split[1:95]))
names(pops) = c("Ind", "Pop", "Row", "Col")
head(pops)
angsd_coords <- merge(pops, coords, by.x="Ind", by.y="Tree")

points <- SpatialPoints(angsd_coords[c("Longitude","Latitude")])
head(angsd_coords)
head(points)

Temperature <- extract(bio10,points)
Precipitation <- extract(bio12, points)

angsd_coords_clim <- cbind.data.frame(angsd_coords,Precipitation)
str(angsd_coords_clim)

# Make the climate PCA: we need to chose the best variable 

clim_PCA = PCA(angsd_coords_clim[,15:12], graph=T)

# Get a screeplot of cliamte PCA eigenvalues

fviz_eig(clim_PCA)
tiff('temp and prec', units="in", width=12, height=8, res=300)
# What is the climate PCA space our red spruce pops occupy?
fviz_pca_biplot(clim_PCA, 
                geom.ind="point",
                col.ind = angsd_coords_clim$Latitude, 
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                title="Climate PCA (Bioclim)",
                legend.title="Latitude")

dev.off()
# Which variables show the strongest correlation on the first 2 climate PC axes?

dimdesc(clim_PCA)[1:2]

# Replace "XX" with your bio variable most significant on climate PC1:

write.table(scale(angsd_coords_clim["bio10"]),
            "allRS_bio10.txt",
            sep="\t",
            quote=F,
            row.names = F,
            col.names=F)

# Replace "YY" with your bio variable most significant on climate PC2:  

write.table(scale(angsd_coords_clim["bio11"]),
            "allRS_bio12.txt",
            sep="\t",
            quote=F,
            row.names = F,
            col.names=F)

