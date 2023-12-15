###################################
#  Selection scans for red spruce #
###################################

library(RcppCNPy) # for reading python numpy (.npy) files

setwd("/Users/maryamnouri-aiin/Desktop/Ecological-Genomics/Pop genomics/results/")

list.files()

### read in selection statistics (these are chi^2 distributed)

s<-npyLoad("allRS_poly.selection.npy")
#number of loci
#the cols are the e number, we chose to do 2
head(s)
dim(s)

# convert test statistic to p-value
pval <- as.data.frame(1-pchisq(s,1))
# head for pvalue
head(pval)
dim(pval)

# different from class note, we added one more col/PC
names(pval) = c("p_PC1", "p_PC2")

## read positions
p <- read.table("allRS_poly_mafs.sites",sep="\t",header=T, stringsAsFactors=T)
# number of loci
dim(p)
head(p)

#filter data frame by selecting rows that has 1
p_filtered = p[which(p$kept_sites==1),]
dim(p_filtered)
head(p_filtered)

# How many sites got filtered out when testing for selection? Why?

## make manhattan plot
plot(-log10(pval$p_PC1),
     col=p_filtered$chromo,
     xlab="Position",
     ylab="-log10(p-value)",
     main="Selection outliers: pcANGSD e=1 (K2)")

# We can zoom in if there's something interesting near a position...

plot(-log10(pval$p_PC1[2e05:2.01e05]),
     col=p_filtered$chromo, 
     xlab="Position", 
     ylab="-log10(p-value)", 
     main="Selection outliers: pcANGSD e=1 (K2)")

# get the contig with the lowest p-value for selection
sel_contig <- p_filtered[which(pval==min(pval$p_PC1)),c("chromo","position")]
sel_contig

# get all the outliers with p-values below some cutoff
cutoff=1e-4   # equals a 1 in 5,000 probability
outlier_contig <- p_filtered[which(pval$p_PC1<cutoff),c("chromo","position")]
outlier_contig

outlier_contig[2000-2500, 1:2]

dim(outlier_contig[which(outlier_contig$position>0),])

outlier_contig <-(outlier_contig[which(outlier_contig$position>0),])
dim(outlier_contig)  

# how many unique contigs harbor outlier loci?
length(unique(outlier_contig$chromo))

write.table(unique(outlier_contig$chromo),
            "allRS_poly_PC1_outlier_contig.txt", 
            sep="\t",
            quote=F,
            row.names=F,
            col.names=F)
# how many outlier loci < the cutoff?
dim(outlier_contigs)[1]




