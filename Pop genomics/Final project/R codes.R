setwd("/Users/maryamnouri-aiin/Desktop/CompostMetagenomics/16S/")

require(readxl)
library(openxlsx)
library(tidyverse)
library(dplyr)
library(nloptr)
library(ggplot2)
library(gridExtra)
library(wesanderson)
library(RColorBrewer)
library("viridis")  
library(readr)
library(scales)
library(igraph)
library(ggridges)
library(agricolae)
library(patchwork)
library(qiime2R)
library(phyloseq)
library(BiocManager)
library(vegan)
library("ape")
library("DESeq2")
library(adaptiveGPCA)
library(phyloseqGraphTest)
library(igraph)


####################################Rarefaction and normalisation####################################
#Rawdata
data <- read.xlsx("16S_1.xlsx", sheet = 1, colNames = TRUE, rowNames = TRUE)
# Exclude the specified column "taxonomy
column_to_exclude <- c("taxonomy")
data <- data[, !colnames(data) %in% column_to_exclude]
data <- data[rowSums(data)>180,]
data <- data[,colSums(data)>10000]
OTU <- specnumber(t(data))
# Determine the rarefaction level
raremax <- min(colSums(data))
# Perform rarefaction
ACrare <- rarefy(t(data), raremax)
#Plot rarified data
plot(OTU, ACrare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")

# Customize the rarefaction curve
AC_Curve <- rarecurve(t(data), step = 20, sample = raremax, col = "blue", cex = 0.6)
# Perform relative rarefaction
AC_Data_rare <- rrarefy(t(data), raremax)
#.csv file to use for the rest of analysis
write.csv(t(AC_Data_rare), file = "AC_16S_Rfy.csv")

####################################https://vaulot.github.io/tutorials/Phyloseq_tutorial.html#step-by-step####################################
####################################Phyloseq####################################
#Three tables are needed: OTU, Taxonomy, Samples
otu_mat<- as.data.frame(read_excel("16Sdata.xlsx", sheet = "OTU matrix"))
tax_mat<- read_excel("16Sdata.xlsx", sheet = "Taxonomy table")
samples_df <- as.data.frame(read_excel("16Sdata.xlsx", sheet = "Samples"))

otu_mat <- otu_mat %>%
  tibble::column_to_rownames("ASV #")
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("ASV #")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("sample") 
#Transform into matrixes otu and tax tables (sample table can be left as data frame)
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
#Find the unique phylum and taxa
otumat <- data.frame(otu_mat)
taxmat <- data.frame(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

physeq0 <- phyloseq(OTU, TAX)
physeq <- phyloseq(OTU, TAX, samples)
physeq
unique(taxmat$Phylum)

#Rarefy_even_depth: Resample an OTU table such that all samples have the same
physeq.prune = prune_taxa(taxa_sums(physeq) > 180, physeq)
physeq.prune

#Unique phyla
unique(taxmat$Phylum)

#Extract the distance matrix from the phyloseq object
my_distance_matrix <- phyloseq::distance(physeq.prune, method = "bray")
head(my_distance_matrix)

#Assuming your sample data is called "sample_data"
adonis_result <- adonis2(my_distance_matrix ~ samples$Treatment +samples$PlantType + samples$SampleType, data = as(samples, "data.frame"))
adonis_result

#Keep only specific phylum
physeqtax <- subset_taxa(physeq.prune, Phylum %in% c("Acidobacteriota", "Actinobacteriota", "Bacteroidota", 
                                                    "Chloroflexi", "Gemmatimonadota", "Planctomycetota", 
                                                    "Verrucomicrobiota", "Proteobacteria", "Nitrospirota", 
                                                    "Firmicutes", "Cynobacteria"))
#Exclude classes
physeqtax <- subset_taxa(physeqtax, !(Class %in% c("BD7-11", "FFCH5909")))
# Exclude NA values for Genus and Species
physeqtax <- subset_taxa(physeqtax, !is.na(Genus) & !is.na(Species))
# Exclude Genus 
physeqtax <- subset_taxa(physeqtax, !(Genus %in% c("JG30-KF-CM45", "CPla-3 termite group", "67-14", "SRB2", "BIrii41", "env.OPS 17", 
                                                   "WWH38", "AKYG1722", "Amb-16S-1323", "A4b", "TRA3-20", "Unknown Family", "NA")))
# Exclude Species 
physeqtax <- subset_taxa(physeqtax, !(Species %in% c("Ellin516", "DEV008",  "SH-PL14", "JGI 0001001-H03", "BIyi10", "JCM 18997", 
                                                     "Ellin517", "NA")))
physeqtax


###Create violin plot of alpha diversity
alphadiversity_treatment_sampletype <- physeqtax %>%                                                            
  plot_richness(
    x = "Treatment",
    measures = c("Shannon")) +                      
  geom_violin(aes(fill = Treatment), show.legend = TRUE) +         
  geom_boxplot(width=0.1) +  
  facet_grid(PlantType~SampleType, ) +                               
  theme_light()+                                                   
  xlab(NULL)+ 
  
  theme(axis.text.y.left = element_text(size = 15),                 
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),           
        axis.title.y = element_text(size = 15))+                     
  theme(strip.text = element_text(face = "bold", size = 20))+        
  scale_fill_manual(values=c("lightgoldenrod","darkorange", "darkred"))+   
  ggtitle("") +                                        #add title
  theme(plot.title=element_text(size = 20, face = "bold", hjust = 0.5))
alphadiversity_treatment_sampletype

#Making a PCoA plot

all_pcoa <- ordinate(
  physeq = physeqtax, 
  method = "PCoA", 
  distance = "jaccard"
)
#Plot
pca_sampletype <- plot_ordination(
  physeq = physeqtax,
  ordination = all_pcoa) +
  geom_point(aes(fill = Treatment, shape = SampleType, color = PlantType), size = 15, stroke=4, alpha =1) +  
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_manual(values = c("lightgoldenrod", "darkorange", "#CC3333")) +
  scale_color_manual(values = c("#673770","#5F7FC7","darkseagreen")) +  
  scale_size_manual(values = c(3, 4, 5)) +  
  theme_classic() +
  theme(
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    axis.text.y.left = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme(legend.key.size = unit(2, "lines"))
pca_sampletype

pca_planttype <- plot_ordination(
  physeq = physeqtax,
  ordination = all_pcoa) +
  geom_point(aes(fill = Treatment, shape = PlantType, color=SampleType), size = 15, stroke=4, alpha =1) +  
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_manual(values = c("lightgoldenrod", "darkorange", "#CC3333")) +
  scale_color_manual(values = c("#673770","#5F7FC7","darkseagreen")) +  
  scale_size_manual(values = c(3, 4, 5)) +  
  theme_classic() +
  theme(
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    axis.text.y.left = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme(legend.key.size = unit(2, "lines"))
pca_planttype

pca_treatment <- plot_ordination(
  physeq = physeqtax,
  ordination = all_pcoa
) +
  geom_point(aes(color = Treatment), size = 15) + 
  scale_color_manual(values = c("lightgoldenrod", "darkorange", "#CC3333")) +
  theme_classic() +
  theme(
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    axis.text.y.left = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  ) +
  theme(legend.key.size = unit(2, "lines"))

pca_treatment


pca_SampleType <- plot_ordination(
  physeq = physeqtax,
  ordination = all_pcoa
) +
  geom_point(aes(color = SampleType), size = 15) + 
  scale_color_manual(values = c("lightgoldenrod", "darkorange", "#CC3333")) +
  theme_classic() +
  theme(
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    axis.text.y.left = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  ) +
  theme(legend.key.size = unit(2, "lines"))

pca_SampleType

pca_PlantType <- plot_ordination(
  physeq = physeqtax,
  ordination = all_pcoa
) +
  geom_point(aes(color = PlantType), size = 15) + 
  scale_color_manual(values = c("lightgoldenrod", "darkorange", "#CC3333")) +
  theme_classic() +
  theme(
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    axis.text.y.left = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  ) +
  theme(legend.key.size = unit(2, "lines"))

pca_PlantType

#Filter out eukaryotes and mitochondria
justbacteria <- physeqtax %>%
  subset_taxa(
    Kingdom == "Bacteria" 
  )
justbacteria

#Convert phyloseq object into a data frame with relative abundance counts
phylumabundance <- justbacteria %>%
  tax_glom(taxrank = "Genus") %>%                     
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>%                                        
  arrange(Phylum) 
head(phylumabundance)

#Summarize the most abundant taxa, filter to only include taxa with > 1% abundance
avg.abundance <- phylumabundance %>%
  select(Phylum, Treatment, Abundance, PlantType, SampleType) %>%
  group_by(Phylum, Treatment, PlantType, SampleType) %>%
  summarize(
    avg_abundance = mean(Abundance), .groups = 'drop'
  ) %>%
  filter(avg_abundance > 0.0003)
head(avg.abundance)
unique(avg.abundance$Phylum)
dim(avg.abundance)

###Save phylum colors
phylum_colors <- c(
  "chocolate1", "cyan3",  "#5F7FC7", "orange","#DA5724", "#508578",
  "#AD6F3B", "#673770", "#652926", "#C84248",  
  "#D1A33D",  "aquamarine3"
)

# Create a factor for SampleType to control the display of x-axis labels
avg.abundance$SampleType <- factor(avg.abundance$SampleType, levels = unique(avg.abundance$SampleType))
avg.abundance$PlantType <- factor(avg.abundance$PlantType, levels = unique(avg.abundance$PlantType))

phylum_avgabundance_sampletype <- ggplot(avg.abundance, aes(x = Treatment, y = avg_abundance, fill = Phylum)) +
  geom_col(position = "fill", show.legend = TRUE) +
  ylab("Proportion of Community") +
  scale_fill_manual(values = phylum_colors) +
  xlab(NULL) +
  facet_grid(PlantType ~ SampleType) +   
  theme(strip.text.x = element_blank()) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 1),
    axis.title.y = element_text(size = 10),
    strip.text.x = element_text(size = 10),  # Adjust facet label size
    strip.placement = "outside",  # Place the facet labels outside the plot area
    title = element_text(size = 10)
  ) +
  labs(x = "Treatment", subtitle = "")  # Adding a subtitle for PlantType

phylum_avgabundance_sampletype

family_colors <- c(
  "darkolivegreen", "chocolate1", "cyan3", "#652926", "gold2", "mediumaquamarine", 
  "#5F7FC7", "orange", "pink3", 
  "#DA5724", "cadetblue", "#508578",
  "#AD6F3B", "#673770", "cadetblue1", "#C84248",  "#D1A33D", 
  "darkolivegreen3", "bisque4", "cyan1", "coral", 
  "darkgoldenrod", "darkorange", "darkseagreen",  "darkslategray3","salmon3",  
  "indianred", "darkseagreen1","sienna1", 
  "chocolate3","aquamarine3", 
  "deepskyblue3", "darkorchid4", "darkseagreen1", "lightsalmon3", 
  "mistyrose3", "tan1", "yellow3",
  "powderblue", "plum4", "peachpuff3", "darkcyan",
  "antiquewhite4", "darkgoldenrod1", "gold4","plum3",
  "darksalmon", "aquamarine","burlywood4", "cornflowerblue",
  "burlywood1", "goldenrod4", "antiquewhite3", "chocolate4", "mediumpurple1")


#Calcultae proportion of each phyla
phylumabundance <- as.data.frame(phylumabundance)
phylumabundance$SampleType <- factor(phylumabundance$SampleType, levels = unique(phylumabundance$SampleType))
phylumabundance$PlantType <- factor(phylumabundance$PlantType, levels = unique(phylumabundance$PlantType))

#e.g. proportion of the community "Acidobacteriota"
Acidobacteriota <- subset(phylumabundance, Phylum =="Acidobacteriota")
Acidobacteriota_avgabundance <- ggplot(Acidobacteriota) +
  geom_col(mapping = aes(x = Treatment, y = Abundance, fill = Genus), position = "fill", show.legend = TRUE) +
  ylab("Proportion of Community") +
  scale_fill_manual(values = phylum_colors) +
  xlab(NULL) +
  facet_grid(PlantType ~ SampleType) +   
  theme(strip.text.x = element_blank()) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 1),
    axis.title.y = element_text(size = 10),
    strip.text.x = element_text(size = 10),  # Adjust facet label size
    strip.placement = "outside",  # Place the facet labels outside the plot area
    title = element_text(size = 10)
  ) +
  labs(x = "Treatment", subtitle = "")  # Adding a subtitle for PlantType


print(Acidobacteriota_avgabundance)

#proportion of the community "Actinobacteriota"
Actinobacteriota <- subset(phylumabundance, Phylum =="Actinobacteriota")
Actinobacteriota_avgabundance <- ggplot(Actinobacteriota) +
  geom_col(mapping = aes(x = Treatment, y = Abundance, fill = Genus), position = "fill", show.legend = TRUE) +
  ylab("Proportion of Community") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL) +
  ggtitle("Actinobacteriota") +  # Add the title here
  theme_minimal() +
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 20))) +  # Position the title at the top and in the middle
  facet_grid(PlantType ~ SampleType) 
print(Actinobacteriota_avgabundance)

#proportion of the community "Firmicutes"
firmicutes <- subset(phylumabundance, Phylum =="Firmicutes")
firmicutes_avgabundance <- ggplot(firmicutes) +
  geom_col(mapping = aes(x = Treatment, y = Abundance, fill = Genus), position = "fill", show.legend = TRUE) +
  ylab("Proportion of Community") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL) +
  ggtitle("Firmicutes") +  # Add the title here
  theme_minimal() +
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 20))) +  # Position the title at the top and in the middle
  facet_grid(PlantType ~ SampleType) 

print(firmicutes_avgabundance)



#proportion of the community "Bacteroidota"
bacteroidota <- subset(phylumabundance, Phylum =="Bacteroidota")
bacteroidota_avgabundance <- ggplot(bacteroidota) +
  geom_col(mapping = aes(x = Treatment, y = Abundance, fill = Genus), position = "fill", show.legend = TRUE) +
  ylab("Proportion of Community") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL) +
  ggtitle("Bacteroidota") +  # Add the title here
  theme_minimal() +
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 20))) +  # Position the title at the top and in the middle
  facet_grid(PlantType ~ SampleType) 

print(bacteroidota_avgabundance)
#proportion of the community "Chloroflexi" 
Chloroflexi <- subset(phylumabundance, Phylum =="Chloroflexi")
Chloroflexi_avgabundance <- ggplot(Chloroflexi) +
  geom_col(mapping = aes(x = Treatment, y = Abundance, fill = Genus), position = "fill", show.legend = TRUE) +
  ylab("Proportion of Community") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL) +
  ggtitle("Chloroflexi") +  # Add the title here
  theme_minimal() +
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 20))) +  # Position the title at the top and in the middle
  facet_grid(PlantType ~ SampleType) 

print(Chloroflexi_avgabundance)

#proportion of the community "Gemmatimonadota"
Gemmatimonadota <- subset(phylumabundance, Phylum =="Gemmatimonadota")
Gemmatimonadota_avgabundance <- ggplot(Gemmatimonadota) +
  geom_col(mapping = aes(x = Treatment, y = Abundance, fill = Genus), position = "fill", show.legend = TRUE) +
  ylab("Proportion of Community") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL) +
  ggtitle("Gemmatimonadota") +  # Add the title here
  theme_minimal() +
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 20))) +  # Position the title at the top and in the middle
  facet_grid(PlantType ~ SampleType) 

print(Gemmatimonadota_avgabundance)

#proportion of the community "Planctomycetota"
Planctomycetota <- subset(phylumabundance, Phylum =="Planctomycetota")
Planctomycetota_avgabundance <- ggplot(Planctomycetota) +
  geom_col(mapping = aes(x = Treatment, y = Abundance, fill = Genus), position = "fill", show.legend = TRUE) +
  ylab("Proportion of Community") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL) +
  ggtitle("Planctomycetota") +  # Add the title here
  theme_minimal() +
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 20))) +  # Position the title at the top and in the middle
  facet_grid(PlantType ~ SampleType) 

print(Planctomycetota_avgabundance)

#proportion of the community "proteobacteria"
proteobacteria <- subset(phylumabundance, Phylum =="Proteobacteria")
proteobacteria_avgabundance <- ggplot(proteobacteria) +
  geom_col(mapping = aes(x = Treatment, y = Abundance, fill = Genus), position = "fill", show.legend = TRUE) +
  ylab("Proportion of Community") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL) +
  ggtitle("proteobacteria") +  # Add the title here
  theme_minimal() +
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 20))) +  # Position the title at the top and in the middle
  facet_grid(PlantType ~ SampleType) 

print(proteobacteria_avgabundance)

#proportion of the community "Verrucomicrobiota"
Verrucomicrobiota <- subset(phylumabundance, Phylum =="Verrucomicrobiota")
Verrucomicrobiota_avgabundance <- ggplot(Verrucomicrobiota) +
  geom_col(mapping = aes(x = Treatment, y = Abundance, fill = Genus), position = "fill", show.legend = TRUE) +
  ylab("Proportion of Community") +
  scale_fill_manual(values = genus_colors) +
  xlab(NULL) +
  ggtitle("Verrucomicrobiota") +  # Add the title here
  theme_minimal() +
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5, margin = margin(b = 20))) +  # Position the title at the top and in the middle
  facet_grid(PlantType ~ SampleType) 

print(Verrucomicrobiota_avgabundance)

