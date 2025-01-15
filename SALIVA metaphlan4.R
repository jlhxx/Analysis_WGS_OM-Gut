# Perform microbiome data analysis using the phyloseq package in R.
setwd("dir...")

library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library("microbiome") # data analysis and visualisation
library("phyloseq") # also the basis of data object. Data analysis and visualisation
library("microbiomeutilities") # some utility tools 
library("RColorBrewer") # nice color options
library("ggpubr") # publication quality figures, based on ggplot2
library("DT") # interactive tables in html and markdown
library("data.table") # alternative to data.frame
library("dplyr") # data handling  
library("vegan")


#IMPORT DATA TO PHYLOSEQ
otu_mat <- read_excel("gut_saliva_16.xlsx", sheet = "saliva_otu")
tax_mat<- read_excel("gut_saliva_16.xlsx", sheet = "saliva_taxonomy")
samples_df <- read_excel("gut_saliva_16.xlsx", sheet = "saliva_sample")

#make otu mat numeric
otu_mat<-as.data.frame(otu_mat)
rownames(otu_mat)<-otu_mat[,1]
otu_mat<-otu_mat[,-1]
a <- rownames(otu_mat)
otu_mat <- as.data.frame(lapply(otu_mat, as.numeric))
rownames(otu_mat) <- a # Reassign rownames

#Phyloseq objects need to have row.names
#define the row names from the otu column
#otu_mat <- otu_mat %>%
#  tibble::column_to_rownames("otu") 

#Idem for the two other matrixes
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("otu")

samples_df <- samples_df %>% 
  tibble::column_to_rownames("sample") 

#Transform into matrixes otu and tax tables (sample table can be left as data frame)

otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

#Transform to phyloseq objects

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

#make phyloseq obj
Saliva <- phyloseq(OTU, TAX, samples)
###########################################################
#REMOVE GGB if necessary
tax_table_saliva<-tax_table(Saliva)
# Extract the taxa names for species
species_names <- tax_table_saliva[, "Species"]

# Identify the taxa that start with "s__GGB"
taxa_to_remove <- taxa_names(Saliva)[species_names %in% grep("^s__GGB", species_names, value = TRUE)]

# Prune the identified taxa
saliva_filtered <- prune_taxa(!taxa_names(Saliva) %in% taxa_to_remove, Saliva)
# Check the number of taxa before and after filtering
cat("Number of taxa before filtering:", ntaxa(Saliva), "\n")
cat("Number of taxa after filtering:", ntaxa(saliva_filtered), "\n")

################################################

# Calculate the number of samples
sample_count <- nsamples(saliva_filtered)

# Determine the minimum number of samples an OTU must be present in
min_samples <- ceiling(0.05 * sample_count)

# Define a prevalence filter function
prevalence_filter <- function(otu) {
  return(sum(otu > 0) >= min_samples)
}

# Apply the prevalence filter
saliva_filtered_prevalence <- prune_taxa(apply(otu_table(saliva_filtered), 1, prevalence_filter), saliva_filtered)

# Convert to relative abundance
saliva_rel_abundance <- transform_sample_counts(saliva_filtered_prevalence, function(x) x / sum(x))

# Define an abundance filter function
abundance_filter <- function(otu) {
  return(any(otu >= 0.01))
}

# Apply the abundance filter
saliva_filtered_abundance <- prune_taxa(apply(otu_table(saliva_rel_abundance), 1, abundance_filter), saliva_filtered_prevalence)

# The final filtered phyloseq object
saliva_filtered_abundance
#write.table(otu_table(saliva_filtered_abundance), file="clipboard-16384", sep="\t", col.names=NA)
#write.table(tax_table(saliva_filtered_abundance), file="clipboard-16384", sep="\t", col.names=NA)


########################################################
#remove any OTUs that present only one time (singletons).
#Saliva.prune = prune_taxa(taxa_sums(Saliva) > 1, Saliva)
Saliva.prune<-saliva_filtered_abundance
#Let’s check whether our data are pruned or not.
Saliva
Saliva.prune
ps0.rar <- Saliva.prune

# Check read count
readcount = data.table(as(sample_data(Saliva.prune), "data.frame"),
                       TotalReads = sample_sums(Saliva.prune), 
                       keep.rownames = TRUE)
setnames(readcount, "rn", "SampleID")
ggplot(readcount, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
#see the distribution of our samples readcounts.

#check samples with low number of reads, “order()” can be used to sort “TotalReads” column.
head(readcount[order(readcount$TotalReads), c("SampleID", "TotalReads")])

#Generate rarefaction curve, rarefaction curve could be used to determined whether the sequencing depth cover microbial diversity of the sample.
#only work for count data
#otu.rare = otu_table(Saliva.prune)
#otu.rare = as.data.frame(t(otu.rare))
#sample_names = rownames(otu.rare)

# we will use vegan rarecurve #only for counts
#otu.rarecurve = rarecurve(otu.rare, step = 10000, label = T)

############################################

set.seed(9242)  #reproducing the filtering and nomalisation. 
summary(sample_sums(Saliva.prune))

#barplot(sample_sums(ps0.rar), las =2)
# quick check taxa prevalence
#barplot(sample_sums(ps0.rar), las =2)
# quick check taxa prevalence
ps0.rar <- Saliva.prune

p.rar <- plot_taxa_prevalence(ps0.rar, "Phylum")

p.rar

#Non-phylogenetic diversities
hmp.div <- alpha(ps0.rar, index = "all")
#datatable(hmp.div)

# get the metadata
hmp.meta <- meta(ps0.rar)

# Add the rownames as a new colum for easy integration later.
hmp.meta$sam_name <- rownames(hmp.meta)

# Add the rownames to diversity table
hmp.div$sam_name <- rownames(hmp.div)

# merge these two data frames into one
div.df <- merge(hmp.div,hmp.meta, by = "sam_name")

# check the tables
colnames(div.df)

###############################################################
#ALTERNATIVE WAY
colnames(hmp.div)

# convert phyloseq object into a long data format.

div.df2 <- div.df[, c("Group", "diversity_shannon", "chao1")]

# the names are not pretty. we can replace them

colnames(div.df2) <- c("Group", "Shannon", "chao1")

# check
colnames(div.df2)

div_df_melt <- reshape2::melt(div.df2)
head(div_df_melt)
# Now use this data frame to plot 
p <- ggboxplot(div_df_melt, x = "Group", y = "value",
               fill = "Group", 
               palette = "jco", 
               legend= "right",
               facet.by = "variable", 
               scales = "free")

p <- p + rotate_x_text() 
# we will remove the x axis lables

p <- p + rremove("x.text")
p

div_df_melt$Group<-as.factor(div_df_melt$Group)
lev <- levels(div_df_melt$Group) # get the variables

# make a pairwise list that we want to compare.
L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])

pval <- list(
  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
  symbols = c("****", "***", "**", "*", "n.s")
)

p2 <- p + stat_compare_means(
  comparisons = L.pairs,
  paired = TRUE,
  na.rm = TRUE,
  label = "p.signif",
  symnum.args = list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("****", "***", "**", "*", "n.s")
  )
)

print(p2)
#downlaod raw alpha datatable
write.table(div_df_melt, file="clipboard-16384", sep="\t", col.names=NA)

#############################################

#Barplot counts
ps1.com <- saliva_filtered_abundance

# if you have dada2/deblur output and sequences as taxa names, then you can change them as follows
taxa_names(ps1.com) <- paste0("ASV_", rownames(tax_table(ps1.com)))

#Set Palette
taxic <- as.data.frame(ps1.com@tax_table) # this will help in setting large color options

# colourCount = length(unique(taxic$Species))  #define number of variable colors based on number of Species (change the level accordingly to Species/Species/Species)
# getPalette = colorRampPalette(brewer.pal(12, "Paired"))  # change the palette as well as the number of colors will change according to palette.

taxic$OTU <- rownames(taxic) # Add the OTU/ASV ids from OTU table into the taxa table at the end.
colnames(taxic) # now have extra taxonomy levels.

taxmat <- as.matrix(taxic) # convert it into a matrix.
new.tax <- tax_table(taxmat) # convert into phyloseq compatible file.
tax_table(ps1.com) <- new.tax # incroporate into phyloseq Object


# Edit the unSpeciesified taxa
tax_table(ps1.com)[tax_table(ps1.com)[, "Species"] == "", "Species"] <- "UnSpeciesified Species"

# Taxonomic names in italics

guide_italics <- guides(fill = guide_legend(label.theme = element_text(
  size = 15,
  face = "italic", colour = "Black", angle = 0
)))


## Plot at Species level - select top 10 highest rel abundance
ps1.com@phy_tree <- NULL

# Merge at Species level
ps1.com.fam <- microbiomeutilities::aggregate_top_taxa2(ps1.com, "Species", top = 10)

plot.composition.COuntAbun <- plot_composition(ps1.com.fam) + theme(legend.position = "bottom") +
  scale_fill_brewer("Species", palette = "Paired") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Relative abundance") + guide_italics + theme(legend.title = element_text(size = 18))

plot.composition.COuntAbun

#################
# Barplot relative abundance

# the previous pseq object ps1.com.fam is only counts.

# Use transform function of microbiome package to convert it to rel abun.
ps1.com.fam.rel <- microbiome::transform(ps1.com.fam, "compositional")

plot.composition.relAbun <- plot_composition(ps1.com.fam.rel,
                                             sample.sort = "Group",
                                             x.label = "env_material") 
plot.composition.relAbun <- plot.composition.relAbun + theme(legend.position = "bottom") 
plot.composition.relAbun <- plot.composition.relAbun + scale_fill_brewer("Species", palette = "Paired") + theme_bw() 
plot.composition.relAbun <- plot.composition.relAbun + theme(axis.text.x = element_text(angle = 90)) 
plot.composition.relAbun <- plot.composition.relAbun + ggtitle("Relative abundance") + guide_italics + theme(legend.title = element_text(size = 18))

print(plot.composition.relAbun)
write.table(otu_table(ps1.com.fam), file="clipboard-16384", sep="\t", col.names=NA)

##########################################################

#Heatmaps
# base plot
data.com <- plot.composition.relAbun$data
colnames(data.com)

p.heat <- ggplot(data.com, aes(x = Sample, y = Tax)) + geom_tile(aes(fill = Abundance)) 

# Change color
p.heat <- p.heat + scale_fill_distiller("Abundance", palette = "RdYlBu") + theme_bw() 

# Make bacterial names italics
p.heat <- p.heat + theme(axis.text.y = element_text(colour = 'black', 
                                                    size = 10, 
                                                    face = 'italic')) 
# Make seperate samples based on main varaible
p.heat <- p.heat + facet_grid(~xlabel, 
                              scales = "free") + rremove("x.text") 

p.heat <- p.heat + ylab("Species")

#Clean the x-axis
p.heat <- p.heat + theme(axis.title.x=element_blank(),
                         axis.text.x=element_blank(),
                         axis.ticks.x=element_blank()) 

# Clean the facet label box
p.heat <- p.heat + theme(legend.key = element_blank(), 
                         strip.background = element_rect(colour="black", fill="white"))

print(p.heat)#
###########################################
# Manually add the 'Group' column based on your criteria of metadata
data.com$Group <- NA
data.com$Group[data.com$Sample %in% c("S1_baseline_A",	"S3_baseline_A",	"S9_baseline_A",	"S11_baseline_A",	"S12_baseline_A",	"S14_baseline_A",	"S16_baseline_A",	"S18_baseline_A",	"S21_baseline_A",	"S22_baseline_A",	"S26_baseline_A",	"S28_baseline_A",	"S29_baseline_A",	"S30_baseline_A",	"S31_baseline_A",	"S32_baseline_A",	"S10_baseline_A")] <- "baseline A"
data.com$Group[data.com$Sample %in% c("S1_baseline_B",	"S3_baseline_B",	"S10_baseline_B",	"S11_baseline_B",	"S12_baseline_B",	"S14_baseline_B",	"S16_baseline_B",	"S18_baseline_B",	"S21_baseline_B",	"S22_baseline_B",	"S26_baseline_B",	"S28_baseline_B",	"S29_baseline_B",	"S30_baseline_B",	"S31_baseline_B",	"S32_baseline_B",	"S9_baseline_B")] <- "baseline B"
data.com$Group[data.com$Sample %in% c("S1_post_A",	"S3_post_A",	"S10_post_A",	"S11_post_A",	"S12_post_A",	"S14_post_A",	"S16_post_A",	"S18_post_A",	"S21_post_A",	"S22_post_A",	"S26_post_A",	"S28_post_A",	"S29_post_A",	"S30_post_A",	"S31_post_A",	"S32_post_A",	"S9_post_A")] <- "post A"
data.com$Group[data.com$Sample %in% c("S1_post_B",	"S3_post_B",	"S9_post_B",	"S10_post_B",	"S11_post_B",	"S12_post_B",	"S14_post_B",	"S16_post_B",	"S18_post_B",	"S21_post_B",	"S22_post_B",	"S26_post_B",	"S28_post_B",	"S30_post_B",	"S29_post_B",	"S31_post_B",	"S32_post_B")] <- "post B"

# Create heatmap with ggplot2
p.heat <- ggplot(data.com, aes(x = Sample, y = Tax)) + 
  geom_tile(aes(fill = Abundance)) + 
  scale_fill_distiller("Abundance", palette = "RdYlBu") + 
  theme_bw() + 
  theme(axis.text.y = element_text(colour = 'black', 
                                   size = 10, 
                                   face = 'italic')) + 
  facet_wrap(~ Group, scales = "free_x") + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) + 
  ylab("Species") + 
  theme(legend.key = element_blank(), 
        strip.background = element_rect(colour = "black", 
                                        fill = "white"))

print(p.heat)
##########################################
# we use count data at Species level from the barplot for counts

ps_df <- microbiomeutilities::phy_to_ldf(ps1.com.fam, 
                                         transform.counts = "compositional")

colnames(ps_df)
# this data.frame can be used to customize several plots.  

# example boxplot at Species level

p.box <- ggstripchart(ps_df, "Group", "Abundance", 
                      facet.by = "Species", color = "Group",
                      palette = "jco"
)

p.box + rremove("x.text")
#saveRDS(object, file = "my_data.rds")

##########################
#Beta diversity metrics
#Beta-diversity: Measures for differences between samples from different groups to identify if there are differences in the overall community composition and structure.
#Unweighted Unifrac
#Unweighted Unifrac is based on presence/absence of different taxa and abundance is not important.
# it is sensitive to the sequencing depth. 
# if we remove OTUs that are detected at least 10 times in 5% of the samples
#saliva.prune not rareified
#ps0.rar is rareified

ps0.rar.filtered <- core(ps0.rar, detection = 10, prevalence = 0.05)

summarize_phyloseq(ps0.rar.filtered)
ps1 <- Saliva.prune

ps1.rel <- microbiome::transform(ps1, "compositional")

#Population-level Density landscapes
p <- plot_landscape(ps1.rel, 
                    "NMDS", 
                    "bray", 
                    size = 3,
                    col = "Group") +
  labs(title = paste("NMDS / Bray-Curtis"))   

p <- p + scale_color_brewer(palette = "Dark2")+ scale_fill_gradient(low = "#e0ecf4", high = "orange") 
p 
##################################
#Core anlaysis
# convert to relative abundance  

ps1.stool.rel <- microbiome::transform(ps1, "compositional")
print(ps1)
ps1.stool.rel2 <- prune_taxa(taxa_sums(ps1.stool.rel) > 0, ps1.stool.rel)

print(ps1.stool.rel2)

#find core otus
core.taxa.standard <- core_members(ps1.stool.rel2, detection = 0.001, prevalence = 50/100)
print(core.taxa.standard)

# Extract the taxonomy table
taxonomy <- as.data.frame(tax_table(ps1.stool.rel2))

# Subset this taxonomy table to include only core OTUs  
core_taxa_id <- subset(taxonomy, rownames(taxonomy) %in% core.taxa.standard)

DT::datatable(core_taxa_id)

###################################################################
#Core abundance and diversity
core.abundance <- sample_sums(core(ps1.stool.rel2, detection = 0.001, prevalence = 50/100))

DT::datatable(as.data.frame(core.abundance))

#Core visualization
#Core heatmaps
# Core with compositionals:
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)

# define gray color palette
gray <- gray(seq(0,1,length=5))
p.core <- plot_core(ps1.stool.rel2, 
                    plot.type = "heatmap", 
                    colours = gray,
                    prevalences = prevalences, 
                    detections = detections, 
                    min.prevalence = .5) +
  xlab("Detection Threshold (Relative Abundance (%))")
print(p.core)   
                 
# Same with the viridis color palette
# color-blind friendly and uniform
# options: viridis, magma, plasma, inferno
# https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html
# Also discrete=TRUE versions available
library(viridis)
print(p.core + scale_fill_viridis())

# Core with compositionals:
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)

# Also define gray color palette

p.core <- plot_core(ps1.stool.rel2, 
                    plot.type = "heatmap", 
                    colours = rev(brewer.pal(5, "Spectral")),
                    prevalences = prevalences, 
                    detections = detections, 
                    min.prevalence = .5) + 
  xlab("Detection Threshold (Relative Abundance (%))")

print(p.core) 
ps1.stool.rel2.f <- microbiomeutilities::format_to_besthit(ps1.stool.rel2)
p.core <- plot_core(ps1.stool.rel2.f, 
                    plot.type = "heatmap", 
                    colours = rev(brewer.pal(5, "Spectral")),
                    prevalences = prevalences, 
                    detections = detections, 
                    min.prevalence = .5) + 
  xlab("Detection Threshold (Relative Abundance (%))")

p.core + theme(axis.text.y = element_text(face="italic"))
print(p.core)
