##### Script for making a hierarchical clustering heatmap, the actual heatmap script is adapted from Jessica's existing script

library(RColorBrewer)

#-----------------------Heatmap classes and methylation values -------#

# In my data, GRSet_beta_known is a dataframe of beta values for all the 450k CpG sites across my samples with known subtype
# GRSet_beta_unknown is a dataframe with all the samples of unknown diagnostic subtype
# XXX_class_vector is a vector of sample types matching the order of GRSet_beta_known ie c("XXX", "not XXX", "XXX", "not XXX") 
# Change the subtype name in "Add predicted samples..." below before making heat map.

##### Load data needed #####

load("GRSet_beta_datasets.Rdata")
# GRSet_beta_unknown
# GRSet_beta_known

load("AML_cpg.freqs.new.Rdata") # Info about which CpG sites to include ($index has rownumbers of GRSet_beta_known), need to be subset!
# cpg.freq.MLL.new
# cpg.freq.t821.new
# cpg.freq.inv16.new
# cpg.freq.mono7.new

#load("subtype_annotations_final.Rdata") # Info about which CpG sites to include ($Name)
# cpg.freq.t821_annotation_final
# cpg.freq.inv16_annotation_final
# cpg.freq.mono7_annotation_final
# cpg.freq.MLL_annotation_final

load("AML_class_vectors.Rdata") # Info about classes for the known data
# t821_class_vector
# inv16_class_vector
# mono7_class_vector
# MLL_class_vector

load("subtypes_predicted.Rdata") # Vectors with indexes of the positively predicted samples of the unknown data
# t821.predicted
# inv16.predicted
# mono7.predicted
# MLL.predicted

load("AML_pheno_info.Rdata") # Contains sample ids
# aml.pheno_clean_nonrelapse_known
# aml.pheno_clean_nonrelapse_unknown
# phenoData_clean_nonrelapse_known
# phenoData_clean_nonrelapse_unknown

##### Add the predicted samples to class vector #####
# Do once for each subtype (t821 is omitted because no samples were classified as t821)
# only do this if making class-specific heatmaps

index <- length(inv16_class_vector) + 1
inv16_class_vector_all <- inv16_class_vector # This vector will include the predicted samples as a new class
while (index <= length(inv16_class_vector) + length(inv16.predicted)) { # Adds the class "inv(16) pred" to the class vector
  inv16_class_vector_all[index] <- "inv(16) pred"
  index <- index + 1
}
index <- length(mono7_class_vector) + 1
mono7_class_vector_all <- mono7_class_vector # This vector will include the predicted samples as a new class
while (index <= length(mono7_class_vector) + length(mono7.predicted)) { # Adds the class "mono 7 pred" to the class vector
  mono7_class_vector_all[index] <- "mono 7 pred"
  index <- index + 1
}
index <- length(MLL_class_vector) + 1
MLL_class_vector_all <- MLL_class_vector # This vector will include the predicted samples as a new class
while (index <= length(MLL_class_vector) + length(MLL.predicted)) { # Adds the class "XXX pred" to the class vector
  MLL_class_vector_all[index] <- "MLL pred"
  index <- index + 1
}

##### Add the predicted samples' beta values to GRSet_beta_known #####
# Only do this for class-specific heatmaps

GRSet_beta_result_inv16 <- cbind(GRSet_beta_known, GRSet_beta_unknown[, inv16.predicted]) # inv16.predicted are the new samples found by the classifier
GRSet_beta_result_mono7 <- cbind(GRSet_beta_known, GRSet_beta_unknown[, mono7.predicted]) # mono7.predicted are the new samples found by the classifier
GRSet_beta_result_MLL <- cbind(GRSet_beta_known, GRSet_beta_unknown[, MLL.predicted]) # MLL.predicted are the new samples found by the classifier

# Subset the data to the chosen CpG sites
data.inv16 <- GRSet_beta_result_inv16[cpg.freq.inv16$index,] # Use all chosen CpG sites for the heatmap
data.mono7 <- GRSet_beta_result_mono7[cpg.freq.mono7$index[1:13],] 
data.MLL <- GRSet_beta_result_MLL[cpg.freq.MLL$index[1:13],] 

##### Make vectors with sample ids #####
# Only use this for making class specific heatmaps

id.vector.inv16 <- c(phenoData_clean_nonrelapse_known$Sample_Name, phenoData_clean_nonrelapse_unknown$Sample_Name[inv16.predicted])
id.vector.mono7 <- c(phenoData_clean_nonrelapse_known$Sample_Name, phenoData_clean_nonrelapse_unknown$Sample_Name[mono7.predicted])
id.vector.MLL <- c(phenoData_clean_nonrelapse_known$Sample_Name, phenoData_clean_nonrelapse_unknown$Sample_Name[MLL.predicted])

##### Adding new sample ids to the phenotype files #####

# 77 samples of known subtype. Will call these AML001-AML077
# 58 samples of unknown subtype. Will call these AML078-AML135

# First, make character vector with "AML001" to "AML135"
newids <- data.frame(matrix(ncol = 1, nrow = 135))
colnames(newids) = "newid"
i <- 1
while (i <= 135) {
  if (i < 10) {
    sample <- paste0("AML00", i)
    newids$newid[i] <- sample
    i <- i+1
  } else if (i < 100) {
    sample <- paste0("AML0", i)
    newids$newid[i] <- sample
    i <- i+1
  } else {
    sample <- paste0("AML", i)
    newids$newid[i] <- sample
    i <- i+1
  }
}

# Put both known and unknown samples in the same object
aml.pheno_all <- rbind(aml.pheno_clean_nonrelapse_known, aml.pheno_clean_nonrelapse_unknown)
# Add a column with new sample ids
aml.pheno_newids <- cbind(newids, aml.pheno_all)
rownames(aml.pheno_newids) = newids$newid
# Split the dataframe again into known and unknown samples
aml.pheno_newids_known <- aml.pheno_newids[1:77,]
aml.pheno_newids_unknown <- aml.pheno_newids[78:135,]

save(aml.pheno_newids_known, aml.pheno_newids_unknown, file = "aml_pheno_newids.Rdata")
#load(file = "aml_pheno_newids.Rdata")

##### Make input objects for all data together #####

# The 'data' object with beta values
samples.unknown.heatmap <- c(inv16.predicted, mono7.predicted, MLL.predicted)
# Make id.vector for all samples in one:
id.vector.all <- c(aml.pheno_newids_known$newid, aml.pheno_newids_unknown$newid[samples.unknown.heatmap])
id.vector.known <- aml.pheno_newids_known$newid
GRSet_beta_unknown_heatmap <- GRSet_beta_unknown[, samples.unknown.heatmap] # subset to include only predicted samples
GRSet_beta_all <- cbind(GRSet_beta_known, GRSet_beta_unknown_heatmap) # add known and predicted unknown data together
cpgs.for.heatmap <- c(cpg.freq.t821.new$index, cpg.freq.MLL.new$index[1:13], cpg.freq.inv16.new$index[1:6], cpg.freq.mono7.new$index)

data.all <- GRSet_beta_all[cpgs.for.heatmap, ] # beta values for the chosen cpg sites and samples
data.known <- GRSet_beta_known[cpgs.for.heatmap, ]

# The 'class vector' object with class info

load(file = "AML_clean_nonrelapse.Rdata")

unknown_indices <- which(aml.pheno_clean_nonrelapse$genotype %in% c("normal", "other clon abn", "no result"))
subtype_vector_known <- aml.pheno_clean_nonrelapse$genotype[-unknown_indices]
class_vector_known <- as.character(subtype_vector_known) # make character instead of factor, for replacement
class_vector_known[class_vector_known %in% c("other 11q23/MLL", "t(9;11)", "t(10;11)", "t(11;19)")] <- "MLL"
class_vector_known[!class_vector_known %in% c("MLL", "t(8;21)", "inv(16)", "mono 7")] <- "other"

class_vector_all <- class_vector_known # This vector will include the predicted samples as new classes
index <- length(class_vector_all) + 1
while (index <= length(class_vector_known) + length(inv16.predicted)) { # Adds the class "inv(16) pred" to the class vector
  class_vector_all[index] <- "inv(16) pred"
  index <- index + 1
}
n <- length(class_vector_all)
index <- length(class_vector_all) + 1
while (index <= n + length(mono7.predicted)) { # Adds the class "mono 7 pred" to the class vector
  class_vector_all[index] <- "mono 7 pred"
  index <- index + 1
}
n <- length(class_vector_all)
index <- length(class_vector_all) + 1
while (index <= n + length(MLL.predicted)) { # Adds the class "MLL pred" to the class vector
  class_vector_all[index] <- "MLL pred"
  index <- index + 1
} 
class_vector_all_other <- class_vector_all
# class_vector_all_other now has all the classes needed for the heatmap and the rest named "other"

# Also subset the pheno objects accordingly (3 inv16, 4 mono7, 7 MLL)
phenoData_heatmap <- phenoData_clean_nonrelapse_unknown[samples.unknown.heatmap, ]
aml.info_heatmap <- aml.info_clean_nonrelapse_unknown[samples.unknown.heatmap, ]
aml.pheno_heatmap <- aml.pheno_newids_unknown[samples.unknown.heatmap, ]

aml.pheno_heatmap.all <- rbind(aml.pheno_newids_known, aml.pheno_heatmap)
aml.info_heatmap.all <- rbind(aml.info_clean_nonrelapse_known, aml.info_heatmap)
phenoData_heatmap.all <- rbind(phenoData_clean_nonrelapse_known, phenoData_heatmap)

##### Make variables for _only_ samples with t821, inv16, mono7, MLL + predicted unknown samples #####

# Removal of known samples of other classes from the "known" variables is performed in Data_for_pamr.R
# Two sets of input will be made: .4classes.known and .4classes.all, the latter including the predicted samples
# First, load the R object with the subset beta values and class vectors:

load(file = "for_heatmaps_4classes.Rdata")

# Subset the pheno file (including new ids AMLXXX):
aml.pheno_newids_known_4classes <- aml.pheno_newids_known[indices_4classes, ]
aml.pheno_newids_all_4classes <- rbind(aml.pheno_newids_known[indices_4classes, ], aml.pheno_newids_unknown[samples.unknown.heatmap, ])

# Now make the input arguments for heatmaps with 1. only known samples and 2. predicted samples added

#samples.unknown.heatmap are the indices of the predicted samples
# Make id.vectors:
id.vector.all.4classes <- c(aml.pheno_newids_known_4classes$newid, aml.pheno_newids_unknown$newid[samples.unknown.heatmap])
id.vector.known.4classes <- aml.pheno_newids_known_4classes$newid
# Make beta value input objects :
#GRSet_beta_unknown_heatmap includes only predicted samples
GRSet_beta_all_4classes <- cbind(GRSet_beta_known_4classes, GRSet_beta_unknown_heatmap) # add known and predicted unknown data together
#cpgs.for.heatmap are the CpG sites chosen for the four final classifiers

data.all.4classes <- GRSet_beta_all_4classes[cpgs.for.heatmap, ] # beta values for the chosen cpg sites and samples
data.known.4classes <- GRSet_beta_known_4classes[cpgs.for.heatmap, ]

# The 'class vector' object with class info

unknown_indices <- which(aml.pheno_clean_nonrelapse$genotype %in% c("normal", "other clon abn", "no result"))
subtype_vector_known <- aml.pheno_clean_nonrelapse$genotype[-unknown_indices]
class_vector_known <- as.character(subtype_vector_known) # make character instead of factor, for replacement
class_vector_known[class_vector_known %in% c("other 11q23/MLL", "t(9;11)", "t(10;11)", "t(11;19)")] <- "MLL"
class_vector_known[!class_vector_known %in% c("MLL", "t(8;21)", "inv(16)", "mono 7")] <- "other"

# class_vector_known has the four classes + class "other"
# remove the "other" class using the indices4_classes vector:
class_vector_known_4classes <- class_vector_known[indices_4classes]

# Add the predicted samples as classes to the class vector:
class_vector_all <- class_vector_known_4classes # This vector will include the predicted samples as new classes
index <- length(class_vector_all) + 1
while (index <= length(class_vector_known_4classes) + length(inv16.predicted)) { # Adds the class "inv(16) pred" to the class vector
  class_vector_all[index] <- "inv(16) pred"
  index <- index + 1
}
n <- length(class_vector_all)
index <- length(class_vector_all) + 1
while (index <= n + length(mono7.predicted)) { # Adds the class "mono 7 pred" to the class vector
  class_vector_all[index] <- "mono 7 pred"
  index <- index + 1
}
n <- length(class_vector_all)
index <- length(class_vector_all) + 1
while (index <= n + length(MLL.predicted)) { # Adds the class "MLL pred" to the class vector
  class_vector_all[index] <- "MLL pred"
  index <- index + 1
} 
class_vector_all_4classes <- class_vector_all
# class_vector_all_4classes now has all the classes needed for the heatmaps with only 4 classes

##### Order CpG sites in class order #####
# To keep class specific CpG sites together in the heatmap.

# Indices of CpG sites are in cpg.freq.XXX$index

t821.cpgs <- cpg.freq.t821.new$index
inv16.cpgs <- cpg.freq.inv16.new$index[1:6]
mono7.cpgs <- cpg.freq.mono7.new$index
MLL.cpgs <- cpg.freq.MLL.new$index[1:13]

# Cluster the CpG sites class-wise
# First, make function that takes a dataframe with beta values, clusters the data, and returns the order of the CpG sites

clust.cpgs <- function(data) {
  d <- as.dist(1 - cor(data, use='pairwise.complete.obs'))
  hca <- hclust(d, method='ward.D2')
  use.data.t <- t(data)
  d2 <- as.dist(1 - cor(use.data.t, use='pairwise.complete.obs'))
  hca2 <- hclust(d2, method='ward.D2')
  return(hca2$order)
}

t821.cpgs.clust <- clust.cpgs(GRSet_beta_known[t821.cpgs, ]) # t821 CpG site indices in clustered order
inv16.cpgs.clust <- clust.cpgs(GRSet_beta_known[inv16.cpgs, ]) # inv16 CpG site indices in clustered order
mono7.cpgs.clust <- clust.cpgs(GRSet_beta_known[mono7.cpgs, ]) # mono7 CpG site indices in clustered order
MLL.cpgs.clust <- clust.cpgs(GRSet_beta_known[MLL.cpgs, ]) # MLL CpG site indices in clustered order

# Reorder cpg sites in the cpg.freq objects created above, then merge them into one object

t821.cpgs.ordered <- t821.cpgs[t821.cpgs.clust]
inv16.cpgs.ordered <- inv16.cpgs[inv16.cpgs.clust]
mono7.cpgs.ordered <- mono7.cpgs[mono7.cpgs.clust]
MLL.cpgs.ordered <- MLL.cpgs[MLL.cpgs.clust]

cpgs_class_ordered <- c(t821.cpgs.ordered, inv16.cpgs.ordered, mono7.cpgs.ordered, MLL.cpgs.ordered) # Use this vector for ordering the cpgs in the heatmaps for figures
data.known.4classes.ordered <- GRSet_beta_known_4classes[cpgs_class_ordered, ] #reordering data input for heatmap
data.known.ordered <- GRSet_beta_known[cpgs_class_ordered, ]
data.all.4classes.ordered <- GRSet_beta_all_4classes[cpgs_class_ordered, ]
data.all.ordered <- GRSet_beta_all[cpgs_class_ordered, ]

# Make "class vector" for class specific CpG sites
cpg.class.vector <- c(rep("t(8;21)", length(t821.cpgs)), rep("inv(16)", length(inv16.cpgs)), rep("mono 7", length(mono7.cpgs)), rep("MLL", length(MLL.cpgs)))

# Save the four sets of beta values:
save(GRSet_beta_all, GRSet_beta_all_4classes, GRSet_beta_known, GRSet_beta_known_4classes, file = "GRSet_beta_forheatmaps.Rdata")

#############################################################################################
##### Set variables for the script below #####
# Do once for each run, run heatmap script below after changing the pdf file name of previous run

data <- data.known.ordered # data.known.ordered, data.all.ordered, data.known.4classes.ordered or data.all.4classes.ordered
class_vector <- class_vector_known #class_vector_known, class_vector_all_other, class_vector_known_4classes or class_vector_all_4classes
id.vector <- id.vector.known # id.vector.known, id.vector.all, id.vector.known.4classes or id.vector.all.4classes
# There is also a cpg.class.vector for use when labeling class specific CpG sites

#-----------------------Make colors ----------------------------------#
# For samples:
GK <- as.character(class_vector)
col.idxG <- rep(NA, length(GK))
ugk <- unique(GK[!is.na(GK)])
ugk <- sort(ugk)
for (i in 1:length(ugk)) {
  col.idxG[GK %in% ugk[i]] <- i
}
col.idxG[is.na(col.idxG)] <- length(ugk) + 1
# Below, use the set of colors corresponding to the input. row 1: all, row 2: all.4classes, row 3: known, row 4: known.4classes
#my.colsG <-c("red", "orange", "green3", "palegreen", "purple", "blue", "pink", "powderblue", "yellow", gray(0))#Need as many colors as classes plus gray at the end
#my.colsG <-c("red", "orange", "green3", "palegreen", "purple", "blue", "powderblue", "yellow", gray(0))
my.colsG <-c("red", "green3", "purple", "pink", "powderblue", "yellow", gray(0))
#my.colsG <-c("red", "green3", "purple", "powderblue", "yellow", gray(0))

# For CpG sites:
GK <- as.character(cpg.class.vector)
row.idxG <- rep(NA, length(GK))
cpg.ugk <- unique(GK[!is.na(GK)])
cpg.ugk <- sort(cpg.ugk)
for (i in 1:length(cpg.ugk)) {
  row.idxG[GK %in% cpg.ugk[i]] <- i
}
row.idxG[is.na(row.idxG)] <- length(cpg.ugk) + 1
# Use the following colors always, independent of n sample classes:
my.rowsG <-c("red", "green3", "purple", "powderblue", "yellow", gray(0))

#----------------------Distance and clustering ------------------------#
d <- as.dist(1 - cor(data, use='pairwise.complete.obs'))
hca <- hclust(d, method='ward.D2')
#use.data.t <- t(data) #comment this row if not clustering the CpG sites
#d2 <- as.dist(1 - cor(use.data.t, use='pairwise.complete.obs')) #comment this row if not clustering the CpG sites
#hca2 <- hclust(d2, method='ward.D2') #comment this row if not clustering the CpG sites

#-----------------------MAKE HEATMAP ----------------------------------#
pdf('AML_heatmap.pdf', paper='special', width=6, height=3,5)
par(mar=c(0,0,0,0), omi=c(0.3,0,0,0), ps=8)

# To label samples, enter a vector  for labCol with the names in the same order as in data, use id.vector
# Use the first heatmap() input for NOT clustering CpG sites, no dendrogram. Also, this one colors the CpG sites.
# Use the second heatmap() input for ordinary clustering of CpG sites

heatmap(data, Colv=as.dendrogram(hca, hang=-1), Rowv=NA,
        ColSideColors=my.colsG[col.idxG], RowSideColors=my.rowsG[row.idxG], labCol= id.vector, labRow= '', 
        scale="none", col= rev(brewer.pal(10, 'RdYlBu')),
        margins=c(1,1), cexCol=1, cexRow = 1, useRaster=T)

#heatmap(data, Colv=as.dendrogram(hca, hang=-1), Rowv=as.dendrogram(hca2, hang=-1), 
#        ColSideColors=my.colsG[col.idxG], labCol= id.vector, labRow= '', 
#        scale="none", col= rev(brewer.pal(10, 'RdYlBu')),
#        margins=c(1,1), cexCol=1, cexRow = 1, useRaster=T)

lc1 <- legend('topleft', fill=my.colsG[-length(my.colsG)], legend=ugk, 
              bty='n', ncol=1)
dev.off()

# Add the subtype to the pdf file name!!!

#------------------MAKE intensity scale-------------------#
# As a separate pdf file
# Do once for all four heatmaps

data <- data.all.ordered # data.known, data.all, data.known.4classes or data.all.4classes
my.colors <- rev(brewer.pal(10, 'RdYlBu')) # same as "col" argument in heatmap() above

pdf('intensity_scale.pdf')
draw.scale(data, rev(my.colors)) # draw.scale() doesn't exist...
dev.off()

#--------------------------------------

##### Make correlation plot #####

library(corrplot)

my.cor <- round(cor(data, use="pairwise.complete.obs", method="pearson"), 2) 
my.names <- class_vector_all_4classes #vector of class names, for example 
  #class_vector_known, class_vector_all_other, class_vector_known_4classes or class_vector_all_4classes
row.names(my.cor) <- my.names 
colnames(my.cor) <- my.names 
pdf(file="plots/Corrplot_2016.pdf", height=7, width=7) 
corrplot.mixed(my.cor, is.corr=FALSE, order="AOE", tl.pos="lt", diag="l")


##### Making csv files for analysis of sample clustering #####

# id vectors:

#id.vector.known
#id.vector.known.4classes
#id.vector.all
#id.vector.all.4classes

# aml.pheno objects (with anonymous ids):

#aml.pheno_newids_known
#aml.pheno_newids_known_4classes
#aml.pheno_heatmap.all
#aml.pheno_newids_all_4classes

# heatmap orders, to reorder id vectors and pheno files

data <- data.all.ordered # Do for data.known.ordered, data.known.4classes.ordered, data.all.ordered, data.all.4classes.ordered
d <- as.dist(1 - cor(data, use='pairwise.complete.obs'))
hca <- hclust(d, method='ward.D2')
heatmap.order.all <- hca$order # change variable name for each heatmap

# Sample orders are in 
#heatmap.order.known.4classes
#heatmap.order.known
#heatmap.order.all
#heatmap.order.all.4classes

id.vector.all.4classes[heatmap.order.all.4classes]
rownames(aml.pheno_newids_all_4classes[heatmap.order.all.4classes,])
# These 2 should look the same (do for all 4 heatmaps)

# Add class vectors as a new column in aml.pheno objects with predicted samples (ordering is done later):
aml.pheno_newids_all_classinfo <- cbind(aml.pheno_heatmap.all, class_vector_all_other)
aml.pheno_newids_all_4classes_classinfo <- cbind(aml.pheno_newids_all_4classes, class_vector_all_4classes)

# Write out the new aml.pheno object with the samples in the correct order, to a csv file
write.csv(aml.pheno_newids_known[heatmap.order.known, ], file = "aml_pheno_heatmap_known.csv")
write.csv(aml.pheno_newids_known_4classes[heatmap.order.known.4classes, ], file = "aml_pheno_heatmap_known_4classes.csv")
write.csv(aml.pheno_newids_all_classinfo[heatmap.order.all, ], file = "aml_pheno_heatmap_all.csv")
write.csv(aml.pheno_newids_all_4classes_classinfo[heatmap.order.all.4classes, ], file = "aml_pheno_heatmap_all_4classes.csv")
