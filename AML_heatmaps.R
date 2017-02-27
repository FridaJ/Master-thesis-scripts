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

load("AML_cpg.freqs.Rdata") # Info about which CpG sites to include ($index has rownumbers of GRSet_beta_known), need to be subset!
# cpg.freq.MLL
# cpg.freq.t821
# cpg.freq.inv16
# cpg.freq.mono7

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

GRSet_beta_result_inv16 <- cbind(GRSet_beta_known, GRSet_beta_unknown[, inv16.predicted]) # inv16.predicted are the new samples found by the classifier
GRSet_beta_result_mono7 <- cbind(GRSet_beta_known, GRSet_beta_unknown[, mono7.predicted]) # mono7.predicted are the new samples found by the classifier
GRSet_beta_result_MLL <- cbind(GRSet_beta_known, GRSet_beta_unknown[, MLL.predicted]) # MLL.predicted are the new samples found by the classifier

# Subset the data to the chosen CpG sites
data.inv16 <- GRSet_beta_result_inv16[cpg.freq.inv16$index,] # Use the chosen CpG sites for the heatmap
data.mono7 <- GRSet_beta_result_mono7[cpg.freq.mono7$index[1:13],] 
data.MLL <- GRSet_beta_result_MLL[cpg.freq.MLL$index[1:13],] 

##### Make vectors with sample ids #####

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

##### Make input objects for all data together #####

# The 'data' object with beta values
samples.unknown.heatmap <- c(inv16.predicted, mono7.predicted, MLL.predicted)
# Make id.vector for all samples in one:
id.vector.all <- c(aml.pheno_newids_known$newid, aml.pheno_newids_unknown$newid[samples.unknown.heatmap])
id.vector.known <- aml.pheno_newids_known$newid
GRSet_beta_unknown_heatmap <- GRSet_beta_unknown[, samples.unknown.heatmap] # subset to include only predicted samples
GRSet_beta_all <- cbind(GRSet_beta_known, GRSet_beta_unknown_heatmap) # add known and predicted unknown data together
cpgs.for.heatmap <- c(cpg.freq.MLL$index[1:13], cpg.freq.inv16$index, cpg.freq.mono7$index[1:13])
GRSet_beta_heatmap <- GRSet_beta_all[cpgs.for.heatmap, ] # beta values with the chosen cpg sites and samples

data.all <- GRSet_beta_heatmap
data.known <- GRSet_beta_known[cpgs.for.heatmap, ]

# The 'class vector' object with class info

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
# class_vector_all now has all the classes needed for the heatmap and the rest named "other"

# Also subset the pheno objects accordingly (3 inv16, 4 mono7, 7 MLL)
phenoData_heatmap <- phenoData_clean_nonrelapse_unknown[samples.unknown.heatmap, ]
aml.info_heatmap <- aml.info_clean_nonrelapse_unknown[samples.unknown.heatmap, ]
aml.pheno_heatmap <- aml.pheno_newids_unknown[samples.unknown.heatmap, ]

aml.pheno_heatmap.all <- rbind(aml.pheno_newids_known, aml.pheno_heatmap)
aml.info_heatmap.all <- rbind(aml.info_clean_nonrelapse_known, aml.info_heatmap)
phenoData_heatmap.all <- rbind(phenoData_clean_nonrelapse_known, phenoData_heatmap)

##### Set variables for the script below #####
# Do once for each subtype, run heatmap script below after changing the pdf file name of previous run

data <- data.all # data.known or data.all
class_vector <- class_vector_all #class_vector_known or class_vector_all
id.vector <- id.vector.all # id.vector.known or id.vector.all

#-----------------------Make colors ----------------------------------#
GK <- as.character(class_vector)
col.idxG <- rep(NA, length(GK))
ugk <- unique(GK[!is.na(GK)])
ugk <- sort(ugk)
for (i in 1:length(ugk)) {
  col.idxG[GK %in% ugk[i]] <- i
}
col.idxG[is.na(col.idxG)] <- length(ugk) + 1
my.colsG <-c("red", "orange", "green3", "palegreen", "purple", "blue", "pink", "powderblue", "yellow", gray(0)) #Need as many colors as classes plus gray at the end
#my.colsG <-c("red", "green3", "purple", "pink", "powderblue", "yellow", gray(0))

#----------------------Distance and clustering ------------------------#
d <- as.dist(1 - cor(data, use='pairwise.complete.obs'))
hca <- hclust(d, method='ward.D2')
use.data.t <- t(data)
d2 <- as.dist(1 - cor(use.data.t, use='pairwise.complete.obs'))
hca2 <- hclust(d2, method='ward.D2')

##OBS here i am only clustering on my chosen CpG sites  
##But you can do the clustering across all the DNA methylation levels on the array here 
##And just plot the meth values for your top 100, 1000 etc.


#-----------------------MAKE HEATMAP ----------------------------------#
pdf('AML_heatmap.pdf', paper='special', width=6, height=3,5)
par(mar=c(0,0,0,0), omi=c(0.3,0,0,0), ps=8)

# To label samples, enter a vector  for labCol with the names in the same order as in data, use id.vector

heatmap(data, Colv=as.dendrogram(hca, hang=-1), Rowv=as.dendrogram(hca2, hang=-1), 
        ColSideColors=my.colsG[col.idxG], labCol= id.vector, labRow= '', 
        scale="none", col= rev(brewer.pal(10, 'RdYlBu')),
        margins=c(1,1), cexCol=1, cexRow = 1, useRaster=T)

lc1 <- legend('topleft', fill=my.colsG[-length(my.colsG)], legend=ugk, 
              bty='n', ncol=2)
dev.off()

# Add the subtype to the pdf file name!!!
#--------------------------------------

##### Making csv file for analysis of sample clustering #####

# See above: aml.pheno_heatmap <- aml.pheno_clean_nonrelapse_unknown[samples.unknown.heatmap, ]
write.csv(aml.pheno_heatmap, file = "aml_pheno_heatmap.csv")
# See above: id.vector.all <- c(phenoData_clean_nonrelapse_known$Sample_Name, phenoData_clean_nonrelapse_unknown$Sample_Name[samples.unknown.heatmap])
write.csv(id.vector[hca$order], file = "heatmap_all_id.csv")
# See above: aml.pheno_heatmap.all <- rbind(aml.pheno_clean_nonrelapse_known, aml.pheno_clean_nonrelapse_unknown[samples.unknown.heatmap, ])
write.csv(aml.pheno_heatmap.all[hca$order, ], file = "aml_pheno_heatmap_ordered.csv")

# After making heatmap only on known samples:

write.csv(id.vector[hca$order], file = "heatmap_id_known.csv")
write.csv(aml.pheno_heatmap.all[hca$order, ], file = "aml_pheno_heatmap_ordered_known.csv")

# I put these files in the same sheet, see file "aml_pheno_heatmap.csv"
# The file has been rearranged to facilitate analysis
