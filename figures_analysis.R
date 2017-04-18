# Script for making figures for analysis and discussion

load("GRSet_beta_datasets.Rdata") # to get GRSet_beta_known and GRSet_beta_unknown
load("AML_class_vectors.Rdata") # to get class vectors, to get the indices for each subtype
load("AML_clean_nonrelapse.Rdata") # to get GRSet_betavalues_nonrelapse
load("subtypes_predicted.Rdata") # Vectors with indexes of the predicted samples of the unknown data

##### Figure for analysis of MLLT1 gene, boxplot #####

# vector with the intensities for the cpg site affecting the gene
# index for the cpg site is 390288 in annotation file, id is cg00766482

which(rownames(GRSet_beta_known == "cg00766482")) # gives 390288, so use this as index for getting the beta values
cg00766482_betas <- GRSet_beta_known[390288,]  # length is 77, 77 samples

# 4 index vectors for the known samples with classifier subtypes

MLL_indices <- which(MLL_class_vector == "MLL")
t821_indices <- which(t821_class_vector == "t(8;21)")
mono7_indices <- which(mono7_class_vector == "mono 7")
inv16_indices <- which(inv16_class_vector == "inv(16)")

# list for input to boxplot(), 4 cpg site intensity vectors, one for each subtype
# range = 0 makes whiskers extend to extreme values, range = 1 makes them extend to quartiles

subtype_betas_cg00766482.lst <- list()
  
subtype_betas_cg00766482.lst[[1]] <- cg00766482_betas[t821_indices]
subtype_betas_cg00766482.lst[[2]] <- cg00766482_betas[inv16_indices]
subtype_betas_cg00766482.lst[[3]] <- cg00766482_betas[mono7_indices]
subtype_betas_cg00766482.lst[[4]] <- cg00766482_betas[MLL_indices]

boxplot(x = subtype_betas_cg00766482.lst, range = 1, names = c("t(8;21)", "inv(16)", "mono 7", "MLL"))

##### Figure for analysis of MLLT1 affected cpg sites, heatmap #####

# get a list of the affected cpg sites by searching annotation file for MLLT1

MLLT1_cpgs_indices <- grep("MLLT1+", annotation.dataframe_for_pamr$UCSC_RefGene_Name)
MLLT1_data <- GRSet_beta_all[MLLT1_cpgs_indices,] #betavalues for the affected cpg sites, for all samples known + predicted

# make a heatmap for the samples of known subtype plus the predicted samples, for only these cpg sites ("all")

data <- MLLT1_data
class_vector <- class_vector_all_other #class_vector_known, class_vector_all_other, class_vector_known_4classes or class_vector_all_4classes

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
my.colsG <-c("red", "orange", "green3", "palegreen", "purple", "blue", "pink", "powderblue", "yellow", gray(0))#Need as many colors as classes plus gray at the end
#my.colsG <-c("red", "orange", "green3", "palegreen", "purple", "blue", "powderblue", "yellow", gray(0))
#my.colsG <-c("red", "green3", "purple", "pink", "powderblue", "yellow", gray(0))
#my.colsG <-c("red", "green3", "purple", "powderblue", "yellow", gray(0))

#----------------------Distance and clustering ------------------------#
d <- as.dist(1 - cor(data, use='pairwise.complete.obs'))
hca <- hclust(d, method='ward.D2')
#use.data.t <- t(data) #comment this row if not clustering the CpG sites
#d2 <- as.dist(1 - cor(use.data.t, use='pairwise.complete.obs')) #comment this row if not clustering the CpG sites
#hca2 <- hclust(d2, method='ward.D2') #comment this row if not clustering the CpG sites

#-----------------------MAKE HEATMAP ----------------------------------#
pdf('MLLT1_heatmap.pdf', paper='special', width=6, height=3,5)
par(mar=c(0,0,0,0), omi=c(0.3,0,0,0), ps=8)

# To label samples, enter a vector  for labCol with the names in the same order as in data, use id.vector
# Use the first heatmap() input for NOT clustering CpG sites, no dendrogram. Also, this one colors the CpG sites.
# Use the second heatmap() input for ordinary clustering of CpG sites

heatmap(data, Colv=as.dendrogram(hca, hang=-1), Rowv=NA,
        ColSideColors=my.colsG[col.idxG], labRow= '', 
        scale="none", col= rev(brewer.pal(10, 'RdYlBu')),
        margins=c(1,1), cexCol=1, cexRow = 1, useRaster=T)

#heatmap(data, Colv=as.dendrogram(hca, hang=-1), Rowv=as.dendrogram(hca2, hang=-1), 
#        ColSideColors=my.colsG[col.idxG], labCol= id.vector, labRow= '', 
#        scale="none", col= rev(brewer.pal(10, 'RdYlBu')),
#        margins=c(1,1), cexCol=1, cexRow = 1, useRaster=T)

lc1 <- legend('topleft', fill=my.colsG[-length(my.colsG)], legend=ugk, 
              bty='n', ncol=1)
dev.off()


##### Figure for analysis of PRMD8 affected cpg sites, heatmap #####

# get a list of the affected cpg sites by searching annotation file for PRDM8

PRDM8_cpgs_indices <- grep("PRDM8+", annotation.dataframe_for_pamr$UCSC_RefGene_Name)
PRDM8_data <- GRSet_beta_all[PRDM8_cpgs_indices,] #betavalues for the affected cpg sites, for all samples known + predicted

# get a list of the PRDM8 cpg site indices chosen for the mono7 classifier

PRDM8_cpgs_mono7 <- cpg.freq.mono7_annotation_final[grep("PRDM8+", cpg.freq.mono7_annotation_final$UCSC_RefGene_Name),]
PRDM8_cpgs_mono7_indices <- PRDM8_cpgs_mono7$index
PRDM8_data_mono7 <- GRSet_beta_all[PRDM8_cpgs_mono7_indices,]

# make a heatmap for the samples of known subtype plus the predicted samples, for only these cpg sites ("all")

data <- PRDM8_data_mono7
class_vector <- class_vector_all_other #class_vector_known, class_vector_all_other, class_vector_known_4classes or class_vector_all_4classes

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
my.colsG <-c("red", "orange", "green3", "palegreen", "purple", "blue", "pink", "powderblue", "yellow", gray(0))#Need as many colors as classes plus gray at the end
#my.colsG <-c("red", "orange", "green3", "palegreen", "purple", "blue", "powderblue", "yellow", gray(0))
#my.colsG <-c("red", "green3", "purple", "pink", "powderblue", "yellow", gray(0))
#my.colsG <-c("red", "green3", "purple", "powderblue", "yellow", gray(0))

#----------------------Distance and clustering ------------------------#
d <- as.dist(1 - cor(data, use='pairwise.complete.obs'))
hca <- hclust(d, method='ward.D2')
#use.data.t <- t(data) #comment this row if not clustering the CpG sites
#d2 <- as.dist(1 - cor(use.data.t, use='pairwise.complete.obs')) #comment this row if not clustering the CpG sites
#hca2 <- hclust(d2, method='ward.D2') #comment this row if not clustering the CpG sites

#-----------------------MAKE HEATMAP ----------------------------------#
pdf('PRDM8_mono7_heatmap.pdf', paper='special', width=6, height=3,5)
par(mar=c(0,0,0,0), omi=c(0.3,0,0,0), ps=8)

# To label samples, enter a vector  for labCol with the names in the same order as in data, use id.vector
# Use the first heatmap() input for NOT clustering CpG sites, no dendrogram. Also, this one colors the CpG sites.
# Use the second heatmap() input for ordinary clustering of CpG sites

heatmap(data, Colv=as.dendrogram(hca, hang=-1), Rowv=NA,
        ColSideColors=my.colsG[col.idxG], labRow= '', 
        scale="none", col= rev(brewer.pal(10, 'RdYlBu')),
        margins=c(1,1), cexCol=1, cexRow = 1, useRaster=T)

#heatmap(data, Colv=as.dendrogram(hca, hang=-1), Rowv=as.dendrogram(hca2, hang=-1), 
#        ColSideColors=my.colsG[col.idxG], labCol= id.vector, labRow= '', 
#        scale="none", col= rev(brewer.pal(10, 'RdYlBu')),
#        margins=c(1,1), cexCol=1, cexRow = 1, useRaster=T)

lc1 <- legend('topleft', fill=my.colsG[-length(my.colsG)], legend=ugk, 
              bty='n', ncol=1)
dev.off()
