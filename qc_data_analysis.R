library(minfi)
library(IlluminaHumanMethylation450kmanifest)

# If working locally
setwd("Documents/R_work")

# If working on Milou
setwd('/pica/v13/b2015118/private/frida-aml')

###### Load files into R environment #####

options(stringsAsFactors = FALSE) #removes the meta level in dataframes

load(file= "GRSet_funnorm.Rdata") # gets the same name as it was saved as
load(file= "RGSet.Rdata")
load(file= "AML450k_pheno_130201_raw.Rdata") 

###### Quality control and data analysis, can be done on local computer #######

# Check that the order of pheno and RG is the same
phenoData <- pData(RGSet) # contains the patient ids for the methylation array data
tmp.names <- aml.info$assay.id # contains the patient sample ids for the phenotype files (aml.info)
tmp.names <- gsub(".AVG_Beta", "", tmp.names) # Removes the extension from the sample ids
match(phenoData$Sample_Name, tmp.names) # gives a vector with the indexes swapped if needed, check this vector
# Note: it appears that everything was in the same order, so no need to change anything. (Oct 7)
# Otherwise:
#df_sorted <- df[match(phenoData$Sample_Name, tmp.names),] #would probably give a sorted file, not tested

# QCplot to check quality of samples
MSet <- preprocessRaw(RGSet)
qc <- getQC(MSet)
plotQC(qc) #index 13 outlier

# Control that the bisulfite conversion has worked
controlStripPlot(RGSet, controls="BISULFITE CONVERSION II")
#To find the index/row number of a specific ID:
which(rownames(phenoData_GR) %in% "7668610011_R01C01") # index 13 outlier

#Density plots to compare sample groups, example
densityPlot(MSet, sampGroups = aml.pheno$sample.type)
#without legend, set legend = FALSE

#Plot predicted sex as X/Y intensities
predictedSex <- getSex(GRSet.funnorm, cutoff = -2)
plotSex(predictedSex)
#If you want to color the male/female groups according to sex given in dataset, exchange column predictedSex
#Instead use sex_dataset:
sex_dataset <- aml.pheno$sex
sex_dataset <- as.character(sex_dataset)
sex_dataset[sex_dataset %in% "female"] <- "F"
sex_dataset[sex_dataset %in% "male"] <- "M"
#Make sex_dataset into a dataframe so cbind can be used
sex_dataset_df <- as.data.frame(sex_dataset)
#Combine the intensities in predictedSex with the sex from dataset
datasetSex <- cbind(predictedSex[,1:2], sex_dataset_df)
colnames(datasetSex)[3] <- "predictedSex" #Colname has to be predictedSex
#Now plot the sex distribution
plotSex(datasetSex)

#Compare predicted sex with sex in dataset

phenoData_GR <- pData(GRSet.funnorm)
sex_predicted <- phenoData_GR$predictedSex
sex_dataset <- aml.pheno$sex
sex_dataset <- as.character(sex_dataset)
sex_dataset[sex_dataset %in% "female"] <- "F"
sex_dataset[sex_dataset %in% "male"] <- "M"

# This gives a summary over mismatches
summary(sex_dataset == sex_predicted)

# This gives a list of which samples were incorrect
sex.matched <- sex_dataset == sex_predicted
sex.matched[is.na(sex.matched)] <- "Unknown" #Need to change NA to "unknown" or NA samples will not be included
idx <- which(sex.matched %in% c("Unknown", "FALSE")) #gives indexes for the non-matching samples

#This shows the pheno data for these samples. Only 02/235 is of concern, others are pos and neg controls so it was expected.


############### Principal Compoment Analysis, PCA ###############

#DO THIS ON UPPMAX!!!

# Setting up for doing a PCA
library(ggplot2)

beta <- getBeta(GRSet.funnorm)
tmp <- is.na(beta[,1])
summary(tmp) # This tells us if we have NAs in our dataset, which we dont

#d <- cor(beta) # shows the correlation between all samples in a matrix
fit <- prcomp(t(beta)) #takes a while
#retx= T, center = T, scale = T) 
attributes(fit)
summary(fit)

phenoData <- pData(GRSet.funnorm)
pca.data <- as.data.frame(phenoData)
pca.data <- cbind(pca.data, fit$x)

#This plot shows the relation between the samples and which position they were run on
# No obvious positional effect
pdf("/proj/b2015118/private/nobackup/frida-aml-nobackup/plots/Position_pca.pdf", width=7, height=7, useDingbats = FALSE)
ggplot(pca.data, aes(x=PC1, y=PC2, color=Array))+
  geom_point()
dev.off()

#This plot shows the relation between the samples and which slide/array they were run on
# No obvious batch effect
pdf("/proj/b2015118/private/nobackup/frida-aml-nobackup/plots/Array_pca.pdf", width=7, height=7, useDingbats = FALSE)
ggplot(pca.data, aes(x=PC1, y=PC2, color=Slide))+
  geom_point()
dev.off()


#Now lets look at the pos and neg controls vs AML samples
pca.data <- cbind(aml.pheno, fit$x)

#This plot shows the relation between the sample types
# No obvious batch effect
pdf("/proj/b2015118/private/nobackup/frida-aml-nobackup/plots/SampleType_pca.pdf", width=7, height=7, useDingbats = FALSE)
ggplot(pca.data, aes(x=PC1, y=PC2, color=sample.type))+
  geom_point()
dev.off()

# If you want to label the data points with index numbers (row numbers), use this
pdf("/proj/b2015118/private/nobackup/frida-aml-nobackup/plots/SampleType_labeled_pca.pdf", width=7, height=7, useDingbats = FALSE)
ggplot(pca.data, aes(x=PC1, y=PC2, color=sample.type, label=1:168))+
  geom_point()+
  geom_text()
dev.off()

##### All of the above PCA but with "clean" data, outliers and controls removed

#First, cleanup of data
GRSet.funnorm_clean <- GRSet.funnorm[, c(1:7,14:35,37:131,133:168)]
RGSet_clean <- RGSet[, c(1:7,14:35,37:131,133:168)]
aml.pheno_clean <- aml.pheno[c(1:7,14:35,37:131,133:168), ]
aml.info_clean <- aml.info[c(1:7,14:35,37:131,133:168), ]
phenoData_clean <- phenoData[c(1:7,14:35,37:131,133:168), ]

save(GRSet.funnorm_clean, RGSet_clean, aml.pheno_clean, aml.info_clean, phenoData_clean, file = "AMLdata_clean_samples.Rdata")

#PCA
beta_clean <- getBeta(GRSet.funnorm_clean)
fit_clean <- prcomp(t(beta_clean))
phenoData_clean <- pData(GRSet.funnorm_clean)
pca.data_clean <- as.data.frame(phenoData_clean)
pca.data_clean <- cbind(pca.data_clean, fit_clean$x)

pdf("/proj/b2015118/private/nobackup/frida-aml-nobackup/plots/Position_pca_clean.pdf", width=7, height=7, useDingbats = FALSE)
ggplot(pca.data_clean, aes(x=PC1, y=PC2, color=Array))+
  geom_point()
dev.off()

pdf("/proj/b2015118/private/nobackup/frida-aml-nobackup/plots/Array_pca_clean.pdf", width=7, height=7, useDingbats = FALSE)
ggplot(pca.data_clean, aes(x=PC1, y=PC2, color=Slide))+
  geom_point()
dev.off()

pca.data_pheno_clean <- cbind(aml.pheno_clean, fit_clean$x)

pdf("/proj/b2015118/private/nobackup/frida-aml-nobackup/plots/SampleType_pca_clean.pdf", width=7, height=7, useDingbats = FALSE)
ggplot(pca.data_pheno_clean, aes(x=PC1, y=PC2, color=sample.type))+
  geom_point()
dev.off()

pdf("/proj/b2015118/private/nobackup/frida-aml-nobackup/plots/SampleType_labeled_pca_clean.pdf", width=7, height=7, useDingbats = FALSE)
ggplot(pca.data_pheno_clean, aes(x=PC1, y=PC2, color=sample.type, label=1:160))+
  geom_point()+
  geom_text()
dev.off()

#Printing a subset of rows from aml.pheno:
aml.pheno[c(14,24,70,86),] #indexes of clustering outliers



