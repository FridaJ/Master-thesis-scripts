library(minfi)
library(IlluminaHumanMethylation450kmanifest)

# If working locally
setwd("Documents/R_work")

###### Load files into R environment #####

options(stringsAsFactors = FALSE) #removes the meta level in dataframes?

load("AMLdata_clean_samples.Rdata") # loads GRSet.funnorm_clean, RGSet_clean, aml.pheno_clean, aml.info_clean, phenoData_clean

####### finding CpG sites with high p values, for CLEAN data #######

detP <- detectionP(RGSet_clean)
failed <- detP > 0.0001
colMeans(failed) # Fraction of failed probes per sample
sum(rowMeans(failed)>0) # How many probes had at least one failed sample?

probes_clean <- which(rowMeans(failed) == 0) #vector of indices with probes to use

# Get info about probes, e.g. chr, type, etc

annotation.dataframe <- getAnnotation(GRSet.funnorm_clean) # make dataframe from annotation data in genomic ratio set

####### removing CpG sites with any failed probes (p > 0.0001) ########

GRSet.funnorm_clean_probes <- GRSet.funnorm_clean[probes_clean, ] # remove bad probes from GRSet
RGSet_clean_probes <- RGSet_clean[probes_clean, ]
annotation.dataframe_clean_probes <- annotation.dataframe[probes_clean, ]

####### removing irrelevant and bad probes as determined by ALL study ########

SupplementalTable1 <- read.csv(file = "SupplementalTable1.txt", sep = "\t", header = T) # Table from ALL paper

# Make vector of CpG sites (col TargetID) with the column "analyzed" == 1

probes_to_keep <- SupplementalTable1[(SupplementalTable1$analyzed == 1), ]
probes_id <- probes_to_keep$TargetID # vector containing the probe ids of the probes to keep

# Go through the CpG sites in AML data files and only keep the CpG sites from the vector above
# This will check if a rowname in the data exists in the vector probes_id and give T or F. The rows with T values will be kept

GRSet_for_pamr <- GRSet.funnorm_clean_probes[rownames(GRSet.funnorm_clean_probes) %in% probes_id, ]
annotation.dataframe_for_pamr <- annotation.dataframe_clean_probes[rownames(annotation.dataframe_clean_probes) %in% probes_id, ]

# get beta values for classifier design

GRSet_betavalues <- getBeta(GRSet_for_pamr)

####### Saving files #######

# save the beta values of GRSet_for_pamr (getBeta)
# also save annotation from GRSet and aml.pheno_clean and aml.info_clean, and phenoData_clean

save(GRSet_betavalues, annotation.dataframe_for_pamr, aml.pheno_clean, aml.info_clean, phenoData_clean, file = "AML_data_for_pamr.Rdata")


