setwd("Documents/R_work")
library(pamr)

load("GRSet_beta_datasets.Rdata")
load("AML_cpg.freqs.new.Rdata")
load("final_classifiers.Rdata")

##### Predict the samples in the known dataset #####
# Save this result for a "confusion table" figure, with all four classes
# Load phenotype data to check the predictions of known samples

load(file="AML_pheno_info.Rdata")

GRSet_beta_known.t821cpgs <- GRSet_beta_known[cpg.freq.t821.new$index, ] # Use all cpg sites found for this classifier
predicted <- pamr.predict(classifier.data.t821, GRSet_beta_known.t821cpgs, threshold = 0) #threshold = 0 to use all cpg sites
summary(predicted)

#not t(8;21)     t(8;21) 
#57          20 

t821.predicted.known <- which(predicted %in% "t(8;21)") # gives a vector of indexes (unknown data) that are predicted t821
aml.pheno_clean_nonrelapse_known$genotype[t821.predicted.known] # only t821 samples

GRSet_beta_known.inv16cpgs <- GRSet_beta_known[cpg.freq.inv16.new$index[1:6], ] # Use 6 cpg sites for this classifier
predicted <- pamr.predict(classifier.data.inv16, GRSet_beta_known.inv16cpgs, threshold = 0) #threshold = 0 to use all cpg sites
summary(predicted)

#inv(16) not inv(16) 
#13          64 

inv16.predicted.known <- which(predicted %in% "inv(16)") # gives a vector of indexes (unknown data) that are predicted inv16
aml.pheno_clean_nonrelapse_known$genotype[inv16.predicted.known] # only inv16 samples

GRSet_beta_known.mono7cpgs <- GRSet_beta_known[cpg.freq.mono7.new$index, ] # Use all cpg sites found for this classifier
predicted <- pamr.predict(classifier.data.mono7, GRSet_beta_known.mono7cpgs, threshold = 0) #threshold = 0 to use all cpg sites
summary(predicted)

#mono 7 not mono 7 
#5         72 

mono7.predicted.known <- which(predicted %in% "mono 7") # gives a vector of indexes (unknown data) that are predicted mono7
aml.pheno_clean_nonrelapse_known$genotype[mono7.predicted.known] # only mono7 samples

GRSet_beta_known.MLLcpgs <- GRSet_beta_known[cpg.freq.MLL.new$index[1:13], ] # Use 13 cpg sites for this classifier
predicted <- pamr.predict(classifier.data.MLL, GRSet_beta_known.MLLcpgs, threshold = 0) #threshold = 0 to use all cpg sites
summary(predicted)

#MLL not MLL 
#30      47 

MLL.predicted.known <- which(predicted %in% "MLL") # gives a vector of indexes (unknown data) that are predicted MLL
MLL.predicted.genotype <- aml.pheno_clean_nonrelapse_known$genotype[MLL.predicted.known] # one sample with genotype "sole+8"!!!
# The sole+8 sample means that one MLL sample is also missed, since there are 30 MLL samples in total

# Check the id of the sole+8 patient:
which(MLL.predicted.genotype == "sole+8")
#[1] 11
MLL.predicted.known[11]
#[1] 16
aml.pheno_clean_nonrelapse_known$id[16]
#[1] "03/324"
aml.pheno_newids_known$id[16]
#[1] "03/324"
aml.pheno_newids_known$newid[16]
#[1] "AML016"

# Also check the patient AML022 (01/053) and see if this is the one not predicted as MLL (seen in the heatmaps)
which(aml.pheno_clean_nonrelapse_known$id == "01/053")
#[1] 22
aml.pheno_newids_known$newid[22]
#[1] "AML022"
aml.pheno_newids_known$id[22]
#[1] "01/053"
MLL.predicted.known
#[1]  1  4  5  6  7  8  9 10 11 12 16 23 28 33 34 35 38 39 40 49 51 52 53 59 61 63 69 70 71 74
# Index 22 is not included, also:
aml.pheno_newids_known$newid[MLL.predicted.known]
#[1] "AML001" "AML004" "AML005" "AML006" "AML007" "AML008" "AML009" "AML010" "AML011" "AML012" "AML016" "AML023"
#[13] "AML028" "AML033" "AML034" "AML035" "AML038" "AML039" "AML040" "AML049" "AML051" "AML052" "AML053" "AML059"
#[25] "AML061" "AML063" "AML069" "AML070" "AML071" "AML074"

##### Using the final classifiers on unknown data #####

# Subset the unknown data to the relevant CpG sites, then use pamr.predict
# Load final_classifiers.Rdata if needed

load(file = "final_classifiers.Rdata")

GRSet_beta_unknown.t821cpgs <- GRSet_beta_unknown[cpg.freq.t821.new$index, ] # Use all cpg sites found for this classifier
predicted <- pamr.predict(classifier.data.t821, GRSet_beta_unknown.t821cpgs, threshold = 0) #threshold = 0 to use all cpg sites
summary(predicted) 

# not t(8;21)     t(8;21) 
# 58           0 

t821.predicted <- which(predicted %in% "t(8;21)") # gives a vector of indexes (unknown data) that are predicted t821

GRSet_beta_unknown.inv16cpgs <- GRSet_beta_unknown[cpg.freq.inv16.new$index[1:6], ] # Use 6 cpg sites for this classifier
predicted <- pamr.predict(classifier.data.inv16, GRSet_beta_unknown.inv16cpgs, threshold = 0) #threshold = 0 to use all cpg sites
summary(predicted) 

# inv(16) not inv(16) 
# 3          55 

inv16.predicted <- which(predicted %in% "inv(16)") # gives a vector of indexes (unknown data) that are predicted inv16

GRSet_beta_unknown.mono7cpgs <- GRSet_beta_unknown[cpg.freq.mono7.new$index, ] # Use all cpg sites found
predicted <- pamr.predict(classifier.data.mono7, GRSet_beta_unknown.mono7cpgs, threshold = 0) #threshold = 0 to use all cpg sites
summary(predicted) 

# mono 7 not mono 7 
# 5         53 

mono7.predicted <- which(predicted %in% "mono 7") # gives a vector of indexes (unknown data) that are predicted mono7

GRSet_beta_unknown.MLLcpgs <- GRSet_beta_unknown[cpg.freq.MLL.new$index[1:13], ] # Use 13 most common cpg sites
predicted <- pamr.predict(classifier.data.MLL, GRSet_beta_unknown.MLLcpgs, threshold = 0) #threshold = 0 to use all cpg sites
summary(predicted) 

# MLL not MLL 
# 7      51 

MLL.predicted <- which(predicted %in% "MLL") # gives a vector of indexes (unknown data) that are predicted MLL

save(t821.predicted, inv16.predicted, mono7.predicted, MLL.predicted, file = "subtypes_predicted.Rdata")

