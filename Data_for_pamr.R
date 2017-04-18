library(pamr)
setwd("Documents/R_work")

###### Load files into R environment #####

options(stringsAsFactors = FALSE) #removes the meta level in dataframes?
load("AML_data_for_pamr.Rdata") # loads GRSet_betavalues, annotation.dataframe_for_pamr, aml.pheno_clean, aml.info_clean and phenoData_clean

###### Remove relapse samples from data ######

relapse_indices <- which(aml.pheno_clean$sample.type %in% "relapse") # vector of indices for relapse samples
# length of vector = 25
GRSet_betavalues_nonrelapse <- GRSet_betavalues[, -relapse_indices] # removes the relapse samples from the data
aml.pheno_clean_nonrelapse <- aml.pheno_clean[-relapse_indices, ]
aml.info_clean_nonrelapse <- aml.info_clean[-relapse_indices, ]
phenoData_clean_nonrelapse <- phenoData_clean[-relapse_indices, ]

save(GRSet_betavalues_nonrelapse, aml.pheno_clean_nonrelapse, aml.info_clean_nonrelapse, phenoData_clean_nonrelapse, file = "AML_clean_nonrelapse.Rdata")

###### Divide AML data into known and unknown subtypes ######

unknown_indices <- which(aml.pheno_clean_nonrelapse$genotype %in% c("normal", "other clon abn", "no result")) 
# gives a vector with indexes for unknown samples
# length = 58
GRSet_beta_unknown <- GRSet_betavalues_nonrelapse[, unknown_indices] # only keeps samples with index included in vector 
# dim() gives 58 columns
GRSet_beta_known <- GRSet_betavalues_nonrelapse[, -unknown_indices] # only keeps samples with index NOT in vector
# dim() gives 77 columns

subtype_vector_known <- aml.pheno_clean_nonrelapse$genotype[-unknown_indices] # vector with subtypes for the samples with known subtype
# length = 77
subtype_vector_unknown <- aml.pheno_clean_nonrelapse$genotype[unknown_indices] # vector with subtypes for the samples with unknown subtype
# length = 58

aml.pheno_clean_nonrelapse_known <- aml.pheno_clean_nonrelapse[-unknown_indices, ] # Phenotype data for known samples
aml.pheno_clean_nonrelapse_unknown <- aml.pheno_clean_nonrelapse[unknown_indices, ] # Phenotype data for unknown samples
aml.info_clean_nonrelapse_known <- aml.info_clean_nonrelapse[-unknown_indices, ] # Info for known samples
aml.info_clean_nonrelapse_unknown <- aml.info_clean_nonrelapse[unknown_indices, ] # Info for unknown samples
phenoData_clean_nonrelapse_known <- phenoData_clean_nonrelapse[-unknown_indices, ]
phenoData_clean_nonrelapse_unknown <- phenoData_clean_nonrelapse[unknown_indices, ]

save(aml.pheno_clean_nonrelapse_known, aml.pheno_clean_nonrelapse_unknown, aml.info_clean_nonrelapse_known, aml.info_clean_nonrelapse_unknown, phenoData_clean_nonrelapse_known,  phenoData_clean_nonrelapse_unknown, file = "AML_pheno_info.Rdata")

complete <- length(complete.cases(GRSet_beta_known)) == length(rownames(GRSet_beta_known)) # TRUE if all probes are complete
cat("All samples complete? ", complete)
# If all are not complete, remove the incomplete

save(GRSet_beta_unknown, GRSet_beta_known, file = "GRSet_beta_datasets.Rdata")

###### Change the class names of the NOT X classes ######

# In subtype_vector_known, change all !X to "not X"

t821_class_vector <- as.character(subtype_vector_known) # make character instead of factor, for replacement
t821_class_vector[!t821_class_vector %in% "t(8;21)"] <- "not t(8;21)" # replaces the other classes with "not t(8;21)"

MLL_class_vector <- as.character(subtype_vector_known)
MLL_class_vector[MLL_class_vector %in% c("other 11q23/MLL", "t(9;11)", "t(10;11)", "t(11;19)")] <- "MLL"
MLL_class_vector[!MLL_class_vector %in% "MLL"] <- "not MLL"

inv16_class_vector <- as.character(subtype_vector_known)
inv16_class_vector[!inv16_class_vector %in% "inv(16)"] <- "not inv(16)"

mono7_class_vector <- as.character(subtype_vector_known)
mono7_class_vector[!mono7_class_vector %in% "mono 7"] <- "not mono 7"

###### Making lists for pamr.tain() input ######

t821.list <- list(x = GRSet_beta_known, y = t821_class_vector) 
MLL.list <- list(x = GRSet_beta_known, y = MLL_class_vector) 
inv16.list <- list(x = GRSet_beta_known, y = inv16_class_vector) 
mono7.list <- list(x = GRSet_beta_known, y = mono7_class_vector) 

save(t821.list, MLL.list, inv16.list, mono7.list, file = "AML_subtype_lists.Rdata")
save(t821_class_vector, MLL_class_vector, inv16_class_vector, mono7_class_vector, file = "AML_class_vectors.Rdata")

##### For making heatmaps with only the 4 classes with classifiers #####

# This will remove other classes from class vectors and beta value data objects

indices_4classes <- which(subtype_vector_known %in% c("t(8;21)", "inv(16)", "mono 7", "other 11q23/MLL", "t(9;11)", "t(10;11)", "t(11;19)"))
subtype_vector_known_4classes <- subtype_vector_known[indices_4classes]

GRSet_beta_known_4classes <- GRSet_beta_known[,indices_4classes] # Use this in AML_heatmaps.R in place of GRSet_beta_known

# Now subset the class vectors, removing the samples that are not of the four relevant subtypes:

t821_class_vector_4classes <- as.character(subtype_vector_known_4classes) # make character instead of factor, for replacement
t821_class_vector_4classes[!t821_class_vector_4classes %in% "t(8;21)"] <- "not t(8;21)" # replaces the other classes with "not t(8;21)"

MLL_class_vector_4classes <- as.character(subtype_vector_known_4classes)
MLL_class_vector_4classes[MLL_class_vector_4classes %in% c("other 11q23/MLL", "t(9;11)", "t(10;11)", "t(11;19)")] <- "MLL"
MLL_class_vector_4classes[!MLL_class_vector_4classes %in% "MLL"] <- "not MLL"

inv16_class_vector_4classes <- as.character(subtype_vector_known_4classes)
inv16_class_vector_4classes[!inv16_class_vector_4classes %in% "inv(16)"] <- "not inv(16)"

mono7_class_vector_4classes <- as.character(subtype_vector_known_4classes)
mono7_class_vector_4classes[!mono7_class_vector_4classes %in% "mono 7"] <- "not mono 7"

# save the beta values and class vectors in an R object to load in the script AML_heatmaps.R

save(indices_4classes, GRSet_beta_known_4classes, t821_class_vector_4classes, MLL_class_vector_4classes, inv16_class_vector_4classes, mono7_class_vector_4classes, file = "for_heatmaps_4classes.Rdata")
