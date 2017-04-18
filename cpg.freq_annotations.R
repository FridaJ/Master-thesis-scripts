##### Preparing for analysis of CpG sites in classifiers, saving relevant data in dataframes #####

load("AML_data_for_pamr.Rdata") # to get annotation.dataframe_for_pamr

# use cbind to bind cpg.freq to the subsetted annotation file with the CpG site information
# the indexes of cpg.freq correspond to row numbers in the annotation dataframe
# cbind the following colnames: chr, pos, Name, Relation_to_Island, UCSC_RefGene_Name, UCSC_RefGene_Accession, UCSC_RefGene_Group
# do this for all the diagnostic subtypes

cpg.freq.t821_annotation_all <- cbind(cpg.freq.t821.new, annotation.dataframe_for_pamr[cpg.freq.t821.new$index, c("Name", "chr", "pos", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
cpg.freq.inv16_annotation_all <- cbind(cpg.freq.inv16.new, annotation.dataframe_for_pamr[cpg.freq.inv16.new$index, c("Name", "chr", "pos", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
cpg.freq.mono7_annotation_all <- cbind(cpg.freq.mono7.new, annotation.dataframe_for_pamr[cpg.freq.mono7.new$index, c("Name", "chr", "pos", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
cpg.freq.MLL_annotation_all <- cbind(cpg.freq.MLL.new, annotation.dataframe_for_pamr[cpg.freq.MLL.new$index, c("Name", "chr", "pos", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])

# Make dataframes with only the information about the cpg sites chosen for the final classifiers:

cpg.freq.t821_annotation_final <- cpg.freq.t821_annotation_all
cpg.freq.inv16_annotation_final <- cpg.freq.inv16_annotation_all[1:6,]
cpg.freq.mono7_annotation_final <- cpg.freq.mono7_annotation_all
cpg.freq.MLL_annotation_final <- cpg.freq.MLL_annotation_all[1:13,]

save(cpg.freq.t821_annotation_final, cpg.freq.inv16_annotation_final, cpg.freq.mono7_annotation_final, cpg.freq.MLL_annotation_final, file = "subtype_annotations_final.Rdata")
