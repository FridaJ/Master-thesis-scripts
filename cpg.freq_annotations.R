##### Preparing for analysis of CpG sites in classifiers, saving relevant data in dataframes #####

# use cbind to bind cpg.freq to the subsetted annotation file with the CpG site information
# the indexes of cpg.freq correspond to row numbers in the annotation dataframe
# cbind the following colnames: chr, pos, Name, Relation_to_Island, UCSC_RefGene_Name, UCSC_RefGene_Accession, UCSC_RefGene_Group
# do this for all the diagnostic subtypes

cpg.freq.t821_annotation_all <- cbind(cpg.freq.t821, annotation.dataframe_for_pamr[cpg.freq.t821$index, c("Name", "chr", "pos", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
cpg.freq.inv16_annotation_all <- cbind(cpg.freq.inv16, annotation.dataframe_for_pamr[cpg.freq.inv16$index, c("Name", "chr", "pos", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
cpg.freq.mono7_annotation_all <- cbind(cpg.freq.mono7, annotation.dataframe_for_pamr[cpg.freq.mono7$index, c("Name", "chr", "pos", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])
cpg.freq.MLL_annotation_all <- cbind(cpg.freq.MLL, annotation.dataframe_for_pamr[cpg.freq.MLL$index, c("Name", "chr", "pos", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")])

# Make dataframes with only the information about the cpg sites chosen for the final classifiers:

cpg.freq.t821_annotation_final <- cpg.freq.t821_annotation_all
cpg.freq.inv16_annotation_final <- cpg.freq.inv16_annotation_all
cpg.freq.mono7_annotation_final <- cpg.freq.mono7_annotation_all[1:13,]
cpg.freq.MLL_annotation_final <- cpg.freq.MLL_annotation_all[1:13,]

save(cpg.freq.t821_annotation_final, cpg.freq.inv16_annotation_final, cpg.freq.mono7_annotation_final, cpg.freq.MLL_annotation_final, file = "subtype_annotations_final.Rdata")
