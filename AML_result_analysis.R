# load classifier.data.XXX
load(file = "GRSet_beta_datasets.Rdata") # GRSet_beta_unknown, GRSet_beta_known


##### Getting the prediction scores for the classifiers' predictions #####

# In pamr.predict, set the argument type to "posterior" to get the posterior probabilities for each sample

GRSet_beta_unknown.inv16cpgs <- GRSet_beta_unknown[cpg.freq.inv16$index, ] # Use all cpg sites for inv16
GRSet_beta_known.inv16cpgs <- GRSet_beta_known[cpg.freq.inv16$index, ] # Use all cpg sites for inv16
GRSet_beta_all.inv16cpgs <- cbind(GRSet_beta_unknown.inv16cpgs, GRSet_beta_known.inv16cpgs)

predicted.posterior.inv16 <- pamr.predict(classifier.data.inv16, GRSet_beta_known.inv16cpgs, type = "posterior", threshold = 0) #threshold = 0 to use all cpg sites

GRSet_beta_unknown.mono7cpgs <- GRSet_beta_unknown[cpg.freq.mono7$index[1:13], ] # Use 13 most common cpg sites
GRSet_beta_known.mono7cpgs <- GRSet_beta_known[cpg.freq.mono7$index[1:13], ] # Use 13 most common cpg sites
GRSet_beta_all.mono7cpgs <- cbind(GRSet_beta_unknown.mono7cpgs, GRSet_beta_known.mono7cpgs)

predicted.posterior.mono7 <- pamr.predict(classifier.data.mono7, GRSet_beta_known.mono7cpgs, type = "posterior", threshold = 0) #threshold = 0 to use all cpg sites

GRSet_beta_unknown.MLLcpgs <- GRSet_beta_unknown[cpg.freq.MLL$index[1:13], ] # Use 13 most common cpg sites
GRSet_beta_known.MLLcpgs <- GRSet_beta_known[cpg.freq.MLL$index[1:13], ] # Use 13 most common cpg sites
GRSet_beta_all.MLLcpgs <- cbind(GRSet_beta_unknown.MLLcpgs, GRSet_beta_known.MLLcpgs)

predicted.posterior.MLL <- pamr.predict(classifier.data.MLL, GRSet_beta_known.MLLcpgs, type = "posterior", threshold = 0) #threshold = 0 to use all cpg sites

