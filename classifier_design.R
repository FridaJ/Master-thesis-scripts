##### Classifier design, corrected version using test sets properly #####

##### About the script #####

# This script will take an AML subtype as input and return a dataframe df.results.subtype with 
# the results from a 5x5-fold cross validated NSC job, to create a classifier for the specified 
# subtype. It will also return a list of cpg sites with their frequency of occurrence in the 25 
# runs. Specify the AML subtype below under "Preparing to run the script".

# The results from the script will be stored in two variables: df.results.subtype and cpg.freq.subtype
# df.results will have 25 rows (one for each run) and 6 columns:
# "threshold" is the threshold used for NSC shrinking of centroids, which gave the least amount of errors when predicting class of the test set
# "N_CpGs" is the number of CpG sites used to classify the subtype using the chosen threshold
# "error_predict" is the number of errors that arise when testing the classifier on the test set
# "CpG_IDs" is a list of indices for the chosen CpG sites
# "test_set" is a list of indices for the samples that are included in the test set for the run
#
# cpg.freq.subtype is an ordered dataframe with the indices of all chosen CpG sites for all 25 runs
# and the frequency with which they were chosen. 25 would been chosen in all runs, 8 would mean
# chosen in 8 runs out of 25.

##### Preparing to run the script #####

library(pamr)
library(emil)
setwd("Documents/R_work")

options(stringsAsFactors = FALSE) #removes the meta level in dataframes?
load("AML_subtype_lists.Rdata") # loads lists that are input for pamr.train() one for each subtype

subtype <- "MLL" # Specify subtype here, "t821", "inv16", "mono7", or "MLL"
subtype.list <- MLL.list # Specify corresponding subtype list to be used

##### Initalizing result dataframe #####

df.results <- data.frame(matrix(ncol = 3, nrow = 25))
colnames(df.results) = c("threshold", "N_CpGs", "error_predict")
df.results$CpG_IDs <- list(NA) # adds a new column "CpG_IDs" with class "list" 
df.results$test_set <- list(NA) # adds a new column "test_set" with class "list"

##### Subsetting for 5x5-fold cross validation #####

# Function for storing the indices for the test set in df.results

store.testset <- function(df.results, test.indices, row.number) {
  if (!missing(test.indices)) {
    df.results$test_set[[row.number]] <- test.indices # stores a numerical vector in row row.number
  }
  return(df.results)
}

# Do the 5x5-fold subsetting into test sets and training sets, specify subtype in beginning

class.vector <- as.factor(subtype.list$y)
resample_result <- resample_crossvalidation(y = class.vector, subset = c(1:length(class.vector))) #from emil package

# To store the test set indices in the result dataframe:
for (i in c(1:ncol(resample_result$fold_set))) {
  test.set <- which(resample_result$fold_set[, i] == F)
  df.results <- store.testset(df.results, test.set, i)
}

##### Functions needed for getting and storing NSC data #####

get.thresholds.ncpgs <- function(cv.results) {
  
  # 30 thresholds will be analyzed with respect to number of CpG sites after shrinking
  # Those thresholds that give at least 5 and less than 100 sites will be stored in a list
  # The corresponding number of CpG sites will be stored in another list
  thresholds <- list()
  ncpgs <- list()
  row <- length(cv.results$size) # to get the starting index for the loop below
  for (i in c(row:1)) {
    if ((i == length(cv.results$size)) & (cv.results$size[i] >= 100)) {
      print("Too many CpG sites chosen, please use other data.")
    } else if ((cv.results$size[i] >= 5) & (cv.results$size[i] < 100)) {
      thresholds[i] <- cv.results$threshold[i]
      ncpgs[i] <- cv.results$size[i]
    }
  }
  # thresholds is now a list containing the thresholds to be used for predicting class of the test sets
  # ncpgs is a list containing the corresponding number of CpG sites, indexed in the same way
  thresholds <- unlist(thresholds) # removes NULL elements
  ncpgs <- unlist(ncpgs)
  return(list(thresholds = thresholds, ncpgs = ncpgs))
}

store.cpg_ids <- function(df.results, cpg_ids, row.number) {
  if (!missing(cpg_ids)) {
    df.results$CpG_IDs[[row.number]] <- cpg_ids # stores a vector of indices in the CpG_IDs 
  }
  return(df.results)
}

store.numeric.results <- function(df.results, threshold, N_CpGs, error.predict = NULL, row.number){
  df.results$threshold[row.number] <- threshold
  df.results$N_CpGs[row.number] <- N_CpGs
  df.results$error_predict[row.number] <- error.predict
  return(df.results)
}

make.training.list <- function(df.results, data.list, row.number) {
  training.data <- data.list$x[, -df.results$test_set[[row.number]]]
  training.classes <- data.list$y[-df.results$test_set[[row.number]]]
  training.list <- list(x = training.data, y = training.classes)
  return(training.list)
}  

make.test.list <- function(df.results, data.list, row.number) {
  test.data <- data.list$x[, df.results$test_set[[row.number]]]
  test.classes <- data.list$y[df.results$test_set[[row.number]]]
  test.list <- list(x = test.data, y = test.classes)
  return(test.list)
} 

# Function for pamr usage and data storage for 1 run of NSC

train.and.evaluate <- function(df.results, training.list, test.list){
  
  trained.data <- pamr.train(training.list)
  cv.results <- pamr.cv(trained.data, training.list)
  # Here, get all thresholds from cv.results that have at least 5 cpg sites but less than 100, store in a list
  thresholds.ncpgs.list <- get.thresholds.ncpgs(cv.results)
  thresholds.list <- thresholds.ncpgs.list$thresholds
  ncpgs.list <- thresholds.ncpgs.list$ncpgs
  # Then, do pamr.predict() and calculate the predicted error for each of these thresholds using the test set, store in a list
  error.predict.list <- list()
  for (i in c(1:length(thresholds.list))) {
    predicted <- pamr.predict(trained.data, test.list$x, threshold = thresholds.list[i]) # predict the class of the samples in test.list
    summary.predicted <- as.matrix(summary(predicted == test.list$y)) # have to be same number of LEVELS!!
    # compare classes in predict.X vector to classes in class vector prepared previously to evaluate
    # store as matrix to make data accessible
    if (!("FALSE" %in% dimnames(summary.predicted)[[1]])) {
      incorrect <- as.integer(summary.predicted[3]) # number of incorrect predictions (only NAs)
      all <- as.integer(summary.predicted[2])+as.integer(summary.predicted[3]) # total number of predictions
      error.predict <- incorrect/all
    } else if ("FALSE" %in% dimnames(summary.predicted)[[1]]) {
      incorrect <- as.integer(summary.predicted[2]) + as.integer(summary.predicted[4]) # number of incorrect predictions including NAs
      all <- as.integer(summary.predicted[2])+as.integer(summary.predicted[3]) + as.integer(summary.predicted[4]) # total number of predictions
      error.predict <- incorrect/all
    }
    error.predict.list[i] <- error.predict
  }
  # error.predict.list now contains the prediction errors corresponding to the thresholds in thresholds.list
  # Get the threshold for the lowest error in that list (with lowest number of CpG sites (= highest threshold))
  # Reverse the lists to be able to use which.min() and get the correct index (if multiple indices with same error)
  error.predict.rev <- unlist(rev(error.predict.list))
  thresholds.rev <- rev(thresholds.list)
  ncpgs.rev <- rev(ncpgs.list)
  index.lowest.error <- which.min(error.predict.rev)
  # Use this index to get input for the results dataframe (df.results)
  error.predict <- error.predict.rev[index.lowest.error]
  threshold <- thresholds.rev[index.lowest.error]
  N_CpGs <- ncpgs.rev[index.lowest.error]
  # to get the indices of the cpg sites for threshold = threshold:
  cpg_ids <- pamr.predict(trained.data, training.list$x, threshold = threshold, type = "nonzero") # a list of indices for the CpG sites used
  df.results <- store.cpg_ids(df.results, cpg_ids, row.number) # stores the indices of the non-shrunken cpg sites
  df.results <- store.numeric.results(df.results, threshold, N_CpGs, error.predict, row.number)
  
  return(df.results)
}

##### Script that calls previous functions to do NSC and stores results #####

for (i in c(1:ncol(resample_result$fold_set))) {
#for (i in c(1:2)) { #only for testing
  row.number <- i
  print(paste0("Run ", row.number))
  test.list <- make.test.list(df.results, data.list = subtype.list, row.number = i)
  training.list <- make.training.list(df.results, data.list = subtype.list, row.number = i)
  df.results <- train.and.evaluate(df.results, training.list = training.list, test.list = test.list)
}

assign(paste0("df.results", ".", subtype), df.results)

##### Count the frequency of cpg ids in a list #####

get.cpg.freq <- function(lst) {
  freq.list <- list()
  for (i in c(1:length(lst))) {
    if (!lst[i] %in% freq.list$index) { # do only if the element is not already in freq.list
      freq.list$index <- c(freq.list$index, lst[i]) # appends cpg site index to freq.list
      freq.list$freq <- c(freq.list$freq, sum(lst == lst[i])) # appends number of cpg site occurances 
    }
  }
  freq.list <- as.data.frame(freq.list) # make into dataframe to be able to use order()
  freq.list.ordered <- freq.list[order(-freq.list$freq), ] # order list according to frequency
  return(freq.list.ordered)
}

start.list <- unlist(df.results$CpG_IDs) # make one long list of all cpg sites found in all runs
cpg.freq <- get.cpg.freq(start.list) # a list of cpg sites and the frequency they were chosen
assign(paste0("cpg.freq", ".", subtype, ".new"), cpg.freq)

##### Removing temp variables #####

rm(df.results)
rm(start.list)
rm(cpg.freq)
rm(row.number)
rm(test.list)
rm(training.list)
rm(resample_result)
rm(i)
rm(test.set)

##### Saving objects #####

# After getting cpg.freq for all subtypes, save them for later use:
save(cpg.freq.t821.new, cpg.freq.inv16.new, cpg.freq.mono7.new, cpg.freq.MLL.new, file = "AML_cpg.freqs.new.Rdata")

# Save the df.results objects
save(df.results.t821, df.results.inv16, df.results.mono7, df.results.MLL, file = "AML_df.results.Rdata")
