##### About the script #####

# This script will take an AML subtype as input and return a dataframe df.results.subtype with 
# the results from a 5x5-fold cross validated NSC job, to create a classifier for the specified 
# subtype. It will also return a list of cpg sites with their frequency of occurrence in the 25 
# runs. Specify the AML subtype below under "Preparing to run the script".

# The results from the script will be stored in two variables: df.results.subtype and cpg.freq.subtype
# df.results will have 25 rows (one for each run) and 6 columns:
# "threshold" is the threshold used for NSC shrinking of centroids
# "N_CpGs" is the number of CpG sites used to classify the subtype using the chosen threshold
# "error_rate_conf" is the number of errors seen in the confusion table when predicting class
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

subtype <- "inv16" # Specify subtype here, "t821", "inv16", "mono7", or "MLL"
subtype.list <- inv16.list # Specify corresponding subtype list to be used

##### Initalizing result dataframe #####

df.results <- data.frame(matrix(ncol = 4, nrow = 25))
colnames(df.results) = c("threshold", "N_CpGs", "error_rate_conf", "error_predict")
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

get.threshold.ncpgs <- function(cv.results) {
  
  # initializing variables
  lowest_error <- which.min(cv.results$error) # gives index (min row number) of the lowest error
  all_lowest <- which(cv.results$error == cv.results$error[lowest_error]) # gives a list of all indices with this error value
  row.number.min <- all_lowest[which.min(all_lowest)] # gives the lowest row number, at the top of the list
  row.number.max <- all_lowest[which.max(all_lowest)] # gives the start row number, at the bottom of the list with this error rate
 
  # check if size >= 5 exists for this error
    # if false, choose the error value for index row.number.min - 1 and try again
  # check if size < 100 exists for this error
    # if false, choose the error value for index row.number.max + 1 and try again
  # else, find the threshold for lowest size >= 5 
  
  while (!any(cv.results$size[all_lowest] >= 5)) {
    all_lowest <- which(cv.results$error[1:row.number.min-1] == cv.results$error[row.number.min - 1]) # gives a list of all indices with this error value
    row.number.min <- all_lowest[which.min(all_lowest)] # gives the lowest row number, at the top of the list
    row.number.max <- all_lowest[which.max(all_lowest)] # gives the start row number, at the bottom of the list with this error rate
  }
  while (!any(cv.results$size[all_lowest] < 100)) {
    lowest_error <- which.max(cv.results$error == cv.results$error[row.number.max + 1]) # Obs! testing a boolean list, which.max gives the first T
    all_lowest <- which(cv.results$error[(row.number.max+1):length(cv.results$error)] == cv.results$error[lowest_error]) + row.number.max # gives a list of all indices with the new error value
    row.number.min <- all_lowest[which.min(all_lowest)] # gives the lowest row number, at the top of the list
    row.number.max <- all_lowest[which.max(all_lowest)] # gives the start row number, at the bottom of the list with this error rate
  }
  while (cv.results$size[row.number.max] < 5) {
    row.number.max <- row.number.max -1
  }
  # row.number.max now holds the index for the threshold to use
  threshold <- cv.results$threshold[row.number.max]
  ncpgs <- cv.results$size[row.number.max]
  list.results <- list(threshold = threshold, ncpgs = ncpgs)
  return(list.results)
}
  
store.cpg_ids <- function(df.results, cpg_ids, row.number) {
  if (!missing(cpg_ids)) {
    df.results$CpG_IDs[[row.number]] <- cpg_ids # stores a vector of indices in the CpG_IDs 
  }
  return(df.results)
}

store.numeric.results <- function(df.results, threshold, N_CpGs, error.rate.conf, error.predict = NULL, row.number){
  df.results$threshold[row.number] <- threshold
  df.results$N_CpGs[row.number] <- N_CpGs
  df.results$error_rate_conf[row.number] <- error.rate.conf
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
  list.results <- get.threshold.ncpgs(cv.results)
  threshold <- list.results$threshold
  N_CpGs <- list.results$ncpgs
  # to get the indices of the cpg sites for threshold = threshold:
  cpg_ids <- pamr.predict(trained.data, training.list$x, threshold = threshold, type = "nonzero") # a list of indices for the CpG sites used
  df.results <- store.cpg_ids(df.results, cpg_ids, row.number) # stores the indices of the cpg sites
  conf.table <- as.data.frame(pamr.confusion(cv.results, threshold = threshold, extra = F))
  # the "extra" argument allows for output to be stored in variable if set to false
  error.rate.conf <- sum(conf.table$Freq[2:3])/sum(conf.table$Freq) # errors divided by total samples
  
  predicted <- pamr.predict(trained.data, test.list$x, threshold = threshold)
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
  
  df.results <- store.numeric.results(df.results, threshold, N_CpGs, error.rate.conf, error.predict, row.number)
  
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
assign(paste0("cpg.freq", ".", subtype), cpg.freq)

# After getting cpg.freq for all subtypes, save them for later use:
save(cpg.freq.t821, cpg.freq.inv16, cpg.freq.mono7, cpg.freq.MLL, file = "AML_cpg.freqs.Rdata")

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
