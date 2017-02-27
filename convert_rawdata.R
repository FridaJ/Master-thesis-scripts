##### Do this in linux environment on UPPMAX #####

#change to working directory
cd /pica/v13/b2015118/private/frida-aml

#start an interactive session with 8 cores (64 GB), 1 hour
interactive -A b2015118 -n 8 -t 1:00:00

#start R
R

##### Do this in R on UPPMAX #########

#Set working directory in R
setwd("/pica/v13/b2015118/private/frida-aml")

#If packages need to be downloaded:
# source("http://bioconductor.org/biocLite.R")

#If packages need to be installed:
# biocLite("minfi")
# biocLite("IlluminaHumanMethylation450kmanifest")
# biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19") 
# "May need to remove dir /pica/h1/fridaj/R/x86_64-redhat-linux-gnu-library/3.3/00LOCK-survival" (rm -rf)
# still didnt work
# Warning messages:
#1: In install.packages(update[instlib == l, "Package"], l, contriburl = contriburl,  :
#installation of package ‘lattice’ had non-zero exit status
#2: In install.packages(update[instlib == l, "Package"], l, contriburl = contriburl,  :
#installation of package ‘Matrix’ had non-zero exit status

#Open packages at the beginning of a session
library(minfi)

#Load data file(s) if they already exist
load(file = "GRSet.Rdata")

#This is done to create a system path to the IDAT file
baseDir <- "/proj/b2015118/private/raw_450k/GLM1"

###### Read the user created sample sheet located in baseDir ######

#creates a data.frame that will be used to read in the data
targets <- read.metharray.sheet(baseDir, pattern="GML1_130115_AML1_samplesheet.csv") # to read in the data in the specified dir?
dim(targets) #to see the number of lines (patients) in the dataset (aml, 168)

###### Create Red/green-channel set ######

RGSet <- read.metharray.exp(targets = targets) #(3 min)
dim(RGSet) #check to see if the data frame has correct dimensions (aml = 168 rows)
save(RGSet, file= "RGSet.Rdata") #save the data to file at this point (1 min)

###### Preprocessing with funnorm to create a genomic/ratio dataset ######

GRSet.funnorm <- preprocessFunnorm(RGSet) #uses the IlluminaHumanMethylation450k packages described above (3 min)
save(GRSet.funnorm, file= "GRSet.Rdata") #saves the data to file to be copied to local computer (30 sec)




