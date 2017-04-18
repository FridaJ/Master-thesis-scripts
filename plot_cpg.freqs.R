# This script will make line charts of the cpg.freq results for the diagnostic subgroups. 
# Set the variables subtype and aml.class to make a plot

setwd("Documents/R_work")
load(file = "AML_cpg.freqs.new.Rdata")

# Define one subtype frequency vector and a string for class
subtype <- cpg.freq.t821.new$freq
aml.class <- "t821"

# Graph subtype using blue points overlayed by a line 
plot(subtype, type="o", col="blue", ylim = c(0,25), axes=FALSE, ann=FALSE)

# Create a title with a red, bold/italic font
title(main=paste0("CpG sites in 5x5 CV, ", aml.class), col.main="red", font.main=4)

# Make x axis using CpG site labels
#axis(1, at=1:length(subtype), lab=cpg.freq.MLL$index)
axis(1, las=1)

# Make y axis with horizontal labels.
axis(2, las=1)

# Create box around plot
box()

# Label the x and y axes with dark green text
title(xlab="CpG sites chosen for classifier", col.lab=rgb(0,0.5,0))
title(ylab="No of CV runs including the CpG site", col.lab=rgb(0,0.5,0))
