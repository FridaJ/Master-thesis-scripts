# This script will make line charts of the cpg.freq results for the diagnostic subgroups. 

setwd("Documents/R_work")

# Define one subtype vector
MLL <- cpg.freq.MLL$freq

# Graph MLL using blue points overlayed by a line 
plot(MLL, type="o", col="blue", ylim = c(0,25), axes=FALSE, ann=FALSE)

# Create a title with a red, bold/italic font
title(main="CpG sites in 5x5 CV, MLL", col.main="red", font.main=4)

# Make x axis using CpG site labels
#axis(1, at=1:nrow(cpg.freq.MLL), lab=cpg.freq.MLL$index)
axis(1, las=1)

# Make y axis with horizontal labels.
axis(2, las=1)

# Create box around plot
box()

# Label the x and y axes with dark green text
title(xlab="CpG sites chosen for classifier", col.lab=rgb(0,0.5,0))
title(ylab="No of CV runs including the CpG site", col.lab=rgb(0,0.5,0))
