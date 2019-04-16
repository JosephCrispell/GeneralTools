#### Read in the data ####

# Set the path variable
path <- "/home/josephcrispell/Desktop/Research/RepublicOfIreland/Fastqs_ALL_15-03-19/Alicia/FASTQC/"

# Read in the FASTQC file summary table
summaryFile <- paste0(path, "Fastqc_summary_04-04-19.txt")
summary <- read.table(summaryFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Look at GC distribution
hist(summary$GC, breaks=20, las=1)

# Look at the left and right trimming
boxplot(summary$LeftTrim, summary$RightTrim, names=c("LEFT", "RIGHT"), las=1, frame=FALSE,
        main="Trimming suggestions", outcol=rgb(0,0,0, 0.1), pch=19)

# Check adapter flag
table(summary$AdapterContentFlag)
