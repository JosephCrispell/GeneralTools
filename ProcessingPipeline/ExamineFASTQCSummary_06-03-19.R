#### Read in the data ####

# Set the path variable
path <- "/home/josephcrispell/storage/Research/RepublicOfIreland/Mbovis/Monaghan/Fastqs_16-12-19/FASTQC/"

# Read in the FASTQC file summary table
summaryFile <- paste0(path, "Fastqc_summary_16-12-19.txt")
summary <- read.table(summaryFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Sort the file names
summary <- summary[order(summary$FileName), ]

# Look at GC distribution
hist(summary$GC, breaks=20, las=1)
summary[summary$GC != 65, c("FileName", "GC")]

# Look at the left and right trimming
boxplot(summary$LeftTrim, summary$RightTrim, names=c("LEFT", "RIGHT"), las=1, frame=FALSE,
        main="Trimming suggestions", outcol=rgb(0,0,0, 0.1), pch=19)

# Check adapter flag
table(summary$AdapterContentFlag)

# Look at the number of reads distribution
hist(summary$NumberReads, breaks=100, las=1)
summary[summary$NumberReads < 250000, c("FileName", "NumberReads")]

# Look at the read length distribution
table(summary$ReadLength)
