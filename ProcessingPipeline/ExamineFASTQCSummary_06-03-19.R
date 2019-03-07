#### Read in the data ####

# Set the path variable
path <- "/home/josephcrispell/Desktop/Research/Cumbria/FASTQs_06-03-19/FastQC/"

# Read in the FASTQC file summary table
summaryFile <- paste0(path, "Fastqc_summary_06-03-19.txt")
summary <- read.table(summaryFile, header=TRUE, sep="\t")

# Look at GC distribution
hist(summary$GC, breaks=20, las=1)

# Look at the left and right trimming
boxplot(summary$LeftTrim, summary$RightTrim, names=c("LEFT", "RIGHT"), las=1, frame=FALSE, main="Trimming suggestions")

# Check adapter flag
table(summary$AdapterContentFlag)
