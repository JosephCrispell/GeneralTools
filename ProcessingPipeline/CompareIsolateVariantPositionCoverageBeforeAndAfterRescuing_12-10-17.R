############
# Set path #
############

path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_13-07-17/"

###############################################################
# Read in variant position coverage - before and after rescue #
###############################################################

# Read in the genome coverage file
file <- paste(path, "vcfFiles/IsolateVariantPositionCoverage_27-09-2017.txt", sep="")
before <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

file <- paste(path, "vcfFiles/IsolateVariantPositionCoverage_RESCUED_27-09-2017.txt", sep="")
after <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Parse the Isolate columns
before$Isolate <- parseIsolateColumn(before$Isolate)
after$Isolate <- parseIsolateColumn(after$Isolate)

#############################
# Plot the isolate coverage #
#############################

# Add species column
after$Species <- "COW"
after$Species[grepl(x=after$Isolate, pattern="WB") == TRUE] <- "BADGER"

# Open a pdf
file <- paste(substr(file, 1, nchar(file) - 4), ".pdf", sep="")
pdf(file)

produceSummaryPlots(before, after)

dev.off()