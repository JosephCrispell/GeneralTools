############
# Set path #
############

path <- "/home/josephcrispell/Desktop/Research/Cumbria/vcfFiles/"

##########################
# Read in coverage files #
##########################

# # Read in the genome coverage file
# genomeCoverageFile <- paste(path, "genomeCoverageSummary_DP-20_21-07-2017.txt", sep="")
# genomeCoverage <- read.table(genomeCoverageFile, header=TRUE, stringsAsFactors=FALSE)

# Read in the isolate coverage file
isolateCoverageFile <- paste(path, "isolateCoverageSummary_DP-20_06-03-2019.txt", sep="")
isolateCoverage <- read.table(isolateCoverageFile, header=TRUE, stringsAsFactors=FALSE)

# Parse the Isolate column
isolateCoverage$IsolateID <- parseIsolateColumn(isolateCoverage$IsolateID)

#############################
# Plot the isolate coverage #
#############################

# Open a pdf
file <- paste(substr(isolateCoverageFile, 1, nchar(isolateCoverageFile) - 4), ".pdf", sep="")
pdf(file)

plot(y=isolateCoverage$PercentageCoverage,
     x=isolateCoverage$MeanDepth,
     las=1, ylab="Proportion", main="Proportion of M. bovis Genome with >19 mapped reads",
     xlab="Mean Read Depth", pch=16, cex=3,
     col=ifelse(grepl(x=isolateCoverage$IsolateID, pattern="WB"), rgb(1,0,0, 0.5),
                rgb(0,0,1, 0.5)))

text(y=isolateCoverage$PercentageCoverage,
     x=isolateCoverage$MeanDepth,
     labels = isolateCoverage$IsolateID, cex=1, pos=4,
     col=ifelse(isolateCoverage$PercentageCoverage < 0.8,
                rgb(0,0,0, 1), rgb(0,0,0, 0)))

# legend("bottomright", legend=c("BADGER", "COW"), text.col=c("red", "blue"), bty="n")


# subset <- isolateCoverage[which(grepl(pattern="-", isolateCoverage$IsolateID)), ]
# 
# plot(y=subset$PercentageCoverage,
#      x=subset$MeanDepth,
#      las=1, ylab="Proportion", 
#      main="Proportion of New Cattle with >19 mapped reads",
#      xlab="Mean Read Depth", pch=16, cex=3,
#      col=ifelse(grepl(x=subset$IsolateID, pattern="WB"), rgb(1,0,0, 0.5),
#                 rgb(0,0,1, 0.5)))
# 
# text(y=subset$PercentageCoverage,
#      x=subset$MeanDepth,
#      labels = subset$IsolateID, cex=1, pos=4,
#      col=ifelse(subset$PercentageCoverage < 0.8,
#                 rgb(0,0,0, 1), rgb(0,0,0, 0)))
# 
# legend("bottomright", legend=c("BADGER", "COW"), text.col=c("red", "blue"), bty="n")


dev.off()

############################
# Plot the genome coverage #
############################

# # Note reference genome size
# MbovisSize <- 4349904
# 
# # Open a pdf
# file <- paste(substr(genomeCoverageFile, 1, nchar(genomeCoverageFile) - 4), ".pdf", sep="")
# pdf(file, height=7, width=7)
# 
# plot(y=genomeCoverage$MeanDepth, x=genomeCoverage$POS, type="o",
#      xlim=c(1, MbovisSize), ylim=c(0, max(genomeCoverage$MeanDepth)),
#      bty="n", ylab="Average Read Depth", xlab="Genome Position", las=1,
#      main="Genome Coverage Across Isolates")
# 
#dev.off()


#############
# FUNCTIONS #
#############

parseIsolateColumn <- function(column){
  
  ids <- c()
  for(i in 1:length(column)){
    parts <- strsplit(column[i], split="_")[[1]]
    ids[i] <- paste(parts[1], "_", parts[2], sep="")
  }
  
  return(ids)
}
