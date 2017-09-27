path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/"


########################################
# Read in a list of sequenced isolates #
########################################

file <- paste(path, "allVCFs-IncludingPoor/vcfFiles/", "sequencedIsolates_28-10-16.txt", sep="")
sequencedIsolates <- read.table(file, header=FALSE, stringsAsFactors=FALSE)[, 1]

########################################################
# Get a list of resequenced isolates that were removed #
########################################################

file <- paste(path, "allVCFs-IncludingPoor/vcfFiles/", "removedResequencedIsolates_28-10-16.txt", sep="")
resequencedIsolatesThatWereRemoved <- read.table(file, header=FALSE, stringsAsFactors=FALSE)[, 1]

################################################################
# Get list of 230 badgers and 100 cattle sequenced isolate IDs #
################################################################

sequencedIsolates <- sequencedIsolates[sequencedIsolates %in% resequencedIsolatesThatWereRemoved == FALSE]


############################################
# Read in the Isolate Coverage Information #
############################################

# Get the genome coverage for all the isolates
file <- paste(path, "allVCFs-IncludingPoor/", "isolateGenomeCoverageSummary_07-09-16.txt", sep="")
coverageTable <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep="\t")
coverageTable$IsolateID <- getIsolateIDFromVCFFileNames(coverageTable$IsolateID)

# Keep the sequenced isolate information
coverageTable <- coverageTable[coverageTable$IsolateID %in% sequencedIsolates, ]

# Order the table by isolate ID and genome coverage
coverageTable <- coverageTable[order(coverageTable$PercentageCoverage, decreasing=TRUE), ]
coverageTable <- rbind(coverageTable[grepl(x=coverageTable$IsolateID, pattern="WB") == TRUE, ],
                       coverageTable[grepl(x=coverageTable$IsolateID, pattern="WB") == FALSE, ])

####################################
# Plot the isolate genome coverage #
####################################
file <- paste(path, "allVCFs-IncludingPoor/", "isolateGenomeCoverageSummary_07-09-16.pdf", sep="")
pdf(file)

coverageTable$Species <- "COW"
coverageTable$Species[grepl(x=coverageTable$IsolateID, pattern="WB") == TRUE] <- "BADGER"
coverageTable$Species <- as.factor(coverageTable$Species)

par(mar=c(5,4.1,4.1,2.1)) # Bottom, left, top and right
boxplot(coverageTable$PercentageCoverage ~ coverageTable$Species, 
        ylim=c(0,1), ylab="Proportion", border=c("red", "blue"),
        names=c("Badgers", "Cattle"), outline=FALSE,
        las=1, pch=20, main="Isolate Genome Coverage")

stripchart(coverageTable$PercentageCoverage ~ coverageTable$Species,
           vertical = TRUE, jitter=0.2,
           method = "jitter", add = TRUE, pch = 21,
           col = c(rgb(1,0,0, 0.5), rgb(0,0,1, 0.5)),
           bg=rgb(0.5,0.5,0.5, 0.5))



dev.off()

#############
# FUNCTIONS #
#############

getIsolateIDFromVCFFileNames <- function(vcfFileNames){
  
  output <- c()
  for(i in 1:length(vcfFileNames)){
    
    output[i] <- strsplit(vcfFileNames[i], split="_")[[1]][1]
  }
  
  return(output)
  
}
