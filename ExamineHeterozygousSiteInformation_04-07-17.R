####################
# Read in the data #
####################

# Set directory path
#path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/allVCFs-IncludingPoor/vcfFiles/"
#path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/NewZealand/NewAnalyses_12-05-16/vcfFiles/"
#path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/allVCFs-IncludingPoor/vcfFiles/PoorlyMappedBadgers/"
#path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/allVCFs-IncludingPoor/vcfFiles/MislabelledBadgers/"
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Cumbria/vcfFiles/"
newZealand <- grepl(x=path, pattern="NewZealand")

# Read in the variant site coverage and alternate allele support files
coverageFile <- paste(path, "heterozygousSiteInfo_Coverage_10-09-2017.txt", sep="")
siteInfoCoverage <- read.table(coverageFile, header=TRUE, stringsAsFactors=FALSE,
                               comment.char="~", check.names=FALSE)
altSupportFile <- paste(path, "heterozygousSiteInfo_AltSupport_10-09-2017.txt", sep="")
siteInfoAltSupport <- read.table(altSupportFile, header=TRUE, stringsAsFactors=FALSE,
                                 comment.char="~", check.names=FALSE)

# Read in the genome coverage information
genomeCoverageFile <- paste(path, "isolateCoverageSummary_DP-20_10-09-2017.txt", sep="")
genomeCoverage <- read.table(genomeCoverageFile, header=TRUE, stringsAsFactors=FALSE,
                             comment.char="~", check.names=FALSE)

# Parse the VCF file names in the columns
isolateIDs <- getIDsFromFileNames(colnames(siteInfoAltSupport)[-1], newZealand)
colnames(siteInfoCoverage) <- c("Positive", isolateIDs)
colnames(siteInfoAltSupport) <- c("Positive", isolateIDs)
genomeCoverage$IsolateID <- isolateIDs

#################
# Plot the data #
#################

# Set threshold of acceptance
supportThreshold <- 0.05
depthThreshold <- 4

# Create an array to note how many heterozygous sites there are for each isolate
nHeterozygousSites <- c()

# Open a pdf
file <- paste(path, "heterozygousSiteInfo_25-07-2017.pdf", sep="")
pdf(file)

for(i in 1:length(isolateIDs)){
  
  # Note the ID of the current isolate
  isolate <- isolateIDs[i]

  # Note the rows containing Heterozygous site info for the current isolate
  badRows <- c(1:nrow(siteInfoAltSupport))[is.na(siteInfoAltSupport[, isolate]) == FALSE &
                                           siteInfoAltSupport[, isolate] > 0 + supportThreshold &
                                           siteInfoAltSupport[, isolate] < 1 - supportThreshold &
                                           siteInfoCoverage[, isolate] >= depthThreshold]
  # Note the heterozygous sites found outwith thresholds
  nHeterozygousSites[i] <- length(badRows)
  
  # Plot the heterozygous site info
  plot(y=siteInfoAltSupport[-badRows, isolate], 
       x=siteInfoCoverage[-badRows, isolate],
       pch=20, las=1, bty="n", cex=2, ylim=c(0,1), xlim=c(0, 100),
       col=ifelse(siteInfoCoverage[-badRows, isolate] >= depthThreshold,
                  rgb(0,0,1, 0.5), rgb(0,0,0, 0.5)), 
       ylab="Proportion HQ Reads Supporting Alternate",
       xlab="High Quality Read Depth", main=isolate)
  polygon(x=c(depthThreshold,100,100,depthThreshold), 
          y=c(0, 0, 0 + supportThreshold, 0 + supportThreshold),
          border=rgb(0,0,1, 0.5), col=rgb(0,0,1, 0.5))
  polygon(x=c(depthThreshold,100,100, depthThreshold),
          y=c(1, 1, 1 - supportThreshold, 1 - supportThreshold),
          border=rgb(0,0,1, 0.5), col=rgb(0,0,1, 0.5))
  polygon(x=c(0,0,depthThreshold,depthThreshold), 
          y=c(0, 1, 1, 0),
          border=rgb(0,0,0, 0.5), col=rgb(0,0,0, 0.5))
  points(y=siteInfoAltSupport[badRows, isolate], 
         x=siteInfoCoverage[badRows, isolate],
         pch=20, cex=2, 
         col=rgb(1,0,0, 0.5))
  legend(x=70, y=0.9, bty="n", cex=0.75,
         legend=paste("n = ", nHeterozygousSites[i], "\nCoverage = ", 
                      round(genomeCoverage[which(genomeCoverage$IsolateID == isolate),
                                           "PercentageCoverage"], digits=2), sep=""))
}


plot(x=nHeterozygousSites, y=genomeCoverage$PercentageCoverage,
     pch=20, bty="n", las=1,
     col=ifelse(grepl(x=isolateIDs, pattern="WB") == TRUE, 
                rgb(1,0,0, 0.5), rgb(0,0,1, 0.5)),
     ylab="Proportion of genome with >= 20 aligned reads",
     xlab="Number of heterozygous sites")
text(x=nHeterozygousSites, y=genomeCoverage$MeanDepth,
     labels = genomeCoverage$IsolateID, cex=0.5,
     col=ifelse(nHeterozygousSites > 100, rgb(0,0,0, 1), rgb(0,0,0, 0)))
#if(newZealand == FALSE){
#  legend("topright", legend=c("Badger", "Cow"), text.col=c("red", "blue"), bty="n")
#}

plot(x=nHeterozygousSites, y=genomeCoverage$MeanDepth,
     pch=20, bty="n", las=1,
     col=ifelse(grepl(x=isolateIDs, pattern="WB") == TRUE, 
                rgb(1,0,0, 0.5), rgb(0,0,1, 0.5)),
     ylab="Mean read depth on genome",
     xlab="Number of heterozygous sites")
text(x=nHeterozygousSites, y=genomeCoverage$MeanDepth,
     labels = genomeCoverage$IsolateID, cex=0.5,
     col=ifelse(nHeterozygousSites > 100, rgb(0,0,0, 1), rgb(0,0,0, 0)))
#if(newZealand == FALSE){
#  legend("topright", legend=c("Badger", "Cow"), text.col=c("red", "blue"), bty="n")
#}
  
# Note the species of the isolates
#if(newZealand == FALSE){
#  species <- getSpecies(isolateIDs)
#  
#  boxplot(nHeterozygousSites ~ species, 
#          ylab="Number of heterozygous sites", border=c("red", "blue"),
#          names=c("Badgers", "Cattle"), outline=FALSE,
#          las=1, pch=20, ylim=range(nHeterozygousSites))
#  
#  stripchart(nHeterozygousSites ~ species,
#             vertical = TRUE, jitter=0.2,
#             method = "jitter", add = TRUE, pch = 21,
#             col = c(rgb(1,0,0, 0.5), rgb(0,0,1, 0.5)),
#             bg=rgb(0.5,0.5,0.5, 0.5))
#}

dev.off()

#############
# FUNCTIONS #
#############

getSpecies <- function(isolateIDs){
  species <- rep("COW", length(isolateIDs))
  for(i in 1:length(isolateIDs)){
    
    if(grepl(x=isolateIDs[i], pattern="WB") == TRUE){
      species[i] <- "BADGER"
    }
  }
  
  return(species)
}

getIDsFromFileNames <- function(fileNames, newZealand){
  
  ids <- c()
  for(i in 1:length(fileNames)){
    parts <- strsplit(strsplit(fileNames[i], split=".vcf")[[1]][1], split="_")[[1]]
    ids[i] <- paste(parts[1], parts[3], sep="_")
    if(newZealand == TRUE){
      ids[i] <- paste(parts[1], parts[2], parts[3], sep="_")
    }
  }
  
  return(ids)
}
