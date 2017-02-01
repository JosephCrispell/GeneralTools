# Script to produce a FASTA alignment to be used in BASTA analysis
## Alignment from main WP clade
## Sampling dates of associated animals
## Sampling locations of associated animals
## Single isolates from the badgers involved (from within clade) - random or best quality

#############
# Libraries #
#############

library(ape)
library(geiger)
library(plotrix)
library(lubridate)

# Current date
date <- format(Sys.Date(), "%d-%m-%y")

############################################
# Get a list of the isolates in main clade #
############################################

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/"

# Read in the newick tree
file <- paste(path, "allVCFs-IncludingPoor/vcfFiles/",
              "mlTree_Prox-10_plusRef_rmResequenced_SNPCov-0.1_28-10-16.tree", sep="")
tree <- read.tree(file=file)

# Plot tree with node numbers
#file <-paste(path, "testTree_nodeNumbers.pdf")
#pdf(file, height=20, width=20)
#plot.phylo(tree, type="fan", show.tip.label=FALSE)
#nodelabels()
#dev.off()

# Get a list of the isolates in the clade
node <- 289
cladeTips <- tips(tree, node=node)
subTree <- extract.clade(tree, node=node)

# Convert this array to list
isolatesInClade <- convertVectorToList(cladeTips)

#############################
# Read in Isolate Sequences #
#############################

# Read in the FASTA file
file <- paste(path, "allVCFs-IncludingPoor/vcfFiles/",
              "sequences_Prox-10_plusRef_rmResequenced_SNPCov-0.1_28-10-16.fasta", sep="")
sequences <- readFASTA(file, skip=1)

# Get the sequences for the isolates
sequenceNames <- cladeTips
isolateSequences <- c()
for(i in 1:length(cladeTips)){
  
  # Check if sequence present for isolate
  if(is.null(sequences[cladeTips[i]]) == FALSE){
    
    isolateSequences[i] <- as.character(sequences[cladeTips[i]])
  }else{
    print(paste("Couldn't find sequence for: ", cladeTips[i]))
  }
}

######################################
# Get the Sampling Date and Location #
######################################

# Cattle Isolates
file <- paste(path, "IsolateData/",
              "CattleIsolateInfo_LatLongs_plusID_outbreakSize_Coverage_AddedTB1453-TB1456.csv", sep="")
cattleInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Badger Isolates
file <- paste(path, "IsolateData/",
              "BadgerInfo_08-04-15_LatLongs_XY.csv", sep="")
badgerInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

tipInfo <- as.data.frame(matrix(nrow=length(cladeTips), ncol=5))
colnames(tipInfo) <- c("IsolateID", "SamplingDate", "X", "Y", "AnimalID")
tipInfo[, "IsolateID"] <- cladeTips

for(tipIndex in 1:length(cladeTips)){
  
  # Cattle
  if(grepl(pattern="TB", x=cladeTips[tipIndex]) == TRUE){
    
    # Find index in table
    strainIndex <- which(cattleInfo$StrainId == cladeTips[tipIndex])
    tipInfo[tipIndex, "SamplingDate"] <- strsplit(as.character(cattleInfo[strainIndex, "BreakdownID"]),
                                                  split="-")[[1]][2] # 14082000501-23/02/1999
    tipInfo[tipIndex, "X"] <- cattleInfo[strainIndex, "Mapx"]
    tipInfo[tipIndex, "Y"] <- cattleInfo[strainIndex, "Mapy"]
    tipInfo[tipIndex, "AnimalID"] <- cattleInfo[strainIndex, "Rawtag"] 
    
    # Badgers
  }else if(grepl(pattern="WB", x=cladeTips[tipIndex]) == TRUE){
    
    # Find index in table
    strainIndex <- which(badgerInfo$WB_id == cladeTips[tipIndex])
    tipInfo[tipIndex, "SamplingDate"] <- badgerInfo[strainIndex, "date"] # 12/01/2000
    tipInfo[tipIndex, "X"] <- badgerInfo[strainIndex, "SampledGrpX"]
    tipInfo[tipIndex, "Y"] <- badgerInfo[strainIndex, "SampledGrpY"]
    tipInfo[tipIndex, "AnimalID"] <- badgerInfo[strainIndex, "tattoo"]
  }
}

# Add species column
tipInfo$Species <- rep("BADGER", nrow(tipInfo))
tipInfo[grepl(pattern="TB", tipInfo$IsolateID), "Species"] <- "COW"

###################################
# Get the isolate genome coverage #
###################################

# Read in genome coverage table
file <- paste(path, "allVCFs-IncludingPoor/vcfFiles/",
              "isolateGenomeCoverageSummary_28-10-16.txt", sep="")
coverageInfo <- read.table(file, header=TRUE, stringsAsFactors=FALSE)
coverageInfo$IsolateID <- getIsolateIDFromFileNames(coverageInfo$IsolateID)

# Add coverage column to tip info
tipInfo$Coverage <- rep(0, nrow(tipInfo))

for(row in 1:nrow(coverageInfo)){
  
  # Is this an isolate we're interested in?
  if(is.null(isolatesInClade[[coverageInfo[row, "IsolateID"]]]) == FALSE){
    tipInfo[isolatesInClade[[coverageInfo[row, "IsolateID"]]], "Coverage"] <- coverageInfo[row, "PercentageCoverage"]
  }
}

####################################################
# Pick Isolates for badgers sampled multiple times #
####################################################

# Add sequences to table
tipInfo$Sequences <- isolateSequences

# Initialise a list to record which badger isolates are kept
rowsToKeep <- c()
index <- 0

# Examine each of the sampled badgers
uniqueIDs <- unique(tipInfo$AnimalID)
for(animal in uniqueIDs){
  
  # Get the row indices for sampled animal
  rowIndices <-  which(tipInfo$AnimalID == animal)
  
  # Check if multiple sequences available
  if(length(rowIndices) > 1){
    
    # Choose an isolate from the available
    #chosenIndex <- sample(rowIndices, size=1)
    chosenIndex <- rowIndices[which.max(tipInfo[rowIndices, "Coverage"])]
    index <- index + 1
    rowsToKeep[index] <- chosenIndex
    
  # Keep single isolate for sampled animal
  }else{
    index <- index + 1
    rowsToKeep[index] <- rowIndices[1]
  }
}

# Keep isolate info for those selected
selectedIsolatesInfo <- tipInfo[rowsToKeep, ]

################
# Define Demes #
################

# Define inner circle radius
thresholdDistance <- 3000
outerDistance <- 6500

# Note Mansion Location
mansionX <- 380909
mansionY <- 201377

# Set an expansion
expand <- 7000

# Open a PDF
file <- paste(path, "BASTA/DefinedDemes_", date, ".pdf", sep="")
pdf(file)

# Create empty plot
plot(1, type="n", yaxt="n", xaxt="n",
     xlim=c(mansionX - expand, mansionX + expand),
     ylim=c(mansionY - expand, mansionY + expand),
     xlab=paste((expand * 2) / 1000, "km"), ylab=paste((expand * 2) / 1000, "km"),
     main="Isolate Locations")

legend("topleft", legend=c("Cow", "Badger", "Woodchester Mansion"),
       pch=c(17, 16, 15), cex=0.65, bty='n')

# Add point for Woodchester Mansion
points(x=mansionX, y=mansionY, pch=15, col="black")

# Add Selected Isolate Locations
points(x=selectedIsolatesInfo$X, y=selectedIsolatesInfo$Y, 
       pch=ifelse(selectedIsolatesInfo$Species == "COW", 17, 16),
       col=ifelse(selectedIsolatesInfo$Species == "COW", rgb(0,0,1, 0.5), rgb(1,0,0, 0.5)))

# Find centre point of badger locations
badgerCentre <- c()
badgerCentre[1] <- mean(selectedIsolatesInfo[selectedIsolatesInfo$Species == "BADGER", "X"], na.rm=TRUE)
badgerCentre[2] <- mean(selectedIsolatesInfo[selectedIsolatesInfo$Species == "BADGER", "Y"], na.rm=TRUE)
draw.circle(x=badgerCentre[1], y=badgerCentre[2], radius=thresholdDistance, border="black")
draw.circle(x=badgerCentre[1], y=badgerCentre[2], radius=outerDistance, border="black")
text(x=381938.4, y=198779.4, labels=paste(paste(thresholdDistance / 1000, "km")))
text(x=381938.4, y=194594, labels=paste(paste(outerDistance / 1000, "km")))

# Note all isolates that are within 3km of the badger centre point
selectedIsolatesInfo$Distance <- rep(NA, nrow(selectedIsolatesInfo))
for(row in 1:nrow(selectedIsolatesInfo)){
  
  # Skip isolates with an unknown location
  if(is.na(selectedIsolatesInfo[row, "X"]) == FALSE){
    selectedIsolatesInfo[row, "Distance"] <- euclideanDistance(x1=badgerCentre[1], y1=badgerCentre[2],
                                                               x2=selectedIsolatesInfo[row, "X"],
                                                               y2=selectedIsolatesInfo[row, "Y"])
  }  
}
plot(selectedIsolatesInfo$Distance,
     col=ifelse(selectedIsolatesInfo$Species == "COW", rgb(0,0,1, 0.75),
                rgb(1,0,0, 0.75)),
     xaxt="n", xlab="", ylab="Distance (m)", main="Distance to Badger Centre Point", las=1, pch=19)
legend("top", legend=c("Cow", "Badger"),
       text.col=c("blue", "red"), cex=1, bty='n')
points(x=c(0, nrow(selectedIsolatesInfo)), y=c(thresholdDistance, thresholdDistance), type="l", lty=2, col="black")
dev.off()

# Define demes
# A - anything within Xkm of badger centre
# B - everything else

selectedIsolatesInfo$Deme <- rep("NA", nrow(selectedIsolatesInfo))
for(row in 1:nrow(selectedIsolatesInfo)){
  
  # Is location data available?
  if(is.na(selectedIsolatesInfo[row, "X"]) == FALSE){
    
    # Was the isolate sampled within 3km?
    if(selectedIsolatesInfo[row, "Distance"] <= thresholdDistance){
      
      if(selectedIsolatesInfo[row, "Species"] == "BADGER"){
        selectedIsolatesInfo[row, "Deme"] <- "innerBadger"
      }else{
        selectedIsolatesInfo[row, "Deme"] <- "innerCow"
      }      
    }else{
      selectedIsolatesInfo[row, "Deme"] <- "outerCow"
    }
    
  # Is it a badger?
  }else if(selectedIsolatesInfo[row, "Species"] == "BADGER"){
    selectedIsolatesInfo[row, "Deme"] <- "innerBadger"
  }
}

# Create a plot to demonstrate rates being estimated
states <- c("innerCow", "outerCow", "innerBadger", "outerBadger")
rates <- matrix(data=c(NA, "innerCow->outerCow", "innerCow->innerBadger", "-",
                       "outerCow->innerCow", NA, "-", "outerCow->outerCow",
                       "innerBadger->innerCow", "-", NA, "innerBadger->outerBadger",
                       "-", "outerBadger->outerCow", "outerBadger->innerBadger", NA),
                ncol=4, nrow=4, byrow=TRUE)




######################
# Print out sequence #
######################

#<data id="alignmentVar" dataType="nucleotide">
#  <sequence taxon="human">
#  TAGAACTG
#</sequence>
#  <sequence taxon="chimp">
#  ACAAACTG
#</sequence>
#  </data>
#  <data id='alignment' spec='FilteredAlignment' filter='-' data='@alignmentVar' constantSiteWeights='1005 1043 2049 1003'/>

# Create a vector to store the fileLines
fileLines <- c()

# Attach Sampling Information to name: ID_Species_Deme
fileLines[length(fileLines) + 1] <- "\t<!-- Sequence Alignment -->"
fileLines[length(fileLines) + 1] <- "\t<data id=\"alignmentVar\" dataType=\"nucleotide\">"

for(row in 1:nrow(selectedIsolatesInfo)){
  
  name <- paste(selectedIsolatesInfo[row, "IsolateID"], selectedIsolatesInfo[row, "Species"],
                selectedIsolatesInfo[row, "Deme"], round(selectedIsolatesInfo[row, "Coverage"], digits=2), sep="_")
  fileLines[length(fileLines) + 1] <- paste("\t\t<sequence taxon=\"", name, "\">", sep="")
  fileLines[length(fileLines) + 1] <- paste("\t\t\t", selectedIsolatesInfo[row, "Sequences"], sep="")
  fileLines[length(fileLines) + 1] <- "\t\t</sequence>"
}
fileLines[length(fileLines) + 1] <- "\t</data>"

########################
# Constant Site Counts #
########################

counts <- "700377 1322380 1316990 699795"
fileLines[length(fileLines) + 1] <- ""
fileLines[length(fileLines) + 1] <- "\t<!-- Constant Site Counts -->"
fileLines[length(fileLines) + 1] <- paste("\t<data id='alignment' spec='FilteredAlignment' filter='-' data='@alignmentVar' constantSiteWeights='",
                                          counts, "'/>", sep="")  

######################
# Add Sampling Dates #
######################

# Convert the sampling dates to date objects - 23/02/1999
selectedIsolatesInfo$SamplingDate <- as.Date(selectedIsolatesInfo$SamplingDate,
                                             format="%d/%m/%Y")

# Create a decimal date column
selectedIsolatesInfo$DecimalDate <- decimal_date(selectedIsolatesInfo$SamplingDate)

# Print out tip dates
fileLines[length(fileLines) + 1] <- ""
fileLines[length(fileLines) + 1] <- "\t<!-- Tip Dates -->"
fileLines[length(fileLines) + 1] <- "\t<timeTraitSet spec='beast.evolution.tree.TraitSet' id='timeTraitSet' traitname=\"date-forward\""
output <- "\t\tvalue=\""
for(row in 1:nrow(selectedIsolatesInfo)){
  name <- paste(selectedIsolatesInfo[row, "IsolateID"], selectedIsolatesInfo[row, "Species"],
                selectedIsolatesInfo[row, "Deme"], round(selectedIsolatesInfo[row, "Coverage"], digits=2), sep="_")

  output <- paste(output, name, "=", selectedIsolatesInfo[row, "DecimalDate"], sep="")
  
  if(row < nrow(selectedIsolatesInfo)){
    output <- paste(output, ",", sep="")
  }
}
output <- paste(output, "\">", sep="")
fileLines[length(fileLines) + 1] <- output
fileLines[length(fileLines) + 1] <- "\t\t<taxa spec='TaxonSet' alignment='@alignment'/>"
fileLines[length(fileLines) + 1] <- "\t</timeTraitSet>"

##################################
# Print out the Deme Information #
##################################

# Four demes:
#   0   innerBadger   Sampled
#   1   innerCow      Sampled
#   2   outerCow      Sampled
#   3   outerBadger   UnSampled
fileLines[length(fileLines) + 1] <- ""
fileLines[length(fileLines) + 1] <- "\t<!-- Deme Assignment -->"
output <- "\t<typeTraitSet id=\"typeTraitSet\" spec=\"TraitSet\" traitname=\"type\" value=\""
for(row in 1:nrow(selectedIsolatesInfo)){
  
  name <- paste(selectedIsolatesInfo[row, "IsolateID"], selectedIsolatesInfo[row, "Species"],
                selectedIsolatesInfo[row, "Deme"], round(selectedIsolatesInfo[row, "Coverage"], digits=2), sep="_")
    
  output <- paste(output, name, "=", selectedIsolatesInfo[row, "Deme"], sep="")
  
  if(row < nrow(selectedIsolatesInfo)){
    output <- paste(output, ",", sep="")
  }
}
output <- paste(output, "\">", sep="")
fileLines[length(fileLines) + 1] <- output
fileLines[length(fileLines) + 1] <- "\t\t<taxa spec='TaxonSet' alignment='@alignment'/>"
fileLines[length(fileLines) + 1] <- "\t</typeTraitSet>"
fileLines[length(fileLines) + 1] <- ""

# Open an output file
fileName <- paste(path, "BASTA/", "xmlInput_", date, ".xml", sep="")
fileConnection <- file(fileName)

# Print out file lines
writeLines(fileLines, fileConnection)

# Close the output file
close(fileConnection)

#############
# FUNCTIONS #
#############

euclideanDistance <- function(x1, y1, x2, y2){
  return(sqrt(sum((x1 - x2)^2 + (y1 - y2)^2)))
}

getIsolateIDFromFileNames <- function(fileNames){
  isolates <- c()
  for(i in 1:length(fileNames)){
    isolates[i] <- strsplit(fileNames[i], split="_")[[1]][1]
  }
  
  return(isolates)
}

readFASTA <- function(fileName, skip){

  # Open file and read in lines
  connection  <- file(fileName, open = "r")
  fileLines <- readLines(connection, warn=TRUE) 
  close(connection)
  
  # Initialise a list to store the isolate sequences
  sequences <- list()
  
  # Parse every line in the file
  for(lineIndex in (1+skip):length(fileLines)){
    
    line <- fileLines[lineIndex]
    
    # Have we found a new sequence?
    if(grepl(pattern=">", line) == TRUE){
      
      # Store the previous sequence
      if(lineIndex != (1+skip)){
        sequences[[name]] <- sequence
      }
      
      # Get the isolate name and initialise sequence - >WB167_S3_69.vcf
      name <- substr(line, start=2, stop=nchar(line))
      name <- strsplit(name, spli="_")[[1]][1]
      sequence <- ""
    }else{
      sequence <- paste(sequence, line, sep="")
    }
  }
  
  return(sequences)
}

convertVectorToList <- function(vector){
  output <- list()
  for(index in 1:length(vector)){
    output[[vector[index]]] <- index
  }
  return(output)
}
