#### Load packages ####
library(ape)
library(phangorn)

#### Set the path variable ####
path <- "/home/josephcrispell/Desktop/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/"

#### Read in the FASTA file ####

# Read in the FASTA file
fastaFile <- paste0(path, "vcfFiles/sequences_withoutHomoplasies_27-03-18.fasta")
sequencesDNAbin <- read.dna(fastaFile, format = "fasta", skip=1) # Note skip misses the first line of file with nSites and nSamples in it

alignment <- as.alignment(sequencesDNAbin)

#### Build the distance matrix ####

# Build the genetic distance matrix
distanceMatrix <- as.matrix(dist.dna(sequencesDNAbin, model="F84")) # F84 is equivalent of HKY model

# Note the names of the sequences used
names <- colnames(distanceMatrix)

#### Read in the metadata #### 

# Note the sampling infor file names
cattleInfoFile <- paste(path, "IsolateData/CattleIsolateInfo_AddedNew_TB1484-TB1490_22-03-18.csv", sep="")
badgerInfoFile <- paste(path, "IsolateData/BadgerInfo_08-04-15_LatLongs_XY_Centroids.csv", sep="")

# Get Isolate sampling information
sequenceInfo <- getIsolateSamplingInformation(cattleInfoFile, badgerInfoFile, isolates=names)

#### FUNCTIONS ####

getIsolateSamplingInformation <- function(cattleInfoFile, badgerInfoFile, isolates){
  
  # Read in sampling information
  cattleInfo <- read.table(cattleInfoFile, header=TRUE, sep=",", stringsAsFactors=FALSE)
  badgerInfo <- read.table(badgerInfoFile, header=TRUE, sep=",", stringsAsFactors=FALSE)
  
  # Initialise a table to store the isolate sampling information
  isolateInfo <- as.data.frame(matrix(nrow=length(isolates), ncol=6))
  colnames(isolateInfo) <- c("IsolateID", "SamplingDate", "X", "Y", "AnimalID", "Species")
  isolateInfo[, "IsolateID"] <- isolates
  
  # Fill the table with the sampling information
  for(index in 1:length(isolates)){
    
    # Get the Isolate ID
    id <- strsplit(isolates[index], split="_")[[1]][1]
    
    # Cattle
    if(grepl(pattern="TB|AF-|HI-", x=id) == TRUE){
      
      # Find index in table
      strainIndex <- which(cattleInfo$StrainId == id)
      isolateInfo[index, "SamplingDate"] <- strsplit(as.character(cattleInfo[strainIndex, "BreakdownID"]),
                                                     split="-")[[1]][2] # 14082000501-23/02/1999
      isolateInfo[index, "X"] <- cattleInfo[strainIndex, "Mapx"]
      isolateInfo[index, "Y"] <- cattleInfo[strainIndex, "Mapy"]
      isolateInfo[index, "AnimalID"] <- cattleInfo[strainIndex, "Rawtag"]
      isolateInfo[index, "Species"] <- "COW"
      
      # Badgers
    }else if(grepl(pattern="WB", x=id) == TRUE){
      
      # Find index in table
      strainIndex <- which(badgerInfo$WB_id == id)
      isolateInfo[index, "SamplingDate"] <- badgerInfo[strainIndex, "date"] # 12/01/2000
      if(is.na(badgerInfo[strainIndex, "GroupCentroidX"]) == FALSE){
        isolateInfo[index, "X"] <- badgerInfo[strainIndex, "GroupCentroidX"]
        isolateInfo[index, "Y"] <- badgerInfo[strainIndex, "GroupCentroidY"]
      }else{
        isolateInfo[index, "X"] <- badgerInfo[strainIndex, "SampledGrpX"]
        isolateInfo[index, "Y"] <- badgerInfo[strainIndex, "SampledGrpY"]
      }
      isolateInfo[index, "AnimalID"] <- badgerInfo[strainIndex, "tattoo"]
      isolateInfo[index, "Species"] <- "BADGER"
    }
  }
  
  # Format the sampling dates
  isolateInfo$SamplingDate <- as.Date(isolateInfo$SamplingDate, format="%d/%m/%Y")
  
  return(isolateInfo)
}
