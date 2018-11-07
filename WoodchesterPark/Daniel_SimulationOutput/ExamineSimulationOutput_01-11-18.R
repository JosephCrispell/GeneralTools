#### Load libraries ####

library(rgdal) # Converting from latlongs to eastings and northings -installation requires: sudo apt-get install gdal-bin proj-bin libgdal-dev libproj-dev
library(plotrix) # For drawing circles
library(lubridate) # Convert dates to decimal numbers

#### Read in the simulation output ####

# Current date
date <- format(Sys.Date(), "%d-%m-%y")

# Set the path variable
path <- "/home/josephcrispell/Desktop/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/BASTA/Simulations/"

# NOTE the badger centre
badgerCentre <- c(381761.7, 200964.3)

# NOTE the INNER threshold distance
inner <- 3500

# Read in a simulation output table
simFile <- paste0(path, "05-11-18/", "sim_14-01.csv")
simOutput <- read.table(simFile, header=TRUE, sep=",", stringsAsFactors=FALSE)
simOutput$animal_ID <- gsub("_", "-", simOutput$animal_ID)


# Read in the herd and sett coordinates - NOTE: latlongs are actually eastings and northings
herdFile <- paste0(path, "coordinates_unit-ID_HERDS_01-11-18.csv")
herdCoords <- read.table(herdFile, header=TRUE, sep=",", stringsAsFactors=FALSE)
colnames(herdCoords)[c(3,4)] <- c("X", "Y")
settFile <- paste0(path, "coordinates_unit-ID_SETTS_01-11-18.csv")
settCoords <- read.table(herdFile, header=TRUE, sep=",", stringsAsFactors=FALSE)
colnames(settCoords)[c(3,4)] <- c("X", "Y")

# Add the unit coordinate data into the simulation output table
simOutput <- addUnitCoords(simOutput, herdCoords, settCoords, badgerCentre)

# Convert the date columns to dates
simOutput$DateSampleTaken <- as.Date(simOutput$DateSampleTaken, format="%Y-%m-%d")
simOutput$InfectionDate <- as.Date(simOutput$InfectionDate, format="%Y-%m-%d")
simOutput$DetectionDate <- as.Date(simOutput$DetectionDate, format="%Y-%m-%d")
simOutput$LastSnpGeneration <- as.Date(simOutput$LastSnpGeneration, format="%Y-%m-%d")

#### Examine the simulation output ####

# Open a pdf file
plotsFile <- paste0(substr(simFile,1,nchar(simFile)-3), "pdf")
pdf(plotsFile)

# Set expand distance for plotting - 10km
expand <- 10000

# Briefly summarise the data
summariseSimulationOutput(simOutput)

# Plot the temporal sampling range
plotTemporalRange(simOutput)

# Check that the dates in the simulation output make sense
checkSimulationOutputDates(simOutput)

# Report how many seeds were used in simulation
countSeeds(simOutput)

# Plot the infection dynamics over time
simulationDynamics <- plotInfectionDynamics(simOutput)

# Plot the spatial locations of the herds and setts in simulation
plotSimulationLocations(herdCoords, settCoords, badgerCentre, expand=expand, inner=inner)

# Plot the spatial locations of the recorded herds and setts
plotSpatialLocations(simOutput, badgerCentre, main="Recorded herd and sett locations", expand=expand, inner=inner)

# Plot the spatial locations of the sampled herds and setts
sampled <- simOutput[is.na(simOutput$DetectionDate) == FALSE, ]
plotSpatialLocations(sampled, badgerCentre, main="Sampled herd and sett locations", expand=expand, inner=inner)

#### Build the sequences ####

# Add a nucleotide sequence, based on the SNPs, for each individual
simOutput <- generateSequences(simOutput)

#### Build BASTA xml file ####

# Select the sequences
selected <- randomlySelectBadgerAndCattleSamples(simOutput, nBadgers=50, nCattle=50, inner=inner)
plotSpatialLocations(selected, badgerCentre, main="Selected herd and sett locations", expand=expand, inner=inner)
plotTemporalRange(selected, main="Temporal range of selected samples")

# Set parameters for the BASTA analysis
equalOrVaryingPopSizes <- "varying"
relaxedOrStrict <- "strict"
demeStructures <- c("2Deme", "4Deme")
chainLength <- 300000000

# Examine each deme structure
for(demeStructure in demeStructures){
  
  # Assign the isolates to Demes
  selected <- assignIsolatesToDemes(demeStructure, selected=selected, innerThreshold=inner,
                                    badgerCentre=badgerCentre, verbose=TRUE)
  
  # Build the XML file
  buildXMLFile(demeStructure, equalOrVaryingPopSizes, path, date, selected, chainLength, relaxedOrStrict, 
               sampleFromPrior=FALSE, estimateKappa=FALSE, addingConstantSites=FALSE)
  
}

# Close the pdf
dev.off()

#############
# FUNCTIONS #
#############

# Building the BASTA xml file - most methods directly from PreparingBASTAAlignment_24-01-18.R
# NOTE: there are minor changes - not estimating KAPPA or using constant sites

buildXMLFile <- function(demeStructure, equalOrVaryingPopSizes="varying", path, date, selected, chainLength,
                         relaxedOrStrict="strict", sampleFromPrior=FALSE, estimateKappa=FALSE, 
                         addingConstantSites=FALSE){
  
  # Note the output XML name
  outputFileName <- paste(getDemeInfo(demeStructure, "Name"), "_",
                          equalOrVaryingPopSizes, "_", relaxedOrStrict,
                          "_", date, sep="")
  if(sampleFromPrior == TRUE){
    outputFileName <- paste(getDemeInfo(demeStructure, "Name"), "_",
                            equalOrVaryingPopSizes, "_", relaxedOrStrict,
                            "_PRIOR_", date, sep="")
  }
  
  # Create a directory for the output file
  dir.create(paste(path, "/", outputFileName, sep=""))
  
  # Create an array of file lines and add in each necessary block
  fileLines <- startBuildingOutputFileLines()
  
  # Add Sequenced block
  fileLines <- addSequenceBlock(fileLines, selected, sampleFromPrior, addingConstantSites)
  
  # Add distribution block - if relaxed
  if(relaxedOrStrict == "relaxed"){
    fileLines <- addDistributionBlock(fileLines)
  }
  
  # NO constant sites block added
  
  # Add sampling dates block
  fileLines <- addTipDateBlock(fileLines, selected)
  
  # Add deme assignment block
  fileLines <- addDemeAssignmentBlock(fileLines, selected)
  
  # Add substitution model block
  fileLines <- addSubstitutionModelBlock(fileLines, selected, relaxedOrStrict,
                                         sampleFromPrior, estimateKappa)
  
  # Add branch rates model - if relaxed
  if(relaxedOrStrict == "relaxed"){
    fileLines <- addBranchRateModelBlock(fileLines)
  }
  
  # Add migration model block
  fileLines <- addMigrationModelBlock(fileLines, demeStructure, initialValue=0.1)
  
  # Add Prior distributions block
  fileLines <- addPriorsBlock(fileLines, demeStructure, relaxedOrStrict, estimateKappa)
  
  # Add Tree likelihood block
  fileLines <- addTreeLikelihoodBlock(fileLines, relaxedOrStrict)
  
  # Add Migration rate likelihood block
  fileLines <- addMigrationRateLikelihoodBlock(fileLines)
  
  # Add Beast settings block
  fileLines <- addBeastSettingsBlock(fileLines, demeStructure,
                                     equalOrVaryingPopSizes == "varying", outputFileName,
                                     initialValue = 0.1, chainLength, relaxedOrStrict)
  
  # Add the end to the XML file
  fileLines <- addEnd(fileLines)
  
  # Write the file out
  writeXMLFile(fileLines, paste(path, outputFileName, "/", outputFileName, ".xml", sep=""))
  
}

assignIsolatesToDemes <- function(demeStructure, selected, innerThreshold,
                                  badgerCentre, verbose=FALSE){
  
  # Initialise a column to store each isolate's deme assignment
  selected$Deme <- NA
  
  # Examine each isolate
  for(row in 1:nrow(selected)){
    
    # Check if location data available and note whether inner/outer, east/west, north/south
    inOrOut <- "inner"
    eastOrWest <- sample(c("east", "west"), size=1)
    northOrSouth <- sample(c("north", "south"), size=1)
    if(is.na(selected[row, "X"]) == FALSE){
      inOrOut <- checkIfInner(selected[row, "Distance"], innerThreshold)
      eastOrWest <- checkIfEast(selected[row, "X"], badgerCentre[1])
      northOrSouth <- checkIfNorth(selected[row, "Y"], badgerCentre[2])
    }else if(verbose == TRUE){
      
      cat(paste("Location data not available for isolate:", selected[row, "IsolateID"], "\n"))
    }
    
    # Note species
    badgerOrCow <- "badger"
    if(selected[row, "isCow"]){
      badgerOrCow <- "cow"
    }
    
    # Build deme assignment
    selected[row, "Deme"] <- buildDemeAssignment(inOrOut, eastOrWest,
                                                    northOrSouth, badgerOrCow,
                                                    demeStructure)
  }
  
  return(selected)
}

checkIfNorth <- function(y, badgerCentreY){
  
  output <- "north"
  if(y < badgerCentreY){
    output <- "south"
  }
  
  return(output)
}

checkIfEast <- function(x, badgerCentreX){
  
  output <- "east"
  if(x < badgerCentreX){
    output <- "west"
  }
  
  return(output)
}

checkIfInner <- function(distance, innerThreshold){
  
  output <- "inner"
  if(distance > innerThreshold){
    output <- "outer"
  }
  
  return(output)
}

buildDemeAssignment <- function(inOrOut, eastOrWest, northOrSouth, badgerOrCow, demeStructure){
  
  # Initialise the assignment
  output <- NULL
  
  # Build the assignment with the input information - specific for the deme structure
  if(demeStructure == "2Deme"){
    output <- badgerOrCow
    
  }else if(demeStructure == "3Deme-outerIsBoth"){
    
    if(inOrOut == "inner"){
      output <- paste(badgerOrCow, inOrOut, sep="-")
    }else{
      output <- "outer"
    }
    
  }else if(demeStructure == "3Deme-outerIsBadger"){
    
    if(badgerOrCow == "badger"){
      output <- paste(badgerOrCow, inOrOut, sep="-")
    }else{
      output <- badgerOrCow
    }
    
  }else if(demeStructure == "3Deme-outerIsCattle"){
    
    if(badgerOrCow == "badger"){
      output <- badgerOrCow
    }else{
      output <- paste(badgerOrCow, inOrOut, sep="-")
    }
    
  }else if(demeStructure == "4Deme"){
    
    output <- paste(badgerOrCow, inOrOut, sep="-")
    
  }else if(demeStructure == "6Deme-EastWest"){
    
    if(inOrOut == "inner"){
      output <- paste(badgerOrCow, inOrOut, sep="-")
    }else{
      output <- paste(badgerOrCow, inOrOut, eastOrWest, sep="-")
    }
    
  }else if(demeStructure == "6Deme-NorthSouth"){
    
    if(inOrOut == "inner"){
      output <- paste(badgerOrCow, inOrOut, sep="-")
    }else{
      output <- paste(badgerOrCow, inOrOut, northOrSouth, sep="-")
    }
    
  }else if(demeStructure == "8Deme-EastWest"){
    
    output <- paste(badgerOrCow, inOrOut, eastOrWest, sep="-")
    
  }else if(demeStructure == "8Deme-NorthSouth"){
    
    output <- paste(badgerOrCow, inOrOut, northOrSouth, sep="-")
    
  }else{
    cat(paste("Deme structure name provided isn't recognised:", demeStructure, "\n"))
  }
  
  return(output)
}

getDemeInfo <- function(deme, infoWanted){
  
  # Define the information for each Deme
  demeInfo <- list(
    "2Deme" = list("Name"="2Deme", "NumberDemes"=2, 
                   "MigrationRateMatrix"=matrix(c(NA, 1,
                                                  1, NA), byrow=TRUE, nrow=2)),
    "3Deme-outerIsBoth" = list("Name"="3Deme-outerIsBoth", "NumberDemes"=3, 
                               "MigrationRateMatrix"=matrix(c(NA, 1, 1,
                                                              1, NA, 1,
                                                              1, 1, NA), 
                                                            byrow=TRUE, nrow=3)),
    "3Deme-outerIsBadger" = list("Name"="3Deme-outerIsBadger", "NumberDemes"=3, 
                                 "MigrationRateMatrix"=matrix(c(NA, 1, 1,
                                                                1, NA, 1,
                                                                1, 1, NA), 
                                                              byrow=TRUE, nrow=3)),
    "3Deme-outerIsCattle" = list("Name"="3Deme-outerIsCattle", "NumberDemes"=3, 
                                 "MigrationRateMatrix"=matrix(c(NA, 1, 0,
                                                                1, NA, 1,
                                                                0, 1, NA), 
                                                              byrow=TRUE, nrow=3)),
    "4Deme" = list("Name"="4Deme", "NumberDemes"=4, 
                   "MigrationRateMatrix"=matrix(c(NA, 1, 0, 1,
                                                  1, NA, 1, 0,
                                                  0, 1, NA, 1,
                                                  1, 0, 1, NA),
                                                byrow=TRUE, nrow=4)),
    "6Deme-EastWest" = list("Name"="6Deme-EastWest", "NumberDemes"=6, 
                            "MigrationRateMatrix"=matrix(c(NA, 1, 0, 0, 1, 1,
                                                           1, NA, 1, 1, 0, 0,
                                                           0, 1, NA, 1, 1, 0,
                                                           0, 1, 1, NA, 0, 1,
                                                           1, 0, 1, 0, NA, 1,
                                                           1, 0, 0, 1, 1, NA),
                                                         byrow=TRUE, nrow=6)),
    "6Deme-NorthSouth" = list("Name"="6Deme-NorthSouth", "NumberDemes"=6, 
                              "MigrationRateMatrix"=matrix(c(NA, 1, 0, 0, 1, 1,
                                                             1, NA, 1, 1, 0, 0,
                                                             0, 1, NA, 1, 1, 0,
                                                             0, 1, 1, NA, 0, 1,
                                                             1, 0, 1, 0, NA, 1,
                                                             1, 0, 0, 1, 1, NA),
                                                           byrow=TRUE, nrow=6)),
    "8Deme-EastWest" = list("Name"="8Deme-EastWest", "NumberDemes"=8, 
                            "MigrationRateMatrix"=matrix(c(NA, 1, 1, 0, 0, 0, 1, 0,
                                                           1, NA, 0, 1, 0, 0, 0, 1,
                                                           1, 0, NA, 1, 1, 0, 0, 0,
                                                           0, 1, 1, NA, 0, 1, 0, 0,
                                                           0, 0, 1, 0, NA, 1, 1, 0,
                                                           0, 0, 0, 1, 1, NA, 0, 1,
                                                           1, 0, 0, 0, 1, 0, NA, 1,
                                                           0, 1, 0, 0, 0, 1, 1, NA),
                                                         byrow=TRUE, nrow=8)),
    "8Deme-NorthSouth" = list("Name"="8Deme-NorthSouth", "NumberDemes"=8, 
                              "MigrationRateMatrix"=matrix(c(NA, 1, 1, 0, 0, 0, 1, 0,
                                                             1, NA, 0, 1, 0, 0, 0, 1,
                                                             1, 0, NA, 1, 1, 0, 0, 0,
                                                             0, 1, 1, NA, 0, 1, 0, 0,
                                                             0, 0, 1, 0, NA, 1, 1, 0,
                                                             0, 0, 0, 1, 1, NA, 0, 1,
                                                             1, 0, 0, 0, 1, 0, NA, 1,
                                                             0, 1, 0, 0, 0, 1, 1, NA),
                                                           byrow=TRUE, nrow=8))
  )
  
  return(demeInfo[[deme]][[infoWanted]])
}

startBuildingOutputFileLines <- function(){
  
  # Create a vector to store the fileLines
  fileLines <- c()
  
  # Insert the first two lines
  fileLines[1] <- "<beast version='2.0' namespace='beast.evolution.alignment:beast.core:beast.core.parameter:beast.evolution.tree:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood:beast.evolution.tree:beast.math.distributions:multitypetreeVolz.distributions:multitypetreeVolz.operators:multitypetreeVolz.util'>"
  fileLines[2] <- ""
  
  return(fileLines)
}

addSequenceBlock <- function(fileLines, selected, sampleFromPrior, addingConstantSites=FALSE){
  
  fileLines[length(fileLines) + 1] <- "\t<!-- Sequence Alignment -->"
  if(sampleFromPrior == FALSE && addingConstantSites == TRUE){
    fileLines[length(fileLines) + 1] <- "\t<data id=\"alignmentVar\" dataType=\"nucleotide\">"
  }else{
    fileLines[length(fileLines) + 1] <- "\t<data id=\"alignment\" dataType=\"nucleotide\">"
  }
  
  for(row in 1:nrow(selected)){
    
    name <- paste(selected[row, "animal_ID"], selected[row, "Deme"], sep="_")
    fileLines[length(fileLines) + 1] <- paste("\t\t<sequence taxon=\"", name, "\">", sep="")
    fileLines[length(fileLines) + 1] <- paste("\t\t\t", selected[row, "Sequence"], sep="")
    fileLines[length(fileLines) + 1] <- "\t\t</sequence>"
  }
  fileLines[length(fileLines) + 1] <- "\t</data>"
  
  return(fileLines)
}

addDistributionBlock <- function(fileLines){
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- Initialise distributions to refer to -->"
  fileLines[length(fileLines) + 1] <- "\t<map name=\"Exponential\" >beast.math.distributions.Exponential</map>"
  
  return(fileLines)
}

addTipDateBlock <- function(fileLines, selected){
  
  # Create a decimal date column
  decimalDates <- decimal_date(selected$DateSampleTaken)
  
  # Print out tip dates
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- Tip Dates -->"
  fileLines[length(fileLines) + 1] <- "\t<timeTraitSet spec='beast.evolution.tree.TraitSet' id='timeTraitSet' traitname=\"date-forward\""
  output <- "\t\tvalue=\""
  for(row in 1:nrow(selected)){
    name <- paste(selected[row, "animal_ID"], selected[row, "Deme"], sep="_")
    
    output <- paste(output, name, "=", decimalDates[row], sep="")
    
    if(row < nrow(selected)){
      output <- paste(output, ",", sep="")
    }
  }
  output <- paste(output, "\">", sep="")
  fileLines[length(fileLines) + 1] <- output
  fileLines[length(fileLines) + 1] <- "\t\t<taxa spec='TaxonSet' alignment='@alignment'/>"
  fileLines[length(fileLines) + 1] <- "\t</timeTraitSet>"
  
  return(fileLines)
}

addDemeAssignmentBlock <- function(fileLines, selected){
  
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- Deme Assignment -->"
  output <- "\t<typeTraitSet id=\"typeTraitSet\" spec=\"TraitSet\" traitname=\"type\" value=\""
  for(row in 1:nrow(selected)){
    
    name <- paste(selected[row, "animal_ID"], selected[row, "Deme"], sep="_")
    
    output <- paste(output, name, "=", selected[row, "Deme"], sep="")
    
    if(row < nrow(selected)){
      output <- paste(output, ",", sep="")
    }
  }
  output <- paste(output, "\">", sep="")
  fileLines[length(fileLines) + 1] <- output
  fileLines[length(fileLines) + 1] <- "\t\t<taxa spec='TaxonSet' alignment='@alignment'/>"
  fileLines[length(fileLines) + 1] <- "\t</typeTraitSet>"
  
  return(fileLines)
}

addSubstitutionModelBlock <- function(fileLines, selected, relaxedOrStrict,
                                      sampleFromPrior, estimateKappa=FALSE){
  
  # Calculate the average site frequencies
  baseFrequencies <- paste(c(0.25, 0.25, 0.25, 0.25), collapse=" ")
  if(sampleFromPrior == FALSE){
    baseFrequencies <- calculateAverageBaseFrequencies(selected$Sequence)
    baseFrequencies <- baseFrequencies / sum(baseFrequencies) # Sum to 1
    baseFrequencies <- paste(baseFrequencies, collapse=" ")
  }
  
  # Add the file lines for the block
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- Substitution model (HKY) -->"
  fileLines[length(fileLines) + 1] <- "\t<siteModel spec=\"SiteModel\" id=\"siteModel\">"
  if(relaxedOrStrict == "strict"){
    fileLines[length(fileLines) + 1] <- "\t\t<!-- Strong rate prior -->"
    fileLines[length(fileLines) + 1] <- "\t\t<mutationRate spec='RealParameter' id=\"mutationRate\" value=\"0.00000005\"  upper=\"0.0000003\" lower=\"0.00000001\"/>"
  }
  fileLines[length(fileLines) + 1] <- "\t\t<substModel spec=\"HKY\">"
  if(estimateKappa == TRUE){
    fileLines[length(fileLines) + 1] <- "\t\t\t<kappa spec='RealParameter' id=\"hky.kappa\" value=\"1.0\"/>"
  }else{
    fileLines[length(fileLines) + 1] <- "\t\t\t<kappa estimate=\"false\" spec='RealParameter' id=\"hky.kappa\" value=\"1.0\"/>"
  }
  fileLines[length(fileLines) + 1] <- "\t\t\t<frequencies estimate=\"false\" spec='Frequencies'>"
  fileLines[length(fileLines) + 1] <- paste("\t\t\t\t<frequencies spec='RealParameter' id=\"hky.freq\" value=\"", 
                                            baseFrequencies, "\"/>", sep="")
  fileLines[length(fileLines) + 1] <- "\t\t\t</frequencies>"
  fileLines[length(fileLines) + 1] <- "\t\t</substModel>"
  fileLines[length(fileLines) + 1] <- "\t</siteModel>"
  
  return(fileLines)
}

calculateAverageBaseFrequencies <- function(sequences){
  
  frequencies <- c(0,0,0,0)
  for(sequence in sequences){
    
    nucleotides <- strsplit(sequence, split="")[[1]]
    
    counts <- as.vector(table(nucleotides)[c("A", "C", "G", "T")])
    
    frequencies <- frequencies + (counts / length(nucleotides))
  }
  
  frequencies <- frequencies / length(sequences)
  
  return(frequencies)
}

addBranchRateModelBlock <- function(fileLines){
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- Branch rate model -->"
  fileLines[length(fileLines) + 1] <- "\t<branchRateModel id=\"ExponentialRelaxedClock\" spec=\"beast.evolution.branchratemodel.UCRelaxedClockModel\" numberOfDiscreteRates=\"10\" tree=\"@tree\">"
  fileLines[length(fileLines) + 1] <- "\t\t<parameter name='clock.rate' id='ucedMean' value=\"0.00000005\"  upper=\"0.0000003\" lower=\"0.00000001\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t<Exponential id=\"Exponential\" name=\"distr\">"
  fileLines[length(fileLines) + 1] <- "\t\t\t<parameter id=\"UCExpLambda\" name=\"mean\">1.0</parameter>"
  fileLines[length(fileLines) + 1] <- "\t\t</Exponential>"
  fileLines[length(fileLines) + 1] <- "\t\t<rateCategories spec='parameter.IntegerParameter' id='expRateCategories' value=\"1\" dimension='10' estimate='true'/>"
  fileLines[length(fileLines) + 1] <- "\t</branchRateModel>"
  
  return(fileLines)
}

addMigrationModelBlock <- function(fileLines, demeStructure, initialValue){
  
  # Get the migration rate values
  migrationRateValues <- getMigrationModelValues(demeStructure, initialValue)
  migrationRateValuesString <- paste(migrationRateValues, collapse=" ")
  rateFlags <- paste(tolower(migrationRateValues != 0), collapse=" ")
  
  nDemes <- getDemeInfo(demeStructure, "NumberDemes")
  
  # Add in the file lines
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- Migration model -->"
  fileLines[length(fileLines) + 1] <- paste("\t<migrationModelVolz spec='MigrationModelVolz' nTypes=\"",
                                            nDemes, 
                                            "\" id='migModel'>", sep="")
  fileLines[length(fileLines) + 1] <- paste("\t\t<rateMatrix spec='RealParameter' value=\"",
                                            migrationRateValuesString, 
                                            "\" dimension=\"", 
                                            length(migrationRateValues),
                                            "\" id=\"rateMatrix\"/>", sep="")
  
  fileLines[length(fileLines) + 1] <- "\t\t<!-- BSSVS -->"
  fileLines[length(fileLines) + 1] <- paste("\t\t<rateMatrixFlags spec='BooleanParameter' value=\"",
                                            rateFlags,
                                            "\" dimension=\"",
                                            length(migrationRateValues),
                                            "\" id=\"rateMatrixFlags\"/>", sep="")
  fileLines[length(fileLines) + 1] <- ""      
  fileLines[length(fileLines) + 1] <- "\t\t<!-- Fixed Population Size -->"
  fileLines[length(fileLines) + 1] <- paste("\t\t<popSizes spec='RealParameter' value=\"1.0\" dimension=\"",
                                            nDemes,
                                            "\" id=\"popSizes\"/>", sep="")
  fileLines[length(fileLines) + 1] <- "\t</migrationModelVolz>"
  
  return(fileLines)
}

getMigrationModelValues <- function(deme, initialValue){
  
  rateMatrix <- getDemeInfo(deme, "MigrationRateMatrix")
  
  modelValues <- c(rateMatrix)
  modelValues <- modelValues[is.na(modelValues) == FALSE] * initialValue
  
  return(modelValues)
}

addPriorsBlock <- function(fileLines, demeStructure, relaxedOrStrict, estimateKappa=FALSE){
  
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- Parameter priors -->"
  fileLines[length(fileLines) + 1] <- "\t<input spec='CompoundDistribution' id='parameterPriors'>"
  if(relaxedOrStrict == "strict"){
    fileLines[length(fileLines) + 1] <- "\t\t<distribution spec='beast.math.distributions.Prior' x=\"@mutationRate\">"
    fileLines[length(fileLines) + 1] <- "\t\t\t<distr spec='LogNormalDistributionModel' M=\"-15.0\" S=\"6.0\"/>"
    fileLines[length(fileLines) + 1] <- "\t\t</distribution>"
  }else{
    fileLines[length(fileLines) + 1] <- "\t\t<distribution spec='beast.math.distributions.Prior' id=\"UCMeanRatePrior\" x=\"@ucedMean\">"
    fileLines[length(fileLines) + 1] <- "\t\t\t<distr spec='Exponential' mean=\"0.0000001\"/>"
    fileLines[length(fileLines) + 1] <- "\t\t</distribution>"
  }
  if(estimateKappa == TRUE){
    fileLines[length(fileLines) + 1] <- ""        
    fileLines[length(fileLines) + 1] <- "\t\t<distribution spec='beast.math.distributions.Prior' x=\"@hky.kappa\">"
    fileLines[length(fileLines) + 1] <- "\t\t\t<distr spec='LogNormalDistributionModel' M=\"0.0\" S=\"4.0\"/>"
    fileLines[length(fileLines) + 1] <- "\t\t</distribution>"
  }
  fileLines[length(fileLines) + 1] <- ""            
  fileLines[length(fileLines) + 1] <- "\t\t<distribution spec='beast.math.distributions.Prior' x=\"@rateMatrix\">"
  fileLines[length(fileLines) + 1] <- "\t\t\t<distr spec='Exponential' mean=\"0.01\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t</distribution>"
  fileLines[length(fileLines) + 1] <- ""                  
  fileLines[length(fileLines) + 1] <- "\t\t<distribution spec='beast.math.distributions.Prior' x=\"@popSizes\">"
  fileLines[length(fileLines) + 1] <- "\t\t\t<distr spec=\"LogNormalDistributionModel\"  M=\"0.0\" S=\"1.0\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t</distribution>"
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t\t<!-- BSSVS -->"
  fileLines[length(fileLines) + 1] <- "\t\t<distribution spec='beast.math.distributions.Prior'>"
  fileLines[length(fileLines) + 1] <- "\t\t\t<x spec='Sum' arg=\"@rateMatrixFlags\"/>"
  fileLines[length(fileLines) + 1] <- paste("\t\t\t<distr spec='Poisson' lambda=\"",
                                            getDemeInfo(demeStructure, "NumberDemes"),
                                            "\"/>", sep="")
  fileLines[length(fileLines) + 1] <- "\t\t</distribution>"
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t</input>"
  
  return(fileLines)
}

addTreeLikelihoodBlock <- function(fileLines, relaxedOrStrict){
  
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- Probability of sequence data given tree  -->"
  fileLines[length(fileLines) + 1] <- "\t<input spec='TreeLikelihood' id=\"treeLikelihood1\">"
  fileLines[length(fileLines) + 1] <- "\t\t<data idref=\"alignment\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t<tree idref=\"tree\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t<siteModel idref='siteModel'/>"
  if(relaxedOrStrict == "relaxed"){
    fileLines[length(fileLines) + 1] <- "\t\t<branchRateModel idref='ExponentialRelaxedClock'/>"
  }
  fileLines[length(fileLines) + 1] <- "\t</input>"
  
  return(fileLines)
}

addMigrationRateLikelihoodBlock <- function(fileLines){
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- Probability of tree given migration rates and population sizes -->"
  fileLines[length(fileLines) + 1] <- "\t<input spec='StructuredCoalescentTreeDensityVolz' id='treePrior'>"
  fileLines[length(fileLines) + 1] <- "\t\t<multiTypeTreeVolz idref=\"tree\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t<migrationModelVolz idref=\"migModel\"/>"
  fileLines[length(fileLines) + 1] <- "\t</input>"
  
  return(fileLines)
}

addBeastSettingsBlock <- function(fileLines, demeStructure, varyingPopSizes, outputFileName,
                                  initialValue, chainLength, relaxedOrStrict){
  
  # Get the migration rate values
  migrationRateValues <- getMigrationModelValues(demeStructure, initialValue)
  migrationRateValuesAsString <- paste(migrationRateValues, collapse=" ")
  nValuesEstimated <- length(which(migrationRateValues != 0))
  
  # Note the number of demes
  nDemes <- getDemeInfo(demeStructure, "NumberDemes")
  
  # Note whether the population sizes are equal of varying
  scaleAll <- "True"
  degreesOfFreedom <- 1
  if(varyingPopSizes == TRUE){
    scaleAll <- "False"
    degreesOfFreedom <- nDemes
  }
  
  # Note sampling frequency
  sampleEvery <- chainLength / 10000
  
  # Stop R shortening numbers using exponential notation
  options(scipen=100)
  
  # Write the text for the block
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- BEAST run settings -->"
  fileLines[length(fileLines) + 1] <- paste("\t<run spec=\"MCMC\" id=\"mcmc\" chainLength=\"",
                                            chainLength, "\" storeEvery=\"",
                                            sampleEvery, "\">", sep="")
  fileLines[length(fileLines) + 1] <- ""    
  fileLines[length(fileLines) + 1] <- "\t\t<!-- initialize tree at random according to StCoal  -->"
  fileLines[length(fileLines) + 1] <- "\t\t<init spec='StructuredCoalescentMultiTypeTreeVolz' id='tree'>"
  fileLines[length(fileLines) + 1] <- paste("\t\t\t<migrationModelVolz spec='MigrationModelVolz' nTypes=\"",
                                            nDemes, "\">", sep="")
  fileLines[length(fileLines) + 1] <- paste("\t\t\t\t<rateMatrix spec='RealParameter' value=\"",
                                            migrationRateValuesAsString, "\" dimension=\"",
                                            length(migrationRateValues), "\"/>", sep="")
  fileLines[length(fileLines) + 1] <- paste("\t\t\t\t<popSizes spec='RealParameter' value=\"1.0\" dimension=\"", 
                                            nDemes, "\"/>", sep="")
  fileLines[length(fileLines) + 1] <- "\t\t\t</migrationModelVolz>"
  fileLines[length(fileLines) + 1] <- "\t\t\t<trait idref='typeTraitSet'/>"
  fileLines[length(fileLines) + 1] <- "\t\t\t<trait idref='timeTraitSet'/>"
  fileLines[length(fileLines) + 1] <- "\t\t</init>"
  fileLines[length(fileLines) + 1] <- ""                
  fileLines[length(fileLines) + 1] <- "\t\t<state>"
  fileLines[length(fileLines) + 1] <- "\t\t\t<stateNode idref=\"tree\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t\t<stateNode idref=\"rateMatrix\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t\t<stateNode idref=\"rateMatrixFlags\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t\t<stateNode idref=\"popSizes\"/>"
  if(relaxedOrStrict == "strict"){
    fileLines[length(fileLines) + 1] <- "\t\t\t<stateNode idref=\"mutationRate\"/>"
  }else{
    fileLines[length(fileLines) + 1] <- "\t\t\t<stateNode idref=\"ucedMean\"/>"
    fileLines[length(fileLines) + 1] <- "\t\t\t<stateNode idref=\"expRateCategories\"/>"
  }
  fileLines[length(fileLines) + 1] <- "\t\t\t<stateNode idref=\"hky.kappa\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t\t<stateNode idref=\"hky.freq\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t</state>"
  fileLines[length(fileLines) + 1] <- ""                              
  fileLines[length(fileLines) + 1] <- "\t\t<distribution spec='CompoundDistribution' id='posterior'>"
  fileLines[length(fileLines) + 1] <- "\t\t\t<distribution idref=\"treeLikelihood1\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t\t<distribution idref='treePrior'/>"
  fileLines[length(fileLines) + 1] <- "\t\t\t<distribution idref=\"parameterPriors\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t</distribution>"
  fileLines[length(fileLines) + 1] <- ""                                      
  fileLines[length(fileLines) + 1] <- "\t\t<!-- Parameter scaling operators: HERE YOU CAN CHANGE THE WEIGHTINGS OF PARAMETERS BASED UPON HOW WELL THEY@RE BEING ESTIMATED -->"
  fileLines[length(fileLines) + 1] <- "\t\t<operator spec='ScaleOperator' id='RateScaler'"
  fileLines[length(fileLines) + 1] <- paste("\t\t\tparameter=\"@rateMatrix\" degreesOfFreedom=\"", 
                                            nValuesEstimated, "\"", sep="")
  fileLines[length(fileLines) + 1] <- "\t\t\tscaleFactor=\"0.8\" weight=\"10\">"
  fileLines[length(fileLines) + 1] <- "\t\t</operator>"
  fileLines[length(fileLines) + 1] <- ""                                        
  fileLines[length(fileLines) + 1] <- "\t\t<operator spec=\"ScaleOperator\" id=\"PopSizeScaler\""
  fileLines[length(fileLines) + 1] <- paste("\t\t\tparameter=\"@popSizes\" scaleAll=\"",
                                            scaleAll, "\" degreesOfFreedom=\"",
                                            degreesOfFreedom, "\"", sep="")
  fileLines[length(fileLines) + 1] <- "\t\t\tscaleFactor=\"0.8\" weight=\"2\">"
  fileLines[length(fileLines) + 1] <- "\t\t</operator>"
  fileLines[length(fileLines) + 1] <- ""                                        
  if(relaxedOrStrict == "strict"){
    fileLines[length(fileLines) + 1] <- "\t\t<operator spec=\"ScaleOperator\" id=\"muRateScaler\""
    fileLines[length(fileLines) + 1] <- "\t\t\tparameter=\"@mutationRate\""
    fileLines[length(fileLines) + 1] <- "\t\t\tscaleFactor=\"0.8\" weight=\"2\">"
    fileLines[length(fileLines) + 1] <- "\t\t</operator>"
  }else{
    fileLines[length(fileLines) + 1] <- "\t\t<operator spec=\"ScaleOperator\" id=\"ucedMeanScaler\""
    fileLines[length(fileLines) + 1] <- "\t\t\tparameter=\"@ucedMean\""
    fileLines[length(fileLines) + 1] <- "\t\t\tscaleFactor=\"0.5\" weight=\"2\">"
    fileLines[length(fileLines) + 1] <- "\t\t</operator>"
    fileLines[length(fileLines) + 1] <- ""  
    fileLines[length(fileLines) + 1] <- "\t\t<operator spec=\"IntRandomWalkOperator\" id=\"ExpCategoriesRandomWalk\""
    fileLines[length(fileLines) + 1] <- "\t\t\tparameter=\"@expRateCategories\""
    fileLines[length(fileLines) + 1] <- "\t\t\twindowSize=\"1\" weight=\"5\">"
    fileLines[length(fileLines) + 1] <- "\t\t</operator>"
    fileLines[length(fileLines) + 1] <- ""  
    fileLines[length(fileLines) + 1] <- "\t\t<operator spec=\"SwapOperator\" id=\"ExpCategoriesSwapOperator\""
    fileLines[length(fileLines) + 1] <- "\t\t\tintparameter=\"@expRateCategories\""
    fileLines[length(fileLines) + 1] <- "\t\t\tweight=\"5\">"
    fileLines[length(fileLines) + 1] <- "\t\t</operator>"
    fileLines[length(fileLines) + 1] <- ""  
    fileLines[length(fileLines) + 1] <- "\t\t<operator spec=\"UniformOperator\" id=\"ExpCategoriesUniform\""
    fileLines[length(fileLines) + 1] <- "\t\t\tparameter=\"@expRateCategories\""
    fileLines[length(fileLines) + 1] <- "\t\t\tweight=\"5\">"
    fileLines[length(fileLines) + 1] <- "\t\t</operator>"
  }
  if(estimateKappa == TRUE){
    fileLines[length(fileLines) + 1] <- ""
    fileLines[length(fileLines) + 1] <- "\t\t<operator spec='ScaleOperator' id='kappaScaler'"
    fileLines[length(fileLines) + 1] <- "\t\t\tparameter=\"@hky.kappa\""
    fileLines[length(fileLines) + 1] <- "\t\t\tscaleFactor=\"0.8\" weight=\"0.05\">"
    fileLines[length(fileLines) + 1] <- "\t\t</operator>"
  }
  fileLines[length(fileLines) + 1] <- ""                                        
  fileLines[length(fileLines) + 1] <- "\t\t<operator spec=\"DeltaExchangeOperator\" id=\"freqExchanger\""
  fileLines[length(fileLines) + 1] <- "\t\t\tparameter=\"@hky.freq\""
  fileLines[length(fileLines) + 1] <- "\t\t\tdelta=\"0.01\" weight=\"0.1\">"
  fileLines[length(fileLines) + 1] <- "\t\t</operator>"
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t\t<!-- BSVS -->"
  fileLines[length(fileLines) + 1] <- "\t\t<operator spec='BitFlipOperator' id='bitFlipOperator'"
  fileLines[length(fileLines) + 1] <- "\t\t\tparameter='@rateMatrixFlags' weight=\"1\">"
  fileLines[length(fileLines) + 1] <- "\t\t</operator>"
  fileLines[length(fileLines) + 1] <- ""
  
  fileLines[length(fileLines) + 1] <- "\t\t<!-- Multi-type tree operators -->"
  fileLines[length(fileLines) + 1] <- "\t\t<operator spec='TypedSubtreeExchangeVolz' id='STX' weight=\"5\" multiTypeTree=\"@tree\" migrationModel=\"@migModel\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t<operator spec=\"TypedWilsonBaldingVolz\" id=\"TWB\" weight=\"5\" multiTypeTree=\"@tree\" migrationModel=\"@migModel\" alpha=\"0.2\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t<operator spec=\"MultiTypeUniformVolz\" id=\"MTU\" weight=\"5\" multiTypeTree=\"@tree\" includeRoot=\"true\" rootScaleFactor=\"0.9\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t<operator spec=\"MultiTypeTreeScaleVolz\" id=\"MTTS2\" weight=\"5\" multiTypeTree=\"@tree\" scaleFactor=\"0.98\" useOldTreeScaler=\"true\"/>"
  fileLines[length(fileLines) + 1] <- ""                                       
  fileLines[length(fileLines) + 1] <- "\t\t<!-- Posterior Log to File -->"
  fileLines[length(fileLines) + 1] <- paste("\t\t<logger logEvery=\"",
                                            sampleEvery, "\" fileName=\"", 
                                            outputFileName, ".log", "\">", sep="")
  fileLines[length(fileLines) + 1] <- "\t\t\t<model idref='posterior'/>"
  fileLines[length(fileLines) + 1] <- "\t\t\t<log idref=\"posterior\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t\t<log idref=\"treeLikelihood1\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t\t<log spec='TreeRootTypeLoggerVolz' structuredCoalescentTreeDensityVolz='@treePrior'/>"
  fileLines[length(fileLines) + 1] <- "\t\t\t<log spec='MigrationModelLoggerVolz' migrationModel='@migModel'/>"
  fileLines[length(fileLines) + 1] <- "\t\t\t<log spec='MigrationCountsLoggerVolz' density='@treePrior'/>"
  if(relaxedOrStrict == "strict"){
    fileLines[length(fileLines) + 1] <- "\t\t\t<log idref=\"mutationRate\"/>"
  }else{
    fileLines[length(fileLines) + 1] <- "\t\t\t<log idref=\"ucedMean\"/>"
    fileLines[length(fileLines) + 1] <- "\t\t\t<log id='rateStat' spec='beast.evolution.branchratemodel.RateStatistic' tree='@tree' branchratemodel='@ExponentialRelaxedClock'/>"
  }
  
  fileLines[length(fileLines) + 1] <- "\t\t\t<log idref=\"hky.kappa\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t\t<log spec='beast.evolution.tree.TreeHeightLogger' tree='@tree'/>"
  fileLines[length(fileLines) + 1] <- "\t\t</logger>"
  fileLines[length(fileLines) + 1] <- ""                                        
  fileLines[length(fileLines) + 1] <- "\t\t<!-- Posterior Tree Log to File -->"
  fileLines[length(fileLines) + 1] <- paste("\t\t<logger logEvery=\"",
                                            sampleEvery, "\" fileName=\"",
                                            outputFileName, ".trees", "\" mode=\"tree\">", sep="")
  fileLines[length(fileLines) + 1] <- "\t\t\t<log idref=\"treePrior\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t</logger>"
  fileLines[length(fileLines) + 1] <- ""                                        
  fileLines[length(fileLines) + 1] <- "\t\t<!-- Posterior Log to Screen -->"
  fileLines[length(fileLines) + 1] <- paste("\t\t<logger logEvery=\"",
                                            sampleEvery, "\">", sep="")
  fileLines[length(fileLines) + 1] <- "\t\t\t<model idref='posterior'/>"
  fileLines[length(fileLines) + 1] <- "\t\t\t<log idref=\"posterior\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t\t<log idref=\"treeLikelihood1\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t\t<log spec='TreeRootTypeLoggerVolz' structuredCoalescentTreeDensityVolz='@treePrior'/>"
  if(relaxedOrStrict == "strict"){
    fileLines[length(fileLines) + 1] <- "\t\t\t<log idref=\"mutationRate\"/>"
  }else{
    fileLines[length(fileLines) + 1] <- "\t\t\t<log idref=\"ucedMean\"/>"
  }
  fileLines[length(fileLines) + 1] <- "\t\t\t<log idref=\"hky.kappa\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t\t<ESS spec='ESS' name='log' arg=\"@rateMatrix\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t</logger>"     
  fileLines[length(fileLines) + 1] <- ""                                       
  fileLines[length(fileLines) + 1] <- "\t</run>"
  
  options(scipen=0)
  
  return(fileLines)
}

addEnd <- function(fileLines){
  fileLines[length(fileLines) + 1] <- ""                                       
  fileLines[length(fileLines) + 1] <- "</beast>"
  
  return(fileLines)
}

writeXMLFile <- function(fileLines, outputFileName){
  
  # Note the deme structure names
  demeStructureName <- getDemeInfo(demeStructure, "Name")
  
  # Open an output file
  fileConnection <- file(outputFileName)
  
  # Print out file lines
  writeLines(fileLines, fileConnection)
  
  # Close the output file
  close(fileConnection)
}



# Processing and examining data

randomlySelectBadgerAndCattleSamples <- function(simOutput, nBadgers, nCattle, propCattleInside=0.25, inner=3500){
  
  # Select only the sampled animals
  sampled <- simOutput[is.na(simOutput$DetectionDate) == FALSE, ]
  
  # Randomly sample nBadgers and nCattle
  # Only sample badgers from within threshold inner distance
  # Sample half cattle 25% inside and 75% outside
  selectedBadgerRows <- sample(which(is.na(simOutput$DetectionDate) == FALSE & 
                                     simOutput$isCow == FALSE & simOutput$Distance <= inner),
                               size=nBadgers, replace=FALSE)
  selectedCattleRows <- c(sample(which(is.na(simOutput$DetectionDate) == FALSE & 
                                       simOutput$isCow == TRUE & simOutput$Distance <= inner), 
                                 size=nCattle * propCattleInside, replace=FALSE),
                          sample(which(is.na(simOutput$DetectionDate) == FALSE & 
                                 simOutput$isCow == TRUE & simOutput$Distance > inner), 
                                 size=nCattle * (1-propCattleInside), replace=FALSE))
  
  return(simOutput[c(selectedBadgerRows, selectedCattleRows), ])
}

euclideanDistance <- function(x1, y1, x2, y2){
  return(sqrt((x1 - x2)^2 + (y1 - y2)^2))
}

plotSimulationLocations <- function(herdCoords, settCoords, badgerCentre, expand=10000, inner=3500){

  # Calculate the range of the X and Y axes
  xRange <- c(badgerCentre[1] - expand, badgerCentre[1] + expand)
  yRange <- c(badgerCentre[2] - expand, badgerCentre[2] + expand)
  
  # Create an empty plot
  plot(x=NULL, y=NULL, bty="n", asp=1, xlim=xRange, ylim=yRange, 
       xlab="", xaxt="n", ylab="", yaxt="n", 
       main="Simulated herd and sett locations")
  
  # Add central point
  points(x=badgerCentre[1], y=badgerCentre[2], pch=15, cex=2)
  
  # Add the herds
  points(x=herdCoords$X, y=herdCoords$Y, pch=17, col=rgb(0,0,1, 0.5))
  
  # Add the setts
  points(x=settCoords$X, y=settCoords$Y, pch=19, col=rgb(1,0,0, 0.5))
  
  # Add inner and outer circles
  draw.circle(x=badgerCentre[1], y=badgerCentre[2], radius=inner,
              border="black")
  draw.circle(x=badgerCentre[1], y=badgerCentre[2], radius=expand,
              border="black")
  text(x=badgerCentre[1], y=(badgerCentre[2] - inner) - 400,
       labels=paste("Inner: ", paste(inner / 1000, "km radius", sep="")))
  text(x=badgerCentre[1], y=(badgerCentre[2] - expand) - 400,
       labels=paste("Outer: ", paste(expand / 1000, "km radius", sep="")))
  
  # Add legend
  legend("bottomright", legend=c("Sett", "Herd"), pch=c(19, 17), col=c("red", "blue"),
         text.col=c("red", "blue"), bty="n")
}

addUnitCoords <- function(simOutput, herdCoords, settCoords, badgerCentre){
  
  # Note the coordinates of each herd and sett - convert to list
  herdUnitCoords <- noteCoordinates(herdCoords)
  settUnitCoords <- noteCoordinates(settCoords)
  
  # Add empty columns into the simulation output to store the unit coordinates
  simOutput$X <- NA
  simOutput$Y <- NA
  simOutput$Distance <- NA
  
  # Examine the simulation output
  for(row in seq_len(nrow(simOutput))){
    
    # Get the unit ID of the last unit animal inhabited
    unit <- strsplit(simOutput[row, "unit_ID"], split=";")[[1]]
    unit <- unit[length(unit)]
    
    # Check if cow and unit coordinates exist
    if(simOutput[row, "isCow"] && is.null(herdUnitCoords[[unit]]) == FALSE){
      simOutput[row, c("X", "Y")] <- herdUnitCoords[[unit]]
    }else if(is.null(settUnitCoords[[unit]]) == FALSE){
      simOutput[row, c("X", "Y")] <- settUnitCoords[[unit]]
    }else{
      cat(paste0("Coordinate data not found for unit: ", unit, " in row: ", row, " ID: ", simOutput[row, "animal_ID"], "\n"))
    }
    
    # Calculate the distance to the Badger Centre
    if(is.na(unit)){
      next
    }
    simOutput[row, "Distance"] <- euclideanDistance(simOutput[row, "X"], simOutput[row, "Y"],
                                                    badgerCentre[1], badgerCentre[2])
  }
  
  return(simOutput)
}

noteCoordinates <- function(coords){
  
  # Initialise a list to store the coordinates associated with each unit
  unitCoords <- list()
  
  # Examine the coords
  for(row in seq_len(nrow(coords))){
    
    unitCoords[[as.character(coords[row, "ID"])]] <- c(coords[row, "X"], coords[row, "Y"])
  }
  
  return(unitCoords)
}

countUnits <- function(simOutput){
  
  # Initialise two lists - one for the badgers and one for the cattle
  counts <- list("SettCounts"=list(), "HerdCounts"=list())
  
  # Examine the simulation output
  for(row in seq_len(nrow(simOutput))){
    
    # Skip animal if no location data available
    if(is.na(simOutput[row, "X"])){
      next
    }
    
    # Build a key using the X and Y coordinates
    key <- paste0(simOutput[row, "X"], ":", simOutput[row, "Y"])
    
    # Check if cow
    if(simOutput[row, "isCow"]){
      
      # Check if key exists
      if(is.null(counts[["HerdCounts"]][[key]])){
        counts[["HerdCounts"]][[key]] <- 1
      }else{
        counts[["HerdCounts"]][[key]] <- counts[["HerdCounts"]][[key]] + 1
      }
      
    }else{
      
      # Check if key exists
      if(is.null(counts[["SettCounts"]][[key]])){
        counts[["SettCounts"]][[key]] <- 1
      }else{
        counts[["SettCounts"]][[key]] <- counts[["SettCounts"]][[key]] + 1
      }
    }
  }
  
  return(counts)
}

getRangeOfCounts <- function(counts){
  
  # Initialise a vector to store the maximum counts for the setts and herds
  rangeOfCounts <- c(-Inf, Inf, -Inf, Inf)
  
  # Examine the sett counts
  for(key in names(counts$SettCounts)){
    
    # Check for new maximum
    if(counts$SettCounts[[key]] > rangeOfCounts[1]){
      rangeOfCounts[1] <- counts$SettCounts[[key]]
    }
    
    # Check for new minimum
    if(counts$SettCounts[[key]] < rangeOfCounts[2]){
      rangeOfCounts[2] <- counts$SettCounts[[key]]
    }
  }
  
  # Examine the herd counts
  for(key in names(counts$HerdCounts)){
    
    # Check for new maximum
    if(counts$HerdCounts[[key]] > rangeOfCounts[3]){
      rangeOfCounts[3] <- counts$HerdCounts[[key]]
    }
    
    # Check for new minimum
    if(counts$HerdCounts[[key]] < rangeOfCounts[4]){
      rangeOfCounts[4] <- counts$HerdCounts[[key]]
    }
  }
  
  return(rangeOfCounts)
}

plotSpatialLocations <- function(simOutput, badgerCentre, main, expand=10000, inner=3500){
  
  # Get and set the margins
  currentMargins <- par("mar")
  par(mar=c(1,1,4,1))
  
  # Count the numbers of animals associated with each unit ID
  counts <- countUnits(simOutput)
  
  # Get the range of counts for the setts and herds
  rangeOfCounts <- getRangeOfCounts(counts)
  
  # Calculate the range of the X and Y axes
  xRange <- c(badgerCentre[1] - expand, badgerCentre[1] + expand)
  yRange <- c(badgerCentre[2] - expand, badgerCentre[2] + expand)
  
  # Create an empty plot
  plot(x=NULL, y=NULL, bty="n", asp=1, xlim=xRange, ylim=yRange, 
       xlab="", xaxt="n", ylab="", yaxt="n", main=main)
  
  # Add central point
  points(x=badgerCentre[1], y=badgerCentre[2], pch=15, cex=2)
  
  # Add inner and outer circles
  draw.circle(x=badgerCentre[1], y=badgerCentre[2], radius=inner,
              border="black")
  draw.circle(x=badgerCentre[1], y=badgerCentre[2], radius=expand,
              border="black")
  text(x=badgerCentre[1], y=(badgerCentre[2] - inner) - 400,
       labels=paste("Inner: ", paste(inner / 1000, "km radius", sep="")))
  text(x=badgerCentre[1], y=(badgerCentre[2] - expand) - 400,
       labels=paste("Outer: ", paste(expand / 1000, "km radius", sep="")))
  
  # Add the herds
  for(key in names(counts$HerdCounts)){
    
    coords <- as.numeric(strsplit(key, split=":")[[1]])
    
    points(x=coords[1], y=coords[2], pch=17, col=rgb(0,0,1, 0.5),
           cex=(counts$HerdCounts[[key]] / max(rangeOfCounts)) * 10)
  }
  
  # Add the setts
  for(key in names(counts$SettCounts)){
    coords <- as.numeric(strsplit(key, split=":")[[1]])
    
    points(x=coords[1], y=coords[2], pch=19, col=rgb(1,0,0, 0.5),
           cex=(counts$SettCounts[[key]] / max(rangeOfCounts)) * 10)
  }
  
  # Get the axis limits
  axisLimits <- par("usr") 
  
  # Define positions for legend
  yLength <- axisLimits[4] - axisLimits[3]
  xLength <- axisLimits[2] - axisLimits[1]
  positions <- seq(from=axisLimits[3], to=axisLimits[3] + 0.8* yLength, 
                   by=(0.8*yLength)/5)
  counts <- seq(from=min(rangeOfCounts), to=max(rangeOfCounts), 
                by=(max(rangeOfCounts) - min(rangeOfCounts))/5)
  for(i in seq_along(positions)){
    
    yPosition <- positions[i]
    count <- counts[i]
    points(x=axisLimits[1] + 0.95*xLength, y=yPosition, pch=19, 
           cex=(count / max(rangeOfCounts)) * 10, col=rgb(1,0,0, 0.5), xpd=TRUE)
    text(x=axisLimits[1] + 0.95*xLength, y=yPosition, labels=round(count, digits=0))
    points(x=axisLimits[1] + 0.05*xLength, y=yPosition, pch=17, 
           cex=(count / max(rangeOfCounts)) * 10, col=rgb(0,0,1, 0.5), xpd=TRUE)
    text(x=axisLimits[1] + 0.05*xLength, y=yPosition, labels=round(count, digits=0))
  }
  text(x=c(axisLimits[1] + 0.05*xLength, axisLimits[1] + 0.95*xLength),
       y=c(axisLimits[3] + 0.95* yLength, axisLimits[3] + 0.95* yLength),
       labels=c("Number/herd", "Number/sett"),
       col=c("blue", "red"), xpd=TRUE)
  
  # Reset the margins
  par(mar=currentMargins)

}

calculateCentralLocation <- function(herdXs, herdYs, settXs, settYs){
  
  # Calculate the centre point
  centreCoords <- c(
    mean(c(herdXs, settXs)),
    mean(c(herdYs, settYs))
  )
  
  return(centreCoords)
}

plotInfectionDynamics <- function(simOutput){
  
  # Get an ordered list of all dates in simulation output
  dates <- sort(unique(c(simOutput$InfectionDate, simOutput$DateSampleTaken)))
  
  # Index the dates
  indexedDates <- indexArrayOfDates(dates)
  
  # Initialise a table to store the population dynamics
  simulationDynamics <- data.frame("Date" = dates, 
                                   "NumberBadgersInfected"=0, "NumberBadgersSampled"=0,
                                   "NumberCattleInfected"=0, "NumberCattleSampled"=0)
  
  # Examine the simulation output and count how many infected and sampled animals there are at any given date
  for(row in seq_len(nrow(simOutput))){
    
    # Progress information
    cat(paste0("\rReading row: ", row, " of ", nrow(simOutput)))
    
    # Skip the root
    if(simOutput[row, "animal_ID"] == "ROOT"){
      next
    }
    
    # Get index of infection date
    simulationDynamics <- countDaysInfectedAndSampled(simulationDynamics,
                                                      infectionDate=simOutput[row, "InfectionDate"],
                                                      samplingDate=simOutput[row, "DateSampleTaken"],
                                                      indexedDates, isCow=simOutput[row, "isCow"])
  }
  cat("\rFinished!\t\t\t\t\t\t\t\t\t\t\t\n")
  
  # Plot the dynamics
  plot(x=simulationDynamics$Date, y=simulationDynamics$NumberBadgersInfected, 
       ylim=c(0, max(simulationDynamics[, 2:5])), 
       type="l", col="red", main="Simulation dynamics", ylab="Number", xlab="Date",
       bty="n", las=1)
  points(x=simulationDynamics$Date, y=simulationDynamics$NumberBadgersSampled,
         type="l", lty=2, col="red")
  points(x=simulationDynamics$Date, y=simulationDynamics$NumberCattleInfected,
         type="l", col="blue")
  points(x=simulationDynamics$Date, y=simulationDynamics$NumberCattleSampled,
         type="l", lty=2, col="blue")
  
  # Add a legend
  legend("topleft", legend=c("N. badgers infected", "N. badgers sampled", 
                             "N. cattle infected", "N. cattle sampled"),
         lty=c(1,2,1,2), col=c("red", "red", "blue", "blue"), bty="n",
         text.col=c("red", "red", "blue", "blue"))
  
  return(simulationDynamics)
}

countDaysInfectedAndSampled <- function(simulationDynamics, infectionDate, samplingDate, indexedDates, isCow){
  
  # Get the index of the infection date
  start <- 1 # Seeds don't have infection date but are infected from day 1
  if(is.na(infectionDate) == FALSE){
    start <- indexedDates[[as.character(infectionDate)]]
  }
  
  # Get the index of the sampling date
  end <- nrow(simulationDynamics) # If no sampling date - infected until last date
  if(is.na(samplingDate) == FALSE){
    end <- indexedDates[[as.character(samplingDate)]]
  }
  
  # Increment the days when animal was infected
  if(isCow == TRUE){
    simulationDynamics[start:end, "NumberCattleInfected"] <- simulationDynamics[start:end, "NumberCattleInfected"] + 1
  }else{
    simulationDynamics[start:end, "NumberBadgersInfected"] <- simulationDynamics[start:end, "NumberBadgersInfected"] + 1
  }
  
  # Increment the days when animal was sampled
  if(end != nrow(simulationDynamics) && isCow == TRUE){
    simulationDynamics[end:nrow(simulationDynamics), "NumberCattleSampled"] <- 
      simulationDynamics[end:nrow(simulationDynamics), "NumberCattleSampled"] + 1
  }else if(end != nrow(simulationDynamics) && isCow == FALSE){
    simulationDynamics[end:nrow(simulationDynamics), "NumberBadgersSampled"] <- 
      simulationDynamics[end:nrow(simulationDynamics), "NumberBadgersSampled"] + 1
  }
  
  return(simulationDynamics)
}

indexArrayOfDates <- function(array){
  output <- list()
  for(i in 1:length(array)){
    output[[as.character(array[i])]] <- i
  }
  
  return(output)
}

countSeeds <- function(simOutput){
  
  nBadgerSeeds <- length(which(grepl(simOutput$animal_ID, pattern="Badger_seed") == TRUE))
  nCattleSeeds <- length(which(grepl(simOutput$animal_ID, pattern="Cow_seed") == TRUE))
  
  cat(paste("Number of badger seeds =", nBadgerSeeds, 
            "\nNumber of cattle seeds =", nCattleSeeds, "\n"))
}

checkSimulationOutputDates <- function(simOutput){
  
  # Compare sampling dates and last SNP dates - should be EQUAL
  plot(x=simOutput$DateSampleTaken, y=simOutput$LastSnpGeneration, 
       xlab="Date sample taken", ylab="Date last SNP occurred", las=1, bty="n",
       main="Comparing sampling and last SNP dates")
  points(x=simOutput$DateSampleTaken, y=simOutput$DateSampleTaken, lty=2, col="red", type="l")
  
  # Compare sampling and detection dates - should be EQUAL
  plot(x=simOutput$DateSampleTaken, y=simOutput$DetectionDate, 
       xlab="Date sample taken", ylab="Date infection detected", las=1, bty="n",
       main="Comparing detection and sampling dates")
  points(x=simOutput$DateSampleTaken, y=simOutput$DateSampleTaken, lty=2, col="red", type="l")
  
  # Compare sampling and detection dates - infection should be before sampling
  plot(x=simOutput$DateSampleTaken, y=simOutput$InfectionDate, 
       xlab="Date sample taken", ylab="Date infected", las=1, bty="n",
       main="Comparing infection and sampling dates")
  points(x=simOutput$DateSampleTaken, y=simOutput$DateSampleTaken, lty=2, col="red", type="l")
}

summariseSimulationOutput <- function(simOutput){
  
  # Subset the data
  badgers <- simOutput[simOutput$isCow == "FALSE", ]
  sampledBadgers <- badgers[is.na(badgers$DateSampleTaken) == FALSE, ]
  cattle <- simOutput[simOutput$isCow == "TRUE", ]
  sampledCattle <- cattle[is.na(cattle$DateSampleTaken) == FALSE, ]
  
  # Report number of sampled cattle and badgers
  cat(paste0("Number sampled badgers: ", nrow(sampledBadgers), " of ", nrow(badgers),
             "\nNumber sampled cattle: ", nrow(sampledCattle), " of ", nrow(cattle), "\n"))
  
  # Report the sampling date range
  cat(paste0("Sampled badgers available from: ", min(sampledBadgers$DateSampleTaken), " to ", max(sampledBadgers$DateSampleTaken),
             "\nSampled cattle available from: ", min(sampledCattle$DateSampleTaken), " to ", max(sampledCattle$DateSampleTaken), "\n"))
  
  # Count the number of animals in each infection status category
  statusCounts <- countAnimalStatuses(badgers, cattle, sampledBadgers, sampledCattle)
}

countAnimalStatuses <- function(badgers, cattle, sampledBadgers, sampledCattle, plot=TRUE){
  
  statusCounts <- data.frame("Exposed"=c(NA, NA, NA, NA),
                             "Infectious"=c(NA, NA, NA, NA),
                             "TestSensitive"=c(NA, NA, NA, NA),
                             row.names=c("Badgers", "Cattle", "SampledBadgers", "SampledCattle"))
  statusCounts["Badgers", "Infectious"] <- nrow(badgers)
  statusCounts["SampledBadgers", "Infectious"] <- nrow(sampledBadgers)
  cattleCounts <- table(cattle$InfectionStatus)
  statusCounts["Cattle", "Exposed"] <- cattleCounts["EXPOSED"]
  statusCounts["Cattle", "Infectious"] <- cattleCounts["INFECTIOUS"]
  statusCounts["Cattle", "TestSensitive"] <- cattleCounts["TESTSENSITIVE"]
  sampledCattleCounts <- table(sampledCattle$InfectionStatus)
  statusCounts["SampledCattle", "Exposed"] <- sampledCattleCounts["EXPOSED"]
  statusCounts["SampledCattle", "Infectious"] <- sampledCattleCounts["INFECTIOUS"]
  statusCounts["SampledCattle", "TestSensitive"] <- sampledCattleCounts["TESTSENSITIVE"]
  
  if(plot){
    # Plot the counts
    currentMar <- par("mar")
    par(mar=c(4.1, 10, 4.1, 2.1))
    barplot(c(statusCounts["Badger", "Infectious"], statusCounts["SampledBadger", "Infectious"],
              statusCounts["Cattle", "Exposed"], statusCounts["SampledCattle", "Exposed"],
              statusCounts["Cattle", "Infectious"], statusCounts["SampledCattle", "Infectious"],
              statusCounts["Cattle", "TestSensitive"], statusCounts["SampledCattle", "TestSensitive"]),
            names.arg=c("Badger - I", "Sampled Badger - I",
                        "Cattle - E", "Sampled Cattle - E",
                        "Cattle - I", "Sampled Cattle - I",
                        "Cattle - T", "Sampled Cattle - T"), 
            horiz=TRUE, las=1, xlab="Number recorded", main="Numbers of samples from infection categories",
            col=c("red", rgb(1,0,0, 0.5), "blue", rgb(0,0,1, 0.5), "blue", rgb(0,0,1, 0.5), "blue", rgb(0,0,1, 0.5)))
    par(mar=currentMar)
  }
  
  return(statusCounts)
}

plotTemporalRange <- function(simOutput, main="Temporal range of sampling"){
  
  # Create two histograms of the sampling dates for badgers and cattle
  histBadgers <- hist(as.Date(simOutput[simOutput$isCow == "FALSE", "DateSampleTaken"]), breaks=10, freq=TRUE, plot=FALSE)
  histCattle <- hist(as.Date(simOutput[simOutput$isCow == "TRUE", "DateSampleTaken"]), breaks=10, freq=TRUE, plot=FALSE)
  
  # Get the y axis limits
  yLim <- c(0, max(c(histBadgers$counts, histCattle$counts)))
  
  # Plot the histograms
  # Note I removed the Y axis as it doesn't plot dates properly
  plot(histBadgers, 
       col=rgb(1,0,0, 0.5), las=1, xlab="Year",
       ylim=yLim, main=main, xaxt="n")
  plot(histCattle, col=rgb(0,0,1, 0.5), add=TRUE)
  
  # Add the X axis back in
  Axis(as.Date(simOutput[simOutput$isCow == "FALSE", "DateSampleTaken"]), col="black", side=1)
  
  # Add a legend
  legend("topleft", legend=c("Badgers", "Cattle"), text.col=c("red", "blue"), bty="n")
}

generateSequences <- function(simOutput){
  
  # Get the maximum SNP ID - equal to the number SNPs that occured in the population
  maxSNPID <- getMaxSNPID(simOutput)
  cat(paste0(maxSNPID, " SNPs found in simulated data\n"))
  
  # Create a reference and alternate allele for each SNP
  snpInfo <- createSNPs(maxSNPID)
  
  # Create each individual's sequence
  simOutput$Sequence <- NA
  
  # Examine each individual and create its sequence
  for(row in seq_len(nrow(simOutput))){
    
    # Skip individuals with no SNPs
    if(is.na(simOutput[row, "SNPs_detected"])){
      next
    }
    
    # First assign the reference sequence to the current individual
    sequence <- snpInfo$Reference
    
    # Get the SNPs associated with the current individual
    snps <- as.numeric(strsplit(simOutput[row, "SNPs_detected"], split=";")[[1]])
    
    # Assign each of the current individual's SNPs to its sequence
    for(snp in snps){
      sequence[snp] <- snpInfo$Alternate[snp]
    }
    
    # Convert the nucleotide array to string and store it
    simOutput[row, "Sequence"] <- paste(sequence, collapse="")
  }
  
  return(simOutput)
}

getMaxSNPID <- function(simOutput){
  
  # Initialise a variable to record the maximum SNP ID - the number of SNPs in the population
  maxSNPID <- -Inf
  
  # Examine the simulation data for each animal and store its SNP IDs
  for(row in seq_len(nrow(simOutput))){
    
    # Get the maximum SNP ID from the current individual
    individualMaxSNPID <- max(as.numeric(strsplit(simOutput[row, "SNPs_detected"], split=";")[[1]]))
    
    if(is.na(individualMaxSNPID)){
      cat(paste("No SNPs found in row:", row, "with animal_ID:", simOutput[row, "animal_ID"], "\n"))
      next
    }
    
    # Check if found new maximum
    if(individualMaxSNPID > maxSNPID){
      maxSNPID <- individualMaxSNPID
    }
  }
  
  return(maxSNPID)
}

createSNPs <- function(maxSNPID){
  
  # Initialise a matrix to store the reference alleles and alternates associated with each snp
  snpInfo <- list("Reference"=c(), "Alternate"=c())
  
  # Create the reference allele for each event
  nucleotides <- c("A", "C", "G", "T")
  snpInfo[["Reference"]] <- sample(nucleotides, size=maxSNPID, replace=TRUE)
  
  # Create the alternate allele for each event
  for(i in seq_len(maxSNPID)){
    snpInfo[["Alternate"]][i] <- sample(nucleotides[nucleotides != snpInfo[["Reference"]][i]], size=1)
  }
  
  return(snpInfo)
}