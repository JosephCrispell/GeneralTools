# Script to produce a FASTA alignment to be used in BASTA analysis
## Alignment from main WP clade
## Sampling dates of associated animals
## Sampling locations of associated animals
## Single isolates from the badgers involved (from within clade) - best quality

#############
# Libraries #
#############

library(ape)
library(geiger)
library(plotrix)
library(lubridate)

############################################
# Get a list of the isolates in main clade #
############################################

# Current date
date <- format(Sys.Date(), "%d-%m-%y")

# Set the path
path <- "/home/josephcrispell/Desktop/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/"

# Read in the newick tree
# NOTE!!! Isolate selection for BASTA clade already done in BuildAndPlotPhylogeny_08-02-18.R!!!
file <- paste(path, "vcfFiles/",  "mlTree_BASTAClade_DatedTips_27-03-18.tree", sep="")
#file <- paste(path, "vcfFiles/",  "mlTree_DatedTips_27-03-18.tree", sep="")
isolatesInClade <- getIsolatesFromTree(file)

###############################
# Get the isolate information #
###############################

# Note the sampling infor file names
cattleInfoFile <- paste(path, "IsolateData/",
              "CattleIsolateInfo_AddedNew_TB1484-TB1490_22-03-18.csv",
              sep="")
badgerInfoFile <- paste(path, "IsolateData/",
                        "BadgerInfo_08-04-15_LatLongs_XY_Centroids.csv", sep="")

# Get Isolate sampling information
selectedIsolateInfo <- getIsolateSamplingInformation(cattleInfoFile, badgerInfoFile, 
                                             isolates=names(isolatesInClade))

# Calculate distance of isolates to badger centre
badgerCentre <- c(381761.7, 200964.3)
selectedIsolateInfo <- calculateDistanceToBadgerCentre(badgerCentre, selectedIsolateInfo)

# Add the isolate sequences
file <- paste(path, "vcfFiles/", "sequences_withoutHomoplasies_27-03-18.fasta", sep="")
selectedIsolateInfo$Sequence <- getIsolateSequences(file, isolatesInClade)

##########################################################
# Build XML files according to different deme structures #
##########################################################

## Options

# Note whether sampling from prior
sampleFromPrior <- TRUE
selectedIsolateInfo <- removeSequencesIfSamplingFromPrior(selectedIsolateInfo, 
                                                          sampleFromPrior)

# Note deme structure to use
demeStructures <- c(
  "2Deme",
  "3Deme-outerIsBoth",
  "3Deme-outerIsBadger",
  "3Deme-outerIsCattle",
  "4Deme",
  "6Deme-EastWest",
  "6Deme-NorthSouth",
  "8Deme-EastWest",
  "8Deme-NorthSouth")

# Set options for population size estimate
popSizeEstimation <- c("equal", "varying") # "equal", "varying"

# Set options for the clock model
clockModels <- c("relaxed") # "strict", "relaxed"

# Options
chainLength <- 300000000

# Note the constant site counts file name
constantSiteCountsFile <- paste(path, "vcfFiles/", "constantSiteCounts_24-03-2018.txt", sep="")

# Define inner circle radius
innerDistance <- 3500

## Building XMLs

# Deme structure loop
for(demeStructure in demeStructures){
  
  # Assign the isolates to Demes
  selectedIsolateInfo <- assignIsolatesToDemes(demeStructure=demeStructure, selectedIsolateInfo,
                                               innerDistance, badgerCentre)
  
  # Population size estimation loop
  for(equalOrVaryingPopSizes in popSizeEstimation){
    
    # Clock model estimation loop
    for(relaxedOrStrict in clockModels){
      
      # Build the XML file
      buildXMLFile(demeStructure, equalOrVaryingPopSizes,
                   paste(path, "BASTA/", sep=""), date,
                   constantSiteCountsFile, selectedIsolateInfo, chainLength,
                   relaxedOrStrict, sampleFromPrior)
    }
  }
}


#############
# FUNCTIONS #
#############

removeSequencesIfSamplingFromPrior <- function(isolateInfo, sampleFromPrior){
  
  if(sampleFromPrior == TRUE){
    
    for(row in 1:nrow(isolateInfo)){
      isolateInfo[row, "Sequence"] <- "N"
    }
  }
  
  return(isolateInfo)
}

buildXMLFile <- function(demeStructure, equalOrVaryingPopSizes, path, date, 
                         constantSiteCountsFile, isolateInfo, chainLength,
                         relaxedOrStrict, sampleFromPrior){
  
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
  fileLines <- addSequenceBlock(fileLines, selectedIsolateInfo, sampleFromPrior)
  
  # Add distribution block - if relaxed
  if(relaxedOrStrict == "relaxed"){
    fileLines <- addDistributionBlock(fileLines)
  }
  
  # Add constant site counts block
  if(sampleFromPrior == FALSE){
    fileLines <- addConstantSiteCountsBlock(fileLines, constantSiteCountsFile)
  }
  
  # Add sampling dates block
  fileLines <- addTipDateBlock(fileLines, selectedIsolateInfo)
  
  # Add deme assignment block
  fileLines <- addDemeAssignmentBlock(fileLines, selectedIsolateInfo)
  
  # Add substitution model block
  fileLines <- addSubstitutionModelBlock(fileLines, isolateInfo, relaxedOrStrict,
                                         sampleFromPrior)
  
  # Add branch rates model - if relaxed
  if(relaxedOrStrict == "relaxed"){
    fileLines <- addBranchRateModelBlock(fileLines)
  }
  
  # Add migration model block
  fileLines <- addMigrationModelBlock(fileLines, demeStructure, initialValue=0.1)
  
  # Add Prior distributions block
  fileLines <- addPriorsBlock(fileLines, demeStructure, relaxedOrStrict)
  
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

addEnd <- function(fileLines){
  fileLines[length(fileLines) + 1] <- ""                                       
  fileLines[length(fileLines) + 1] <- "</beast>"
  
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
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t\t<operator spec='ScaleOperator' id='kappaScaler'"
  fileLines[length(fileLines) + 1] <- "\t\t\tparameter=\"@hky.kappa\""
  fileLines[length(fileLines) + 1] <- "\t\t\tscaleFactor=\"0.8\" weight=\"0.05\">"
  fileLines[length(fileLines) + 1] <- "\t\t</operator>"
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

addMigrationRateLikelihoodBlock <- function(fileLines){
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- Probability of tree given migration rates and population sizes -->"
  fileLines[length(fileLines) + 1] <- "\t<input spec='StructuredCoalescentTreeDensityVolz' id='treePrior'>"
  fileLines[length(fileLines) + 1] <- "\t\t<multiTypeTreeVolz idref=\"tree\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t<migrationModelVolz idref=\"migModel\"/>"
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

addPriorsBlock <- function(fileLines, demeStructure, relaxedOrStrict){
  
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
  fileLines[length(fileLines) + 1] <- ""        
  fileLines[length(fileLines) + 1] <- "\t\t<distribution spec='beast.math.distributions.Prior' x=\"@hky.kappa\">"
  fileLines[length(fileLines) + 1] <- "\t\t\t<distr spec='LogNormalDistributionModel' M=\"0.0\" S=\"4.0\"/>"
  fileLines[length(fileLines) + 1] <- "\t\t</distribution>"
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

addSubstitutionModelBlock <- function(fileLines, isolateInfo, relaxedOrStrict,
                                      sampleFromPrior){
  
  # Calculate the average site frequencies
  baseFrequencies <- paste(c(0.25, 0.25, 0.25, 0.25), collapse=" ")
  if(sampleFromPrior == FALSE){
    baseFrequencies <- calculateAverageBaseFrequencies(isolateInfo$Sequence)
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
  fileLines[length(fileLines) + 1] <- "\t\t\t<kappa spec='RealParameter' id=\"hky.kappa\" value=\"1.0\"/>"
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

addDemeAssignmentBlock <- function(fileLines, isolateInfo){
  
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- Deme Assignment -->"
  output <- "\t<typeTraitSet id=\"typeTraitSet\" spec=\"TraitSet\" traitname=\"type\" value=\""
  for(row in 1:nrow(isolateInfo)){
    
    name <- paste(isolateInfo[row, "IsolateID"], isolateInfo[row, "Deme"], sep="_")
    
    output <- paste(output, name, "=", isolateInfo[row, "Deme"], sep="")
    
    if(row < nrow(isolateInfo)){
      output <- paste(output, ",", sep="")
    }
  }
  output <- paste(output, "\">", sep="")
  fileLines[length(fileLines) + 1] <- output
  fileLines[length(fileLines) + 1] <- "\t\t<taxa spec='TaxonSet' alignment='@alignment'/>"
  fileLines[length(fileLines) + 1] <- "\t</typeTraitSet>"
  
  return(fileLines)
}

addTipDateBlock <- function(fileLines, isolateInfo){
  
  # Convert the sampling dates to date objects - 23/02/1999
  isolateInfo$SamplingDate <- as.Date(isolateInfo$SamplingDate,
                                               format="%d/%m/%Y")
  
  # Create a decimal date column
  decimalDates <- decimal_date(isolateInfo$SamplingDate)
  
  # Print out tip dates
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- Tip Dates -->"
  fileLines[length(fileLines) + 1] <- "\t<timeTraitSet spec='beast.evolution.tree.TraitSet' id='timeTraitSet' traitname=\"date-forward\""
  output <- "\t\tvalue=\""
  for(row in 1:nrow(isolateInfo)){
    name <- paste(isolateInfo[row, "IsolateID"], isolateInfo[row, "Deme"], sep="_")
    
    output <- paste(output, name, "=", decimalDates[row], sep="")
    
    if(row < nrow(isolateInfo)){
      output <- paste(output, ",", sep="")
    }
  }
  output <- paste(output, "\">", sep="")
  fileLines[length(fileLines) + 1] <- output
  fileLines[length(fileLines) + 1] <- "\t\t<taxa spec='TaxonSet' alignment='@alignment'/>"
  fileLines[length(fileLines) + 1] <- "\t</timeTraitSet>"
  
  return(fileLines)
}

addConstantSiteCountsBlock <- function(fileLines, constantSiteCountFile){
  
  ## Get constant site counts from file
  
  # Open the file and store the file lines
  connection <- file(constantSiteCountFile, "r")
  lines <- readLines(connection)
  close(connection)
  
  # Calculate the constant site counts
  counts <- c(0,0,0,0)
  for(line in lines){
    
    # Skip lines without counts
    if(grepl(line, pattern="Counts|A|Allele|counts") == TRUE || line == ""){
      next
    }
    
    # Split the line into its parts
    parts <- strsplit(line, "\t")[[1]]
    
    # Add the counts to the tally
    counts <- counts + as.numeric(parts)
  }
  
  ## Print out block
  counts <- paste(counts, collapse=" ")
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- Constant Site Counts -->"
  fileLines[length(fileLines) + 1] <- paste("\t<data id='alignment' spec='FilteredAlignment' filter='-' data='@alignmentVar' constantSiteWeights='",
                                            counts, "'/>", sep="")
  
  return(fileLines)
}

addDistributionBlock <- function(fileLines){
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- Initialise distributions to refer to -->"
  fileLines[length(fileLines) + 1] <- "\t<map name=\"Exponential\" >beast.math.distributions.Exponential</map>"
  
  return(fileLines)
}

addSequenceBlock <- function(fileLines, isolateInfo, sampleFromPrior){
  
  fileLines[length(fileLines) + 1] <- "\t<!-- Sequence Alignment -->"
  if(sampleFromPrior == FALSE){
    fileLines[length(fileLines) + 1] <- "\t<data id=\"alignmentVar\" dataType=\"nucleotide\">"
  }else{
    fileLines[length(fileLines) + 1] <- "\t<data id=\"alignment\" dataType=\"nucleotide\">"
  }
  
  
  for(row in 1:nrow(isolateInfo)){
    
    name <- paste(isolateInfo[row, "IsolateID"], isolateInfo[row, "Deme"], sep="_")
    fileLines[length(fileLines) + 1] <- paste("\t\t<sequence taxon=\"", name, "\">", sep="")
    fileLines[length(fileLines) + 1] <- paste("\t\t\t", isolateInfo[row, "Sequence"], sep="")
    fileLines[length(fileLines) + 1] <- "\t\t</sequence>"
  }
  fileLines[length(fileLines) + 1] <- "\t</data>"
  
  
  return(fileLines)
}

startBuildingOutputFileLines <- function(){
  
  # Create a vector to store the fileLines
  fileLines <- c()
  
  # Insert the first two lines
  fileLines[1] <- "<beast version='2.0' namespace='beast.evolution.alignment:beast.core:beast.core.parameter:beast.evolution.tree:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood:beast.evolution.tree:beast.math.distributions:multitypetreeVolz.distributions:multitypetreeVolz.operators:multitypetreeVolz.util'>"
  fileLines[2] <- ""
  
  return(fileLines)
}

assignIsolatesToDemes <- function(demeStructure, isolateInfo, innerThreshold,
                                  badgerCentre, verbose=FALSE){
  
  # Initialise a column to store each isolate's deme assignment
  isolateInfo$Deme <- rep("NA", nrow(isolateInfo))
  
  # Examine each isolate
  for(row in 1:nrow(isolateInfo)){
    
    # Check if location data available and note whether inner/outer, east/west, north/south
    inOrOut <- "inner"
    eastOrWest <- sample(c("east", "west"), size=1)
    northOrSouth <- sample(c("north", "south"), size=1)
    if(is.na(isolateInfo[row, "X"]) == FALSE){
      inOrOut <- checkIfInner(isolateInfo[row, "Distance"], innerThreshold)
      eastOrWest <- checkIfEast(isolateInfo[row, "X"], badgerCentre[1])
      northOrSouth <- checkIfNorth(isolateInfo[row, "Y"], badgerCentre[2])
    }else if(verbose == TRUE){
      
      cat(paste("Location data not available for isolate:", isolateInfo[row, "IsolateID"], "\n"))
    }
    
    # Note species
    badgerOrCow <- tolower(isolateInfo[row, "Species"])
    
    # Build deme assignment
    isolateInfo[row, "Deme"] <- buildDemeAssignment(inOrOut, eastOrWest,
                                                    northOrSouth, badgerOrCow,
                                                    demeStructure)
  }
  
  return(isolateInfo)
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

calculateDistanceToBadgerCentre <- function(badgerCentre, isolateInfo){
  
  isolateInfo$Distance <- rep(NA, nrow(isolateInfo))
  for(row in 1:nrow(isolateInfo)){
    
    # Skip isolates with an unknown location
    if(is.na(isolateInfo[row, "X"]) == FALSE){
      isolateInfo[row, "Distance"] <- euclideanDistance(x1=badgerCentre[1], y1=badgerCentre[2],
                                                         x2=isolateInfo[row, "X"],
                                                         y2=isolateInfo[row, "Y"])
    }  
  }
  
  return(isolateInfo)
}

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
    
    # Cattle
    if(grepl(pattern="TB|AF-|HI-", x=isolates[index]) == TRUE){
      
      # Find index in table
      strainIndex <- which(cattleInfo$StrainId == isolates[index])
      isolateInfo[index, "SamplingDate"] <- strsplit(as.character(cattleInfo[strainIndex, "BreakdownID"]),
                                                    split="-")[[1]][2] # 14082000501-23/02/1999
      isolateInfo[index, "X"] <- cattleInfo[strainIndex, "Mapx"]
      isolateInfo[index, "Y"] <- cattleInfo[strainIndex, "Mapy"]
      isolateInfo[index, "AnimalID"] <- cattleInfo[strainIndex, "Rawtag"]
      isolateInfo[index, "Species"] <- "COW"
      
      # Badgers
    }else if(grepl(pattern="WB", x=isolates[index]) == TRUE){
      
      # Find index in table
      strainIndex <- which(badgerInfo$WB_id == isolates[index])
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
  
  isolateInfo$SamplingDate <- as.Date(isolateInfo$SamplingDate, format="%d/%m/%Y")
  
  return(isolateInfo)
}

getIsolateSequences <- function(fastaFile, isolatesInClade){

  sequences <- readFASTA(fastaFile, skip=1)
  
  # Get the isolate IDs
  isolates <- names(isolatesInClade)
  
  # Get the sequences for the isolates
  isolateSequences <- c()
  for(i in 1:length(isolates)){
    
    # Check if sequence present for isolate
    if(is.null(sequences[[isolates[i]]]) == FALSE){
      
      isolateSequences[i] <- sequences[[isolates[i]]]
    }else{
      print(paste("Couldn't find sequence for: ", isolates[i]))
    }
  }
  
  return(isolateSequences)
}

getIsolatesFromTree <- function(treeFile){
  
  # Read in the BASTA clade tree
  tree <- read.tree(file=treeFile)

  # Get the tip labels - note that they'll have sampling dates attached to them
  isolates <- tree$tip.label
  for(i in 1:length(isolates)){
    isolates[i] = strsplit(isolates[i], split="_")[[1]][1]
  }

  # Convert this array to list
  isolatesInClade <- convertVectorToList(isolates)
  
  return(isolatesInClade)
}

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
