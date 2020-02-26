#### Preparation ####

# Load the required libraries
library(ape) # Reading in FASTA
library(lubridate) # Converting to decimal dates

# Set the path
path <- file.path("~", "Desktop", "BuildingBASTAXML")

# Set the date
date <- format(Sys.Date(), "%d-%m-%y")

#### Read in FASTA, metadata and constant site counts ####

# Read in the FASTA file
fastaFile <- file.path(path, "example.fasta")
sequences <- read.dna(fastaFile, format="fasta", as.character=TRUE)

# Read in the metadata
metadataFile <- file.path(path, "metadata.csv")
metadata <- read.table(metadataFile, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Convert sampling dates to decimal dates
metadata$DecimalDate <- decimal_date(as.Date(metadata$Date))

# Read in the constant site counts                          A,C,G,T proportions: 0.175, 0.325, 0.325, 0.175
# - I've estimated these: (genomeSize - nSites) split between A, C, G, T with 65% bias to G and C
constantSiteCounts <- generateConstantSiteCounts(nSites=ncol(sequences))

#### Build BASTA xml(s) ####

# Set the options to be used in the BASTA analyses
equalOrVaryingPopSizes <- "equal" # "equal" or "varying" ?
relaxedOrStrict <- "relaxed" # "relaxed" or "strict" ?
chainLength <- 300000000
sampleFromPrior <- FALSE

# Build the BASTA xml - I AM HERE!!!!!!!!!!!
buildXMLFile(sequences, metadata, equalOrVaryingPopSizes, path, date, 
             constantSiteCounts, chainLength, relaxedOrStrict, sampleFromPrior)

#### FUNCTIONS - preparation ####

generateConstantSiteCounts <- function(nSites, genomeSize=4345492, nucleotideProbs=c(0.175, 0.325, 0.325, 0.175)){
  
  # Calculate the number if sites not considered in the FASTA file
  nGenomeSitesNotConsidered <- genomeSize - nSites
  
  # Spread these not considered sites into nucleotide counts accordinag to therii respective probabilities
  nucleotideCounts <- ceiling(nucleotideProbs * nGenomeSitesNotConsidered)
  
  return(nucleotideCounts)
}

generateRandomMetadata <- function(nSequences, path){
  
  # Generate random data for sequences
  metadata <- data.frame("ID"=rownames(sequences), 
                         "Date"=sample(seq(from=as.Date("01-01-2005", format="%d-%m-%Y"),
                                           to=as.Date("01-12-2015", format="%d-%m-%Y"), by=1),
                                       size=nSequences, replace=TRUE),
                         "Species"=sample(c("badger", "cow"), size=nSequences, replace=TRUE))
  
  # Write the data to file
  metadataFile <- file.path(path, "metadata.csv")
  write.table(metadata, file=metadataFile, quote=FALSE, sep=",", row.names=FALSE)
}

#### FUNCTIONS - build BASTA xml ####

writeXMLFile <- function(fileLines, outputFileName){
  
  # Open an output file
  fileConnection <- file(outputFileName)
  
  # Print out file lines
  writeLines(fileLines, fileConnection)
  
  # Close the output file
  close(fileConnection)
}

addEnd <- function(fileLines){
  fileLines[length(fileLines) + 1] <- ""                                       
  fileLines[length(fileLines) + 1] <- "</beast>"
  
  return(fileLines)
}

addBeastSettingsBlock <- function(fileLines, varyingPopSizes, outputFileName,
                                  initialValue, chainLength, relaxedOrStrict, migrationRateMatrix,
                                  nDemes){
  
  # Get the migration rate values
  migrationRateValues <- getMigrationModelValues(migrationRateMatrix, initialValue)
  migrationRateValuesAsString <- paste(migrationRateValues, collapse=" ")
  nValuesEstimated <- length(which(migrationRateValues != 0))
  
  # Note whether the population sizes are equal of varying
  scaleAll <- "True"
  degreesOfFreedom <- 1
  if(varyingPopSizes == TRUE){
    scaleAll <- "False"
    degreesOfFreedom <- nDemes
  }
  
  # Note sampling frequency
  sampleEvery <- floor(chainLength / 10000)
  
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

addPriorsBlock <- function(fileLines, nDemes, relaxedOrStrict){
  
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
                                            nDemes, "\"/>", sep="")
  fileLines[length(fileLines) + 1] <- "\t\t</distribution>"
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t</input>"
  
  return(fileLines)
}

getMigrationModelValues <- function(rateMatrix, initialValue){
  
  modelValues <- c(rateMatrix)
  modelValues <- modelValues[is.na(modelValues) == FALSE] * initialValue
  
  return(modelValues)
}

addMigrationModelBlock <- function(fileLines, initialValue, migrationRateMatrix, nDemes){
  
  # Get the migration rate values
  migrationRateValues <- getMigrationModelValues(migrationRateMatrix, initialValue)
  migrationRateValuesString <- paste(migrationRateValues, collapse=" ")
  rateFlags <- paste(tolower(migrationRateValues != 0), collapse=" ")
  
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

calculateAverageBaseFrequencies <- function(sequences){
  
  # Initialise a vector to store the nucleotide frequencies
  frequencies <- c(0,0,0,0)
  
  # Examine each frequency
  for(row in 1:nrow(sequences)){
    
    # Calculate each number of times each nucleotide appears in current sequence
    counts <- as.vector(table(sequences[row, ])[c("a", "c", "g", "t")])
    
    # Covert counts to proportion of current sequence's length and add to growing proportion sum
    frequencies <- frequencies + (counts / ncol(sequences))
  }
  
  # Covert sum of proportions to mean
  frequencies <- frequencies / nrow(sequences)
  
  return(frequencies)
}

addSubstitutionModelBlock <- function(fileLines, sequences, relaxedOrStrict,
                                      sampleFromPrior){
  
  # Calculate the average site frequencies
  baseFrequencies <- paste(c(0.25, 0.25, 0.25, 0.25), collapse=" ")
  if(sampleFromPrior == FALSE){
    baseFrequencies <- calculateAverageBaseFrequencies(sequences)
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

addDemeAssignmentBlock <- function(fileLines, metadata){
  
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- Deme Assignment -->"
  output <- "\t<typeTraitSet id=\"typeTraitSet\" spec=\"TraitSet\" traitname=\"type\" value=\""
  for(row in 1:nrow(metadata)){
    
    name <- paste0(metadata[row, "ID"], "_", metadata[row, "Species"])
    
    output <- paste0(output, name, "=", metadata[row, "Species"])
    
    if(row < nrow(metadata)){
      output <- paste(output, ",", sep="")
    }
  }
  output <- paste(output, "\">", sep="")
  fileLines[length(fileLines) + 1] <- output
  fileLines[length(fileLines) + 1] <- "\t\t<taxa spec='TaxonSet' alignment='@alignment'/>"
  fileLines[length(fileLines) + 1] <- "\t</typeTraitSet>"
  
  return(fileLines)
}

addTipDateBlock <- function(fileLines, metadata){
  
  # Print out tip dates
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- Tip Dates -->"
  fileLines[length(fileLines) + 1] <- "\t<timeTraitSet spec='beast.evolution.tree.TraitSet' id='timeTraitSet' traitname=\"date-forward\""
  output <- "\t\tvalue=\""
  for(row in 1:nrow(metadata)){
    name <- paste0(metadata[row, "ID"], "_", metadata[row, "Species"])
    
    output <- paste0(output, name, "=", metadata[row, "DecimalDate"])
    
    if(row < nrow(metadata)){
      output <- paste(output, ",", sep="")
    }
  }
  output <- paste(output, "\">", sep="")
  fileLines[length(fileLines) + 1] <- output
  fileLines[length(fileLines) + 1] <- "\t\t<taxa spec='TaxonSet' alignment='@alignment'/>"
  fileLines[length(fileLines) + 1] <- "\t</timeTraitSet>"
  
  return(fileLines)
}

addConstantSiteCountsBlock <- function(fileLines, counts){
  
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

addSequenceBlock <- function(fileLines, sequences, metadata, sampleFromPrior){
  
  # Start the sequence block
  fileLines[length(fileLines) + 1] <- "\t<!-- Sequence Alignment -->"
  if(sampleFromPrior == FALSE){
    fileLines[length(fileLines) + 1] <- "\t<data id=\"alignmentVar\" dataType=\"nucleotide\">"
  }else{
    fileLines[length(fileLines) + 1] <- "\t<data id=\"alignment\" dataType=\"nucleotide\">"
  }
  
  # Add each of the sequences into the block
  for(row in 1:nrow(metadata)){
    
    # Build the sequence name
    name <- paste0(metadata[row, "ID"], "_", metadata[row, "Species"])
    
    # Build the sequence
    sequence <- paste(sequences[as.character(metadata[row, "ID"]), 1:ncol(sequences)], collapse="")
    
    # Build the file lines using the name and sequence
    fileLines[length(fileLines) + 1] <- paste0("\t\t<sequence taxon=\"", name, "\">")
    fileLines[length(fileLines) + 1] <- paste0("\t\t\t", ifelse(sampleFromPrior, "N", sequence))
    fileLines[length(fileLines) + 1] <- "\t\t</sequence>"
  }
  
  # Close the sequence block
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

buildXMLFile <- function(sequences, metadata, equalOrVaryingPopSizes, path, date, 
                         constantSiteCounts, chainLength,
                         relaxedOrStrict, sampleFromPrior){
  
  # Define the deme structure:
  # Demes:
  #   0 badger
  #   1 cow
  #
  #        | badger | cow | 
  # _______|________|_____|
  # badger |   0    |  0  |
  # _______|________|_____|
  # cow    |   0    |  0  |
  # _______|________|_____|
  migrationRateMatrix <- matrix(c(NA, 1,
                                  1, NA),
                                byrow=TRUE, nrow=2)
  
  # Note the number of demes
  nDemes <- 2
  
  # Note the output XML name
  outputFileName <- paste0("BASTA_", equalOrVaryingPopSizes, "_", relaxedOrStrict,
                          "_", date)
  if(sampleFromPrior == TRUE){
    outputFileName <- paste("BASTA_", equalOrVaryingPopSizes, "_", relaxedOrStrict,
                            "_PRIOR_", date)
  }
  
  # Create a directory for the output file
  dir.create(file.path(path, outputFileName), showWarnings=FALSE)
  
  # Create an array of file lines and add in each necessary block
  fileLines <- startBuildingOutputFileLines()
  
  # Add Sequenced block
  fileLines <- addSequenceBlock(fileLines, sequences, metadata, sampleFromPrior)
  
  # Add distribution block - if relaxed
  if(relaxedOrStrict == "relaxed"){
    fileLines <- addDistributionBlock(fileLines)
  }
  
  # Add constant site counts block
  if(sampleFromPrior == FALSE){
    fileLines <- addConstantSiteCountsBlock(fileLines, constantSiteCounts)
  }
  
  # Add sampling dates block
  fileLines <- addTipDateBlock(fileLines, metadata)
  
  # Add deme assignment block
  fileLines <- addDemeAssignmentBlock(fileLines, metadata)
  
  # Add substitution model block
  fileLines <- addSubstitutionModelBlock(fileLines, metadata, relaxedOrStrict,
                                         sampleFromPrior)
  
  # Add branch rates model - if relaxed
  if(relaxedOrStrict == "relaxed"){
    fileLines <- addBranchRateModelBlock(fileLines)
  }
  
  # Add migration model block
  fileLines <- addMigrationModelBlock(fileLines, initialValue=0.1, migrationRateMatrix, 
                                      nDemes=nDemes)
  
  # Add Prior distributions block
  fileLines <- addPriorsBlock(fileLines, nDemes=nDemes, relaxedOrStrict)
  
  # Add Tree likelihood block
  fileLines <- addTreeLikelihoodBlock(fileLines, relaxedOrStrict)
  
  # Add Migration rate likelihood block
  fileLines <- addMigrationRateLikelihoodBlock(fileLines)
  
  # Add Beast settings block
  fileLines <- addBeastSettingsBlock(fileLines, equalOrVaryingPopSizes == "varying", outputFileName,
                                     initialValue = 0.1, chainLength, relaxedOrStrict,
                                     migrationRateMatrix, nDemes=nDemes)
  
  # Add the end to the XML file
  fileLines <- addEnd(fileLines)
  
  # Write the file out
  writeXMLFile(fileLines, paste0(file.path(path, outputFileName, outputFileName), ".xml"))
}

