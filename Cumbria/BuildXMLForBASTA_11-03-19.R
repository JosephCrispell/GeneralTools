#### SET UP ####

# Set the path variable
path <- "/home/josephcrispell/Desktop/Research/Cumbria/"

#### Read in the sample information ####

# Load the sample information
samplingInfoFile <- paste(path, "17z_metadata_040319.csv", sep="")
samplingInfo <- read.table(samplingInfoFile, header=TRUE, stringsAsFactors=FALSE, sep=",")

#### Read in the FASTA file ####

# Read in the FASTA file
fastaFile <- paste(path, "vcfFiles/sequences_Prox-10_07-03-2019.fasta", sep="")
sequences <- readFASTA(fastaFile, skip=1)

# Remove the reference sequence
sequences[["Ref-1997"]] <- NULL

# Note the sampling info for each sequence
sequenceInfo <- getSamplingInfoForEachSequence(sequences, samplingInfo)

#### Prepare the XML file for BASTA ####

# Note the constant site counts file
constantSiteCountsFile <- paste(path, "vcfFiles/", "constantSiteCounts_07-03-2019.txt", sep="")

# Set parameters for the BASTA analysis
equalOrVaryingPopSizes <- "varying"
relaxedOrStrict <- "relaxed"
demeStructure <- "2Deme"
chainLength <- 300000000

# Build the XML file
buildXMLFile(demeStructure, equalOrVaryingPopSizes, path=simulationDirectory, date, sequenceInfo, chainLength, relaxedOrStrict, 
             sampleFromPrior=FALSE, estimateKappa=FALSE, genomeSize=genomeSize)

#### FUNCTIONS ####

buildXMLFile <- function(demeStructure, equalOrVaryingPopSizes="varying", path, date, selected, chainLength,
                         relaxedOrStrict="strict", sampleFromPrior=FALSE, estimateKappa=FALSE, 
                         genomeSize){
  
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
  fileLines <- addSequenceBlock(fileLines, selected, sampleFromPrior)
  
  # Add distribution block - if relaxed
  if(relaxedOrStrict == "relaxed"){
    fileLines <- addDistributionBlock(fileLines)
  }
  
  # Add constant sites block
  fileLines <- addConstantSiteCountsBlock(fileLines, nSNPs=nchar(selected$Sequence[1]), genomeSize)
  
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
                                     initialValue = 0.1, chainLength, relaxedOrStrict, estimateKappa)
  
  # Add the end to the XML file
  fileLines <- addEnd(fileLines)
  
  # Write the file out
  writeXMLFile(fileLines, paste(path, outputFileName, "/", outputFileName, ".xml", sep=""))
}

startBuildingOutputFileLines <- function(){
  
  # Create a vector to store the fileLines
  fileLines <- c()
  
  # Insert the first two lines
  fileLines[1] <- "<beast version='2.0' namespace='beast.evolution.alignment:beast.core:beast.core.parameter:beast.evolution.tree:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood:beast.evolution.tree:beast.math.distributions:multitypetreeVolz.distributions:multitypetreeVolz.operators:multitypetreeVolz.util'>"
  fileLines[2] <- ""
  
  return(fileLines)
}

addSequenceBlock <- function(fileLines, selected, sampleFromPrior){
  
  fileLines[length(fileLines) + 1] <- "\t<!-- Sequence Alignment -->"
  if(sampleFromPrior == FALSE){
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
  decimalDates <- decimal_date(selected$SamplingDate)
  
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
                                  initialValue, chainLength, relaxedOrStrict, estimateKappa){
  
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

addConstantSiteCountsBlock <- function(fileLines, nSNPs, genomeSize){
  
  ## Create the constant site counts
  n <- genomeSize - nSNPs
  counts <- round(c(n/4, n/4, n/4, n/4), digits=0)
  
  ## Print out block
  counts <- paste(counts, collapse=" ")
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- Constant Site Counts -->"
  fileLines[length(fileLines) + 1] <- paste("\t<data id='alignment' spec='FilteredAlignment' filter='-' data='@alignmentVar' constantSiteWeights='",
                                            counts, "'/>", sep="")
  
  return(fileLines)
}

# Preparing data for building XML

getSamplingInfoForEachSequence <- function(sequences, samplingInfo){
  
  # Initialise a data.frame to store species, sampling date and species
  output <- data.frame(IsolateID=NA, SamplingDate=NA, Species=NA, Sequence=NA, stringsAsFactors=FALSE)
  
  # Examine the information for each sample
  for(row in seq_len(nrow(samplingInfo))){
    
    # Get the ID
    id <- samplingInfo[row, "Isolate"]
    
    # Get the sequence for the current sample
    sequence <- sequences[[id]]
    
    # Check that present in sequences
    if(is.null(sequence)){
      cat(paste0("Sequence not present for: ", id, "\n"))
      next
    }
    
    # Get the sampling date
    samplingDate <- as.character(as.Date(paste0("15-", samplingInfo[row, "Date"]), format="%d-%b-%y"))
    
    # Get the species
    species <- ifelse(grepl(samplingInfo[row, "Notes"], pattern="Badger"), "badger", "cow")
    
    # Store the information for the current sample 
    output[row, ] <- c(id, samplingDate, species, sequence)
  }
  
  return(output)
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
  
  # Store the last sequence
  sequences[[name]] <- sequence
  
  return(sequences)
}
