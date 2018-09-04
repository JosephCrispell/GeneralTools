#### Load Libraries ####
library(devtools)
library(homoplasyFinder)
library(ape)

#### Note FASTA and NEWICK file locations ####

# Create path variable
path <- "/home/josephcrispell/Desktop/Research/Grandjean2017ConvergentEvolution/"

# Note FASTA file
fastaFile <- paste0(path, "FASTQs/vcfFiles/sequences_Prox-0_02-09-2018.fasta")

# Note NEWICK file - NOTE: removed ">" from tip labels
newickFile <- paste0(path, "RAxML_02-09-2018/RAxML_bipartitions.RaxML-R_02-09-18_rmGT")

#### Run HomoplasyFinder ####

# Run homoplasyFinder
inconsistentSites <- runHomoplasyFinderInJava(treeFile=newickFile, fastaFile=fastaFile, path=path)

# Read in consistency index report
results <- read.table(paste0(path, "HomoplasyFinderOutput/consistencyIndexReport_03-09-18.txt"), header=TRUE, sep="\t")

#### Identify genome positions of inconsistent sites ####

# Read in the FASTA positions file
fastaPositions <- read.table(paste0(path, "FASTQs/vcfFiles/fastaPositions_Prox-0_02-09-2018.txt"), 
                             header=TRUE)

# Add genome position to results table
results$GenomePosition <- fastaPositions[results$Position, 1]

#### Read in table of homoplasies found by Grandjean et al. ####

# Read in table of most homoplastic polymorphisms
mostHomoplastic <- read.table(paste0(path, "SupplementaryInformation/MostHomoplasticPolymorphisms_28-08-18.csv"), 
                              header=TRUE, sep=",")

# Read in table of phyC
phycHomoplasies <- read.table(paste0(path, "SupplementaryInformation/PhyC_Homoplasies_27-08-18.csv"), 
                              header=TRUE, sep=",")

# How many homoplasies are common in these lists?
length(intersect(phycHomoplasies$Reference.position, mostHomoplastic$Reference.ID))

#### How many of these did homoplasyFinder identify?? ####

# Read in the FASTA file
sequences <- read.dna(fastaFile, format="fasta", as.character=TRUE)

# Note the positions not identified from most homoplastic table
notFoundMostHomoplastic <- setdiff(mostHomoplastic$Reference.ID, results$GenomePosition)
alleleCountsMostHomoplastic <- getAlleleCountsForPositionsNotFound(notFoundMostHomoplastic, fastaPositions, sequences)

# Note the positions not identified from phyc homoplasies
notFoundPhyc <- setdiff(phycHomoplasies$Reference.position, results$GenomePosition)
alleleCountsPhyc <- getAlleleCountsForPositionsNotFound(notFoundPhyc, fastaPositions, sequences)

#### Examine the quality of the isolates with an N at the positions not found ####

# Read in the merged vcf file
mergedTable <- readInMergedVCFsfile(paste0(path, "FASTQs/vcfFiles/merged_02-09-2018.txt"))

# Look at the quality of iolates for particular positions that weren't found
position <- 3935454
isolates <- getIsolatesNamesWithAlleleAtPosition("n", position=position, fastaPositions, sequences)
qualityInfo <- getQualityInformationForIsolatesAtPosition(isolates, position=position, mergedTable)
qualityInfo[grepl(qualityInfo, pattern="C[.]|-------------------") == FALSE]

#### FUNCTIONS ####

readInMergedVCFsfile <- function(file){
  
  # Read in the merged vcf file
  mergedTable <- read.table(file, skip=28, sep=":", header=TRUE, stringsAsFactors=FALSE, comment.char="~", check.names=FALSE)
  
  # Split the isolate names into an array
  isolates <- c(strsplit(colnames(mergedTable)[1], split="\t")[[1]][3], colnames(mergedTable)[-1])
  
  # Create the Position column at start of merged table
  mergedTable$Position <- NA
  colnames(mergedTable) <- c(isolates, "Position")
  mergedTable <- mergedTable[, c(ncol(mergedTable), 1:(ncol(mergedTable) - 1))]

  # Split the first column into CHROM, Position, and Isolate1QualityInfo
  for(row in seq_len(nrow(mergedTable))){
    
    # Print progress
    cat(paste("\rReading row:", row, "of", nrow(mergedTable), sep=" "))
    
    # Split the quality information
    mergedTable[row, c(1,2)] <- strsplit(mergedTable[row, 2], split="\t")[[1]][c(2,3)]
  }
  cat("\rFinished :-)\n")
  
  return(mergedTable)
}

getQualityInformationForIsolatesAtPosition <- function(isolates, position, mergedTable){
  
  # Find row for position
  row <- which(mergedTable$Position == position)
  
  # Return quality information for isolates
  return(as.character(mergedTable[row, isolates]))
}

getIsolatesNamesWithAlleleAtPosition <- function(allele, position, fastaPositions, sequences){
  
  # Get the position in the FASTA file
  fastaPosition <- seq_len(nrow(fastaPositions))[fastaPositions$Position == position]
  return(rownames(sequences)[sequences[, fastaPosition] == "n"])
}

getAlleleCountsForPositionsNotFound <- function(notFound, fastaPositions, sequences){
  
  # Intialise a matrix to store the allele counts at each FASTA position
  alleleCounts <- matrix(0, nrow=length(notFound), ncol=6)
  
  # Add column names
  colnames(alleleCounts) <- c("Position", "A", "C", "G", "T", "N")
  
  # Fill in positions column
  alleleCounts[, 1] <- notFound
  
  # Initialise a list to note the columns of each allele
  alleleColumns <- list("a"=2, "c"=3, "g"=4, "t"=5, "n"=6)
  
  # Examine each of the positions not found
  for(row in seq_along(notFound)){
    
    # Get position in FASTA
    positionInFASTA <- seq_len(nrow(fastaPositions))[fastaPositions$Position == notFound[row]]
    
    # Check if this position is present in FASTA - if not skip
    if(length(positionInFASTA) == 0){
      next
    }
    
    # Count the numbers of As, Cs, Gs, Ts, and Ns
    for(i in seq_len(nrow(sequences))){
      alleleCounts[row, alleleColumns[[sequences[i, positionInFASTA]]]] <- 
        alleleCounts[row, alleleColumns[[sequences[i, positionInFASTA]]]] + 1
    }
  }
  
  return(alleleCounts)
}

