#### Load libraries ####

library(ape) # Calculating genetic distances
library(phangorn) # NJ() function

#### Read in the simulation output ####

# Current date
date <- format(Sys.Date(), "%d-%m-%y")

# Set the path variable
path <- "/home/josephcrispell/Desktop/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/BASTA/Simulations/"

# NOTE the genome size
genomeSize <- 4345492

#### Examine the simulation output ####
fileName <- paste0(path, "sim_119-02/sim_48-01.csv")

# Create a directory specific to the current file
simulationDirectory <- paste0(path, substr(fileName,1,nchar(fileName)-4), "/")
dir.create(simulationDirectory)

# Read in a simulation output table
simFile <- paste0(path, "sequences/", fileName)
simOutput <- read.table(simFile, header=TRUE, sep=",", stringsAsFactors=FALSE)
simOutput$animal_ID <- gsub("_", "-", simOutput$animal_ID)

# REMOVE ROOT individual
simOutput <- simOutput[simOutput$animal_ID != "ROOT", ]

#### Build the sequences ####

# Identify the SNPs associated with each seed
seedSNPs <- getSeedSNPs(simOutput)

# Add a nucleotide sequence, based on the SNPs, for each individual
simOutput <- generateSequences(simOutput)

# Identify the non-seed SNPs
nonSeedSNPs <- c(1:nchar(simOutput[1, "Sequence"]))[-as.numeric(seedSNPs)]

# Create a fasta file

# Examine the genetic distances between the sequences
geneticDistancesWithoutSeedSNPs <- generateGeneticDistances(simOutput, sitesToIgnore=as.numeric(seedSNPs))
geneticDistances <- generateGeneticDistances(simOutput)

# Build and plot a phylogenetic tree
buildPhylogeny(geneticDistances)

#############
# FUNCTIONS #
#############

# Processing and examining data
buildPhylogeny <- function(fastaFile, treeBuildingTool="NeighbourJoining"){
  
  # Get and set the margins
  currentMar <- par("mar")
  par(mar=c(0,0,0,0))
  
  # Initialise the tree output
  tree <- null
  
  # Build neighbour joining tree if requested
  if(treeBuildingTool == "NeighbourJoining"){
    tree <- NJ(geneticDistances)
  }
  
  
  return(treel)
}

getLowerTriangle <- function(distances){
  
  # Initialise a vector to the the distances in the lower triangle
  output <- c()
  
  # Examine each row
  for(row in seq_len(nrow(distances))){
    
    # Examine each column
    for(col in seq_len(ncol(distances))){
      
      # Skip the upper and diagonal
      if(row >= col){
        next
      }
      
      # Store the current distance
      output[length(output) + 1] <- distances[row, col]
    }
  }
  
  return(output)
}

generateGeneticDistances <- function(simOutput, sitesToIgnore=NULL){
  
  # Initialise a list for the sequences
  sequences <- list()
  
  # Put each sequence into the list
  for(row in seq_len(nrow(simOutput))){
    
    if(is.null(sitesToIgnore)){
      sequences[[simOutput[row, "animal_ID"]]] <- strsplit(simOutput[row, "Sequence"], split="")[[1]]
    }else{
      sequences[[simOutput[row, "animal_ID"]]] <- strsplit(simOutput[row, "Sequence"], split="")[[1]][-sitesToIgnore]
    }
    
  }
  
  # Calculate the genetic distances
  geneticDistances <- dist.dna(as.DNAbin.alignment(as.alignment(sequences)), model="N", as.matrix=TRUE)
  
  return(geneticDistances)
}

getSeedSNPs <- function(simOutput){
  
  # Initialise variables to count how many seeds are present
  nBadgerSeeds <- 0
  nCattleSeeds <- 0
  
  # Initialise a vector to the seed SNPs
  seedSNPs <- c()
  
  # Initialise a vector to record how many SNPs each seed had
  nSeedSNPs <- c()
  
  # Examine the information for each simulated animal
  for(row in seq_len(nrow(simOutput))){
    
    # Skip non-seeds
    if(grepl(simOutput[row, "animal_ID"], pattern="seed") != TRUE){
      next
    }
    
    # Add to seed counts
    if(grepl(simOutput[row, "animal_ID"], pattern="Badger")){
      nBadgerSeeds <- nBadgerSeeds + 1
    }else{
      nCattleSeeds <- nCattleSeeds + 1
    }
    
    # Get the set of SNPs for the current seed
    snps <- strsplit(simOutput[row, "SNPs_infected"], split=";")[[1]]
    
    # Count the snps
    nSeedSNPs[row] <- length(snps)
    
    # Examine each snp
    for(snp in snps){
      
      # Add to growing vector
      if(snp %in% seedSNPs == FALSE){
        seedSNPs[length(seedSNPs) + 1] <- snp
      }else{
        cat("ERROR! SNP in seed: ", simOutput[row, "animal_ID"], " already found in a different seed.\n")
      }
    }
  }
  
  # Report the number of seeds
  cat(paste("Number of badger seeds =", nBadgerSeeds, 
            "\nNumber of cattle seeds =", nCattleSeeds, "\n"))
  
  # Plot the number of SNPs found in each seed
  hist(nSeedSNPs, las=1, xlab="Number of SNPs", main="Number of SNPs in seeds at start")
  
  return(seedSNPs)
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