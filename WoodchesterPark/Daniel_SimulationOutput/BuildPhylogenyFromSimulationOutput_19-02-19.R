#### Load libraries ####

library(ape) # Writing FASTA file
library(phangorn) # Tree building

#### Read in the simulation output ####

# Current date
date <- format(Sys.Date(), "%d-%m-%y")

# Set the path variable
path <- "/home/josephcrispell/Desktop/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/BASTA/Simulations/"

# NOTE the genome size
genomeSize <- 4345492

#### Examine the simulation output ####

# Read in a simulation output table
simFile <- paste0(path, "sim_119-02/sim_119-02.csv")
simOutput <- read.table(simFile, header=TRUE, sep=",", stringsAsFactors=FALSE)
simOutput$animal_ID <- gsub("_", "-", simOutput$animal_ID)

# REMOVE ROOT individual
simOutput <- simOutput[simOutput$animal_ID != "ROOT", ]

#### Build the sequences ####

# Identify the SNPs associated with each seed
seedSNPs <- getSeedSNPs(simOutput)

# Add a nucleotide sequence, based on the SNPs, for each individual
simOutput <- generateSequences(simOutput)

# Create a fasta file
fastaFile <- paste0(substr(simFile, 1, nchar(simFile)-4), ".fasta")
alignment <- createFasta(simOutput, fastaFile, seedSNPs=NULL)

#### Build phylogeny ####

# Build a phylogenetic tree using phangorn - neighbour joining and maximum likelihood
njTree <- buildNJTree(alignment)
mlTree <- buildAndBootStrapMLTree(initialTree=njTree, alignment=alignment, nBootstraps=100)

# Build a maximum likelihood phylogenetic tree using RAXML
raxmlTree <- runRAXML(fastaFile, date, nBootstraps=100, nThreads=6, path=paste0(path, "sim_119-02/"), alreadyRun=FALSE)

#### Plot the phylogeny ####

pdf(paste0(path, "sim_119-02/test.pdf"), height=40, width=40)

plot.phylo(raxmlTree, "fan")
nodelabels()

dev.off()

#############
# FUNCTIONS #
#############

buildNJTree <- function(alignment){
  
  # Build the distance matrix
  distanceMatrix <- dist.dna(as.DNAbin(alignment), model="JC69")

  # Build neighbour joining tree
  njTree <- nj(distanceMatrix)
  
  return(njTree)
}

buildAndBootStrapMLTree <- function(initialTree, alignment, nBootstraps=NULL){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #### Formatting alignment ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  # Convert alignment to phyDat object
  sequencesPhyDat <- as.phyDat(alignment)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #### Maximum Likelihood ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  cat("Running Maximum Likelihood tree estimation\n")
  
  # Compute likelihood of tree given sequences
  likelihoodObject <- pml(initialTree, sequencesPhyDat)
  
  # Set maximum likelihood controls
  controls <- pml.control(maxit=100000, trace=0)
  
  # Run maximum likelihood
  fittingOutput <- optim.pml(likelihoodObject, 
                             optNni = TRUE,       # Optimise topology
                             optInv = TRUE,       # Optimise proportion of variable sites
                             model = "JC",       # Use Jukes Cantor substitution model (this is model used to generate sequences - all sites equally likely)
                             rearrangement="NNI", # Nearest Neighbour Interchanges
                             control=controls)
  tree <- fittingOutput$tree
  
  #~~~~~~~~~~~~~~~~~~~~~#
  #### Bootstrapping ####
  #~~~~~~~~~~~~~~~~~~~~~#
  
  if(is.null(nBootstraps) == FALSE){
    cat("Running bootstrapping of Maximum Likelihood tree estimation\n")
    
    # Bootstrap the result of maximum likelihood
    bootstrapResults <- bootstrap.pml(fittingOutput, bs = nBootstraps, optNni = TRUE,
                                      jumble=TRUE)
    
    
    # Get phylogenetic tree with bootstrap values
    cat("Getting tree with bootstrap values. Finished..\n")
    tree <- plotBS(fittingOutput$tree, bootstrapResults, p = 50, type="phylogram")
  }

  return(tree)
}

runRAXML <- function(fastaFile, date, nBootstraps, nThreads, path, alreadyRun=FALSE){
  
  # Installing RAxML on Ubuntu:
  # sudo apt update
  # sudo apt install raxml
  
  # Note the RAxML directory name
  directory <- paste0(path, "RAxML_", date)
  
  # Create the directory
  dir.create(directory)
  
  # Set the Working directory - this will be where the output files are dumped
  setwd(directory)
  
  # Build analysis name
  analysisName <- paste("RaxML-R_", date, sep="")
  
  # Check if already run this
  if(alreadyRun == FALSE){
    
    # Create the RAxML directory for the output file
    suppressWarnings(dir.create(directory))
    
    # Build the command
    model <- "GTRCAT" # No rate heterogenity
    seeds <- sample(1:100000000, size=2, replace=FALSE) # For parsimony tree and boostrapping
    
    command <- paste("raxmlHPC", 
                     " -f a", # Algorithm: Rapid boostrap inference
                     " -N ", nBootstraps,
                     " -T ", nThreads,
                     " -m ", model, " -V", # -V means no rate heterogenity
                     " -p ", seeds[1], " -x ", seeds[2], # Parsimony and boostrapping seeds
                     " -n ", analysisName,
                     " -s ", fastaFile, sep="")
    system(command, intern=TRUE)
  }
  
  # Get the tree and read it in
  treeBS <- getTreeFileWithSupportValues(analysisName)
  
  return(treeBS)
}

getTreeFileWithSupportValues <- function(analysisName){
  
  # Get files in current working directory
  files <- list.files()
  
  # Select the tree file with BS support values
  treeBSFile <- files[grepl(files, pattern=paste("RAxML_bipartitions[.]", analysisName, sep="")) == TRUE]
  
  # Open the file
  treeBS <- read.tree(treeBSFile)
  
  return(treeBS)
}

createFasta <- function(simOutput, fileName, seedSNPs=NULL){
  
  # Initialise a vector to store the FASTA file lines
  fileLines <- c()
  
  # Create a matrix to store the sequences
  sequences <- matrix(NA, nrow=nrow(simOutput), ncol=nchar(simOutput[1, "Sequence"])-length(seedSNPs))
  rownames(sequences) <- simOutput$animal_ID
  
  # Examine each of the sequences and create the FASTA file lines
  for(row in seq_len(nrow(simOutput))){
    
    # Add in the sequence identifier
    fileLines[length(fileLines) + 1] <- paste0(">", simOutput[row, "animal_ID"])
    
    # Check if wanting to ignore seed SNPs
    if(is.null(seedSNPs)){
      
      # Store the current sequence in the matrix
      sequences[row, ] <- strsplit(simOutput[row, "Sequence"], split="")[[1]]
      
      # Add the current sequence
      fileLines[length(fileLines) + 1] <- simOutput[row, "Sequence"]
    }else{
      
      # Store the current sequence in the matrix
      sequences[row, ] <- strsplit(simOutput[row, "Sequence"], split="")[[1]][-seedSNPs]

      # Add the current sequence
      fileLines[length(fileLines) + 1] <- paste(sequences[row, ], collapse="")
    }
  }
  
  # Open the output FASTA file
  fileConnection <- file(fileName)
  
  # Print the file lines to file
  writeLines(fileLines, fileConnection)
  
  # Close the output FASTA file
  close(fileConnection)

  # Convert the sequences matrix to an alignment object
  alignment <- as.alignment(sequences)

  return(alignment)
}

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