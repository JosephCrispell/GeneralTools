# Load libraries
library(basicPlotteR)

# Set the path variable
path <- "/home/josephcrispell/Desktop/Research/RepublicOfIreland/Mbovis/Wicklow/vcfFiles/"

# Read in the coverage table
table <- read.table(paste0(path, "isolateCoverageSummary_DP-20_19-03-2019.txt"), header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Calculate some summary statistics
median(table$PercentageCoverage)
median(table$MeanDepth)
quantile(table$PercentageCoverage, probs=c(0.025, 0.975))

# Read in the FASTA file
sequences <- readFASTA(paste0(path, "sequences_Prox-10_19-03-2019.fasta"))

# Count alleles at each position
alleleCounts <- countAlleles(sequences)

# Count the number of alleles present at each site
length(which(alleleCounts$NumberAlleles > 1))

#### FUNCTIONS ####

countAlleles <- function(sequences){
  
  # Get the sequence ids
  ids <- names(sequences)
  
  # Initialise a tbale ot store the counts
  counts <- data.frame("A"=NA, "C"=NA, "G"=NA, "T"=NA, "N"=NA)
  
  # Examine each position
  for(position in seq_along(sequences[[1]])){
    
    # Initialise a list to count each allele
    alleles <- list("A"=0, "C"=0, "G"=0, "T"=0, "N"=0)
    
    # Examine each of the sequences
    for(id in ids){
      
      # Count the allele at the current position
      alleles[[sequences[[id]][position]]] <- alleles[[sequences[[id]][position]]] + 1
    }
    
    # Store the allele counts
    counts[position, "A"] <- alleles[["A"]]
    counts[position, "C"] <- alleles[["C"]]
    counts[position, "G"] <- alleles[["G"]]
    counts[position, "T"] <- alleles[["T"]]
    counts[position, "N"] <- alleles[["N"]]
    
    # Monitor progress
    progress(position, length(sequences[[1]]))
  }
  
  # Count the number of alleles at each position
  allelePresentAbsent <- counts[, c("A", "C", "G", "T")]
  allelePresentAbsent[allelePresentAbsent > 0] <- 1
  counts$NumberAlleles <- rowSums(allelePresentAbsent)
  
  return(counts)
}

readFASTA <- function(fastaFile, skip=1){
  
  # Open a connection to a file to read (open="r")
  connection <- file(fastaFile, open="r")
  
  # Get all lines from file and store in vector
  fileLines <- readLines(connection)
  
  # Close file connection
  close(connection)
  
  # Skip the first X lines
  fileLines <- fileLines[-skip]
  
  # Initialise a list to store the sequences
  sequences <- list()
  
  # Loop through each of the lines in file
  for(line in fileLines){
    
    # Check if found new sequence
    if(grepl(line, pattern="^>")){
      
      # Get the sequence name
      name <- substr(line, 2, nchar(line))
      
      # Start storing the sequence
      sequences[[name]] <- c()
      
      # Else if haven't found new sequence - continue storing previous sequence
    }else{
      sequences[[name]] <- c(sequences[[name]], strsplit(line, split="")[[1]])
    }
  }
  
  return(sequences)
}
