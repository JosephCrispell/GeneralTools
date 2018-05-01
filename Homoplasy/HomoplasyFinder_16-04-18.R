#-------------------#
#### Preparation ####
#-------------------#

# Packages
library(geiger) # For the tips function
library(ape) # read.dna()

#--------------------------------#
#### Run HomoplasyFinder in R ####
#--------------------------------#

# Read in the tree
pathToTreeFile <- file.path("C:", "Users", "Joseph Crisp", "Desktop", "UbuntuSharedFolder", "Homoplasy", "DataForTesting",
                            "example-AFTER_09-04-18.tree")
tree <- read.tree(pathToTreeFile)

# Read in the FASTA file
pathToFastaFile <- file.path("C:", "Users", "Joseph Crisp", "Desktop", "UbuntuSharedFolder", "Homoplasy", "DataForTesting",
                             "example_09-04-18.fasta")
sequences <- read.dna(pathToFastaFile, format="fasta", skip=1)

# Run HomoplasyFinder
results <- homoplasyFinder(tree, sequences)

#-----------------------------------#
#### Run HomoplasyFinder in Java ####
#-----------------------------------#

# Note the paths to the tree, FASTA, and HomoplasyFinder jar files
pathToTreeFile <- file.path("C:", "Users", "Joseph Crisp", "Desktop", "UbuntuSharedFolder", "Homoplasy",
                            "mlTree_27-03-18.tree")
pathToFastaFile <- file.path("C:", "Users", "Joseph Crisp", "Desktop", "UbuntuSharedFolder", "Homoplasy",
                             "sequences_Prox-10_24-03-2018.fasta")
pathToJarFile <- file.path("C:", "Users", "Joseph Crisp", "Desktop", "UbuntuSharedFolder", "Homoplasy",
                           "HomoplasyFinder_06-03-18.jar")

# Run HomoplasyFinder
results <- runHomoplasyFinderJavaTool(pathToJarFile, pathToFastaFile, pathToTreeFile)

#------------------------#
#### Plotting Results ####
#------------------------#

# Read in the tree
pathToTreeFile <- file.path("C:", "Users", "Joseph Crisp", "Desktop", "UbuntuSharedFolder", "Homoplasy",
                            "mlTree_27-03-18.tree")
tree <- read.tree(pathToTreeFile)

# Read in the FASTA file
pathToFastaFile <- file.path("C:", "Users", "Joseph Crisp", "Desktop", "UbuntuSharedFolder", "Homoplasy",
                             "sequences_Prox-10_24-03-2018.fasta")
sequences <- read.dna(pathToFastaFile, format="fasta", skip=1)

# Plot the homoplasy sites on the phylogeny
plotTreeAndHomoplasySites(tree, results, nSites)

#-----------------#
#### FUNCTIONS ####
#-----------------#

## Plotting results methods
plotTreeAndHomoplasySites <- function(tree, results){

  # Set the plotting margins
  par(mar=c(0,0,1,0.5))
  
  # Plot the phylogeny
  plot.phylo(tree, show.tip.label=TRUE, type="phylogram", align.tip.label=TRUE, tip.color=rgb(0,0,0,0),
             main="Homoplasies Identified")
  
  # Get the axisLimits
  axisLimits <- par("usr")
  
  # Add Scale bar
  xLength <- axisLimits[2] - axisLimits[1]
  yLength <- axisLimits[4] - axisLimits[3]
  points(x=c(axisLimits[1] + 0.1*xLength, axisLimits[1] + 0.2*xLength),
         y=c(axisLimits[3]+0.2*yLength, axisLimits[3]+0.2*yLength), type="l", lwd=3)
  text(x=axisLimits[1] + 0.15*xLength, y=axisLimits[3] +0.18*yLength, 
       labels=round(0.1*xLength, digits=3), cex=1, xpd=TRUE)
  
  # Note the locations of the tips
  tipCoordinates <- getTipCoordinates(tree$tip.label)
  maxCoords <- maxCoordinates(tipCoordinates)
  
  # Calculate width of space for nucleotide
  charWidth <- (axisLimits[2] - maxCoords[1]) / nrow(results)
  
  # Set nucleotide colours
  nucleotideColours <- list("A"="red", "C"="blue", "G"="cyan", "T"="orange")
  
  # Plot FASTA alignment beside tree
  for(homoplasyIndex in 1:nrow(results)){
    
    # Get an array of the nucleotides associated with the current position
    nucleotides = strsplit(results[homoplasyIndex, "Alleles"], split=",")[[1]]
    
    # Get an array of concatenated ID for each nucleotide
    concatenatedIsolates = strsplit(results[homoplasyIndex, "IsolatesForAlleles"], split=",")[[1]]
    
    # Examine each nucleotides
    for(nucleotideIndex in 1:length(nucleotides)){
      
      # Examine each isolate with the current nucleotide
      for(id in strsplit(concatenatedIsolates[nucleotideIndex], split=":")[[1]]){
        
        # Get XY coordinates for tip
        xy <- tipCoordinates[[id]]
        
        # Plot a polygon for the current tip's nucleotide at the current homoplasy's position
        polygon(x=c(maxCoords[1] + ((homoplasyIndex-1) * charWidth),
                    maxCoords[1] + ((homoplasyIndex-1) * charWidth),
                    maxCoords[1] + (homoplasyIndex * charWidth),
                    maxCoords[1] + (homoplasyIndex * charWidth)),
                y=c(xy[2], xy[2] + 1, xy[2] + 1, xy[2]),
                col=nucleotideColours[[nucleotides[nucleotideIndex]]],
                border=rgb(0,0,0,0), xpd=TRUE)
      }
    }
  }
    
  # Add separator lines
  for(homoplasyIndex in 1:nrow(results)){
    points(x=c(maxCoords[1] + ((homoplasyIndex-1) * charWidth),
               maxCoords[1] + ((homoplasyIndex-1) * charWidth)),
           y=c(axisLimits[3], axisLimits[4]),
           type="l", col="white")
  }
  
  # Note the positions of each homoplasy
  for(row in 1:nrow(results)){
    text(x=maxCoords[1] + ((row-1) * charWidth) + (0.5 * charWidth),
         y=maxCoords[2] + 0.02*yLength,
         labels=results[row, "Position"], cex=0.5, srt=90, xpd=TRUE)
  }
  
  # Reset plotting margins
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}

maxCoordinates <- function(tipCoordinates){
  
  max <- c(NA, NA)
  
  for(key in names(tipCoordinates)){
    
    if(is.na(max[1]) == TRUE || tipCoordinates[[key]][1] > max[1]){
      max[1] <- tipCoordinates[[key]][1]
    }
    
    if(is.na(max[2]) == TRUE || tipCoordinates[[key]][2] > max[2]){
      max[2] <- tipCoordinates[[key]][2]
    }
  }
  
  return(max)
}

getTipCoordinates <- function(tipLabels){
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tips <- list()
  for(i in 1:length(tipLabels)){
    tips[[as.character(tipLabels[i])]] <- c(lastPP$xx[i], lastPP$yy[i])
  }
  
  return(tips)
}

## HomoplasyFinder methods
runHomoplasyFinderJavaTool <- function(jarFile, fastaFile, treeFile, verbose=TRUE){
  
  # Get the current date
  date <- format(Sys.Date(), "%d-%m-%y")
  
  # Check for spaces in file paths - these can cause problems. If found surround path in quotes
  jarFile <- checkForSpaces(jarFile)
  fastaFile <- checkForSpaces(fastaFile)
  treeFile <- checkForSpaces(treeFile)
  
  # Check whether verbose requested
  verboseFlag = 0
  if(verbose){
    verboseFlag = 1
  }
  
  # Run the HomoplasyFinder jar file
  system(paste("java -jar", jarFile, verboseFlag, fastaFile, treeFile, sep=" "),
         ignore.stdout=FALSE)
  
  # Retrieve the output from HomoplasyFinder
  file <- paste("homoplasyReport_", date, ".txt", sep="")
  results <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                        check.names=FALSE)
  
  return(results)
}

checkForSpaces <- function(filePath){
  
  if(grepl(filePath, pattern=" ") == TRUE){
    filePath <- paste("\"", filePath, "\"", sep="")
  }
  
  return(filePath)
}

homoplasyFinder <- function(tree, sequencesDNABin, verbose=TRUE){
  
  #### Parse tree and FASTA
  
  # Note the isolates at the tips of each node in the phylogenetic tree
  nodes <- getNodes(tree)
  
  # Convert the DNAbin sequence alignment to alignment class
  sequences <- as.alignment(sequencesDNABin)
  
  # Record the alleles present in FASTA sequences
  alleles <- recordAllelesInPopulation(sequences, verbose)
  
  # Remove the constant sites (only one nucleotide present)
  removeConstantSites(alleles, sequences$nb, verbose)
  
  #### Assign alleles to nodes in phylogeny
  
  # Assign alelles to nodes in the phylogeny - where possible
  notAssigned <- assignAlellesToNodes(nodes, alleles, sequences$nam, verbose)
  
  #### Report the homopalsies identified
  
  results <- reportHomoplasiesIdentified(notAssigned, alleles)

  return(results)
}

removeConstantSites <- function(alleles, nSites, verbose=TRUE){
  
  nSitesRemoved <- 0
  
  nucleotides <- c('a', 'c', 'g', 't')
  for(position in 1:nSites){
    
    # Count the number of alleles at the current position
    count <- 0
    for(nucleotide in nucleotides){
      if(is.null(alleles[[paste(position, nucleotide, sep=":")]]) == FALSE){
        count <- count + 1
      }
    }
    # If only 1 allele present - remove it as it is a constant site (same across all isolates)
    if(count == 1){
      nSitesRemoved <- nSitesRemoved + 1
      alleles[[paste(position, 'N', sep=":")]] <- NULL
      for(nucleotide in nucleotides){
        if(is.null(alleles[[paste(position, nucleotide, sep=":")]]) == FALSE){
          alleles[[paste(position, nucleotide, sep=":")]] <- NULL
        }
      }     
    }
  }
  
  if(verbose){
    cat(paste("Removed", nSitesRemoved, "constant site(s).\n"))
  }
}

reportHomoplasiesIdentified <- function(notAssigned, alleles){
  
  # Get the positions where homoplasies were identified
  positions <- getPositions(notAssigned)
  
  # Initialise a data frame to store the homoplasy information
  output <- data.frame("Position"=names(positions), "Alleles"=NA, "IsolatesForAlleles"=NA, stringsAsFactors=FALSE)
  
  # Report each homoplasy
  for(row in 1:length(positions)){
    
    # Store the position
    output[row, "Position"] <- names(positions)[row]
    
    # Store the nucleotides associated with the current homoplasy's position
    nucleotides <- positions[[names(positions)[row]]]
    output[row, "Alleles"] <- paste(toupper(nucleotides), collapse=",")
    
    # Get the isolates associated with each nucleotide observed at the current site
    isolates <- paste(alleles[[paste(names(positions)[row], nucleotides[1], sep=":")]], collapse=":")
    for(i in 2:length(nucleotides)){
      
      isolates <- paste(isolates, paste(alleles[[paste(names(positions)[row], nucleotides[2], sep=":")]], collapse=":"), sep=",")
    }
    output[row, "IsolatesForAlleles"] <- isolates
  }
  
  return(output)
}

getPositions <- function(alleles){
  
  positions <- list()
  for(allele in alleles){
    
    parts <- strsplit(allele, split=":")[[1]]
    position <- parts[1]
    if(is.null(positions[[position]]) == FALSE){
      positions[[position]] <- c(positions[[position]], parts[2])
    }else{
      positions[[position]] <- c(parts[2])
    }
  }
  
  return(positions)
}

assignAlellesToNodes <- function(nodes, alleles, isolates, verbose){
  
  # Initialise a hashtable to record which alleles were assigned
  nodeAlleles <- list()
  
  # Initialise a vector to store the assigned alleles
  assigned <- c()

  # Examine each node
  for(i in 1:length(nodes)){

    # Progress
    if(verbose){
      cat(paste("\rAssigning alleles to node", i, "of", length(nodes)))
    }

    # Get the isolates at the tips of the current node
    tips <- nodes[[names(nodes)[i]]]
    
    # Get the rest of the isolates
    isolatesBelow <- isolates[isolates %in% tips == FALSE]

    # Check and see if these isolates are associated with an allele
    for(allele in names(alleles)){
      
      # Skip alleles that have already been assigned
      if(allele %in% assigned){
        next
      }
      
      # Get the isolates with an "N" at the current allele's position
      isolatesWithN <- getIsolatesWithNsAtPosition(allele, alleles)
      
      # Remove the isolates with an "N" from those associated with the current node
      tipsWithoutNs <- tips[tips %in% isolatesWithN == FALSE]
      isolatesBelowWithoutNs <- isolatesBelow[isolatesBelow %in% isolatesWithN == FALSE]
      
      # Get the isolates with the current allele
      isolatesWithAllele <- alleles[[allele]]
      
      # Compare the two sets of alleles - if they match exactly then current allele can be assigned to node
      if(areSetsOfIsolatesTheSame(tipsWithoutNs, isolatesWithAllele) == TRUE){

        nodeAlleles[[names(nodes)[i]]] <- c(nodeAlleles[[names(nodes)[i]]], allele)
        assigned[length(assigned) + 1] <- allele
      
      }else if(areSetsOfIsolatesTheSame(isolatesBelowWithoutNs, isolatesWithAllele) == TRUE){

        nodeAlleles[[names(nodes)[i]]] <- c(nodeAlleles[[names(nodes)[i]]], allele)
        assigned[length(assigned) + 1] <- allele
      }
    }
  }
  
  if(verbose){
    cat("\rFinished assigning alleles to nodes.\t\t\t\t\n")
  }
  
  # Note the alleles that weren't assigned
  notAssigned <- names(alleles)[names(alleles) %in% assigned == FALSE]
  
  # Remove the N alleles
  notAssigned <- notAssigned[grepl(notAssigned, pattern=":n") == FALSE]
  
  return(notAssigned)
}

getIsolatesWithNsAtPosition <- function(allele, alleles){
  
  # Build the allele with an N
  key <- paste(strsplit(allele, split=":")[[1]][1], "n", sep=":")
  
  # Check whether any isolates with an N were noted
  isolates <- c()
  if(is.null(alleles[[key]]) == FALSE){
    isolates <- alleles[[key]]
  }
  
  return(isolates)
}

areSetsOfIsolatesTheSame <- function(a, b){
  result <- FALSE
  if(length(a) == length(b) && length(intersect(a, b)) == length(a)){
    result <- TRUE
  }
  return(result)
}

getNodes <- function(tree){
  nodes <- list()
  
  # Nodes number 1:nTips and then onwards: two methods to calculate total number:
  # - Number of edges + 1
  # - Number of tips plus number internal nodes: length(tree$tip.label) + tree$Nnode
  for(node in 1:(length(tree$tip.label) + tree$Nnode)){
    nodes[[as.character(node)]] <- tips(tree, node)
  }
  
  return(nodes)
}

recordAllelesInPopulation <- function(sequences, verbose){
  
  # Initialise a list to store the alleles found and the sequences they are found in
  alleles <- list()
  
  # Examine each of the sequences
  for(i in 1:length(sequences$nam)){
    
    # Progress
    if(verbose){
      cat(paste("\rReading alleles in sequence", i, "of", length(sequences$nam)))
    }
    
    # Split the current sequence into nucleotides
    nucleotides <- strsplit(sequences$seq[i], split="")[[1]]
    
    # Examine each nucleotide in the current sequence
    for(pos in 1:length(nucleotides)){
      
      # Build an allele ID
      id <- paste(pos, nucleotides[pos], sep=":")
      
      # Check if encountered this allele before
      if(is.null(alleles[[id]]) == FALSE){
        
        alleles[[id]] <- c(alleles[[id]], sequences$nam[i])
      }else{
        alleles[[id]] <- c(sequences$nam[i])
      }
    }
  }
  
  if(verbose){
    cat("\rFinished reading alleles from sequences.\t\t\t\t\n")
  }
  
  return(alleles)
}
