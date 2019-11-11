#### Preparation ####

# Load libraries
library(ape) # Reading and phylogeny
library(geiger) # For the tips function in used in defineBranchColoursOfClade
library(basicPlotteR) # Because I made it! :-)

# Set the path 
path <- file.path("~", "Desktop", "Research", "RepublicOfIreland", "Mbovis", "WorkingWithDogSample", "")

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

#### Read in the sequence data ####

# Read in the FASTA file
fastaFile <- file.path(path, "vcfFiles", "sequences_Prox-10_01-11-2019.fasta")
sequences <- read.dna(fastaFile, format="fasta", as.character=TRUE)
nSites <- ncol(sequences)

# Plot the FASTA file
fastaPlotFile <- paste0(substr(fastaFile, 0, nchar(fastaFile) - 6), ".pdf")
plotFASTA(sequences, pdfFileName=fastaPlotFile, pdfHeight=14, pdfWidth=21, lineForSequenceNames=-3)

# Calculate coverage from sequences file
propNs <- data.frame("PropNs"=sapply(1:nrow(sequences), calculatePropNsInSequence, sequences))
rownames(propNs) <- rownames(sequences)

#### Build the phylogeny ####

# Build a phylogeny using RAxML
tree <- runRAXML(fastaFile, date="04-11-19", file.path(path, "vcfFiles", ""), outgroup="\\>Ref-1997")

# Remove Reference
tree <- drop.tip(tree, tree$tip.label[grepl(tree$tip.label, pattern=">Ref-1997")])

# Convert the branch lengths to SNPs
if(mean(tree$edge.length) < 1){
  tree$edge.length <- tree$edge.length * nSites
}

#### Get the tip information ####

# Read in the Wicklow sample IDs
wicklowInfoFile <- file.path(path, "Wicklow_Species_26-04-19.csv")
wicklowInfo <- read.table(wicklowInfoFile, header=TRUE, sep=",", stringsAsFactors=FALSE)
wicklowLinkFile <- file.path(path, "Wicklow_LINK_17-07-18.tsv")
wicklowLinkTable <- read.table(wicklowLinkFile, header=TRUE, stringsAsFactors=FALSE, sep="\t")

# Read in the Monaghan sample IDs
monaghanInfoFile <- file.path(path, "Monaghan_Species_08-08-19.csv")
monaghanInfo <- read.table(monaghanInfoFile, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Note the Wicklow and Monaghan vcf files
wicklowVCFs <- read.table(file.path(path, "Wicklow_vcfFiles.txt"), stringsAsFactors=FALSE)
monaghanVCFs <- read.table(file.path(path, "Monaghan_vcfFiles.txt"), stringsAsFactors=FALSE)

# Get the tip information
tipInfo <- getTipInfo(tree$tip.label, wicklowInfo, wicklowLinkTable, monaghanInfo, wicklowVCFs, monaghanVCFs)
tipInfo$PropNs <- propNs[tree$tip.label, "PropNs"]

# Note the ID of animals outside of Wicklow and Monaghan studies
other <- "1034_1.vcf.gz"
tipInfo[tipInfo$ID == ">1034_1.vcf.gz", c("Region", "Species")] <- "OTHER"

# Set the row names of the tip information
rownames(tipInfo) <- tipInfo$ID

#### Plot the phylogeny ####

################################################################################
# NOTE: NOT SURE WHAT 161 & 182 Wicklow VCFs ARE - IGNORED IN WICKLOW ANALYSES #
################################################################################

# Set the plotting margins
par(mar=c(2,0,0,10))

# Define the tip shapes and colours
tipShapesAndColours <- list("BADGER"=c("red", 19), "COW"=c("blue", 17),
                            "DEER"=c("black", 15), "OTHER"=c("green", 18),
                            "NA"=c("grey", 8))

### Plot the full phylogeny ###

# Define the tip shapes and colours based on species
tipShapes <- getTipShapeOrColourBasedOnSpecies(tipInfo, tipShapesAndColours, which="shape")
tipColours <- getTipShapeOrColourBasedOnSpecies(tipInfo, tipShapesAndColours, which="colour")

# Plot the phylogeny
plot.phylo(tree, show.tip.label=FALSE)

# Add tip shapes
tiplabels(pch=tipShapes, col=tipColours)

# Add region information
tiplabels(text=tipInfo$Region, offset=15, frame="none", cex=0.5, xpd=TRUE)

# Add a SNP scale
addSNPScale(position="bottom", size=10)

### Plot the Wicklow clade ###

# Note the node defining the clade and select it
node <- 200
clade <- extract.clade(tree, node)

# Define the tip shapes and colours
tipShapes <- getTipShapeOrColourBasedOnSpecies(tipInfo[clade$tip.label, ], tipShapesAndColours, which="shape")
tipColours <- getTipShapeOrColourBasedOnSpecies(tipInfo[clade$tip.label, ], tipShapesAndColours, which="colour")

# Plot the phylogeny
plot.phylo(clade, show.tip.label=FALSE)

# Add tip shapes
tiplabels(pch=tipShapes, col=tipColours)

# Add region information
tiplabels(text=tipInfo[clade$tip.label, "Region"], offset=1, frame="none", cex=0.5, xpd=TRUE)

# Add a SNP scale
addSNPScale(position="bottom", size=10)

#### FUNCTIONS ####

calculatePropNsInSequence <- function(sequenceIndex, sequences){
  
  # Count the number of Ns in the sequence
  nucleotideCounts <- table(sequences[sequenceIndex, ])
  names(nucleotideCounts) <- toupper(names(nucleotideCounts))
  
  # Note the number of Ns
  numberMissing <- 0
  if("N" %in% names(nucleotideCounts)){
    numberMissing <- nucleotideCounts[["N"]]
  }
  
  return(numberMissing / ncol(sequences))
}

getTipShapeOrColourBasedOnSpecies <- function(tipInfo, tipShapesAndColours, which,
                                              alpha=1){
  
  # Initialise a vector to store the shapes or colours
  output <- c()
  
  # Examine each tip
  for(row in seq_len(nrow(tipInfo))){
    
    # Check if Species available
    if(is.na(tipInfo[row, "Species"]) || tipInfo[row, "Species"] == "NA"){
      
      # Check if wanting shape or colour
      if(which == "shape"){
        output[row] <- as.numeric(tipShapesAndColours[["NA"]][2])
      }else if(which == "colour"){
        output[row] <- setAlpha(tipShapesAndColours[["NA"]][1], alpha)
      }else{
        cat(paste("Error! Option", which, "not recognised!"))
      }
      
    # If species available assign appropriate colour or shape
    }else{
      
      # Check if wanting shape or colour
      if(which == "shape"){
        output[row] <- as.numeric(tipShapesAndColours[[tipInfo[row, "Species"]]][2])
      }else if(which == "colour"){
        output[row] <- setAlpha(tipShapesAndColours[[tipInfo[row, "Species"]]][1], alpha)
      }else{
        cat(paste("Error! Option", which, "not recognised!"))
      }
    }
  }
  
  return(output)
}

getTipInfo <- function(tipLabels, wicklowInfo, wicklowLinkTable, monaghanInfo, wicklowVCFs, monaghanVCFs){
  
  # Initialise a dataframe to store the tip information
  tipInfo <- data.frame("ID"=tipLabels, "Region"=NA, "Species"=NA, stringsAsFactors=FALSE)
  
  # Examine each of the tip labels
  for(index in seq_along(tipLabels)){
    
    # Remove the ">" prefix
    tipLabel <- substr(tipLabels[index], 2, nchar(tipLabels[index]))
    
    # Check if from Wicklow or Monaghan datasets
    tipInfo[index, "Region"] <- checkIfWicklowOrMonaghan(tipLabel, wicklowVCFs, monaghanVCFs)
    
    # Parse the vcf file name to get tip label
    tipLabel <- strsplit(tipLabel, split="_")[[1]][1]
    
    # Remove a trailing "p"
    tipLabel <- gsub("p", "", tipLabel)
    
    # Check if WICKLOW
    if(tipInfo[index, "Region"] == "WICKLOW"){
      
      tipInfo[index, "Species"] <- getSpeciesWICKLOW(tipLabel, wicklowInfo, wicklowLinkTable)
      
    # Check if MONAGHAN
    }else if(tipInfo[index, "Region"] == "MONAGHAN"){
      
      tipInfo[index, "Species"] <- getSpeciesMONAGHAN(tipLabel, monaghanInfo)
    }
  }
  
  # Convert the species column to uppercase
  tipInfo$Species <- toupper(tipInfo$Species)
  
  return(tipInfo)
}

getSpeciesWICKLOW <- function(tipLabel, wicklowInfo, wicklowLinkTable){
  
  # Initialise a variable to store species
  species <- "NA"
  
  # Initialise a variable to store the aliquot code
  aliquotCode <- NA
  
  # Check if the current tip is associated with the original dataset
  if(grepl(tipLabel, pattern="-MBovis")){
    
    # Get the sequence number from the current tip label
    sequenceNumber <- strsplit(tipLabel, split="-")[[1]][1]
    
    # Find the row in the link table
    row <- which(wicklowLinkTable$Isolate.Code == sequenceNumber)
    
    # Get the current tips aliquot code
    if(length(row) != 0){
      aliquotCode <- wicklowLinkTable[row, "Aliquot"]
    }else{
      warning(paste("Problem for old WICKLOW batch. Couldn't find sequence number: ", sequenceNumber, " ", tipLabel, "\n"))
    }
    
  # Get information from most recent sequencing run
  }else{
    
    # Get the second part of the tip label - looks like an aliuot label without 00s
    aliquotCodePart <- strsplit(tipLabel, split="-")[[1]][2]
    
    # Find row that matches above part
    row <- which(grepl(wicklowInfo$Aliquot, pattern=aliquotCodePart))
    
    # Get the full aliquot code
    if(length(row) == 1){
      aliquotCode <- wicklowInfo[row, "Aliquot"]
    }else{
      cat(paste("Problem for new WICKLOW batch. Couldn't find aliquot part: ", aliquotCodePart, " (found ", length(row), " matches)\n\n"))
    }
  }
  
  # If aliquot code available then get the tip species and sampling year
  if(is.na(aliquotCode) == FALSE){
    
    # Get the row in the metadata table for the current aliquot code
    row <- which(wicklowInfo$Aliquot == aliquotCode)
    
    # Note the species and convert multiple "Cow labels to single "Cow label
    species <- wicklowInfo[row, "Species"]
    if(species %in% c("Heifer", "Steer", "Calf", "Bull", "Bovine")){
      species <- "Cow"
    }
  }
  
  return(species)
}

getSpeciesMONAGHAN <- function(tipLabel, monaghanInfo){

  # Initialise a variable to store species
  species <- "NA"
  
  # Build an aliquot code for the current isolate
  aliquot <- paste0("TB19-", paste(rep(0, 6-nchar(tipLabel)), collapse=""), tipLabel)

  # find the row in the animal info table
  tagRow <- which(monaghanInfo$Aliquot == aliquot)
    
  # Check that row was found
  if(length(tagRow) == 0){
    warning("Unable to find sampling information for: ", tipLabel, "\tAliquot: ", aliquot)
  }else if(length(tagRow) > 1){
    warning("Multiple entries in sampling information for: ", tipLabel, "\tAliquot: ", aliquot)
  }else{

    # Store the species
    species <-ifelse(grepl(monaghanInfo[tagRow, "Animal.ID"], pattern="^RR"),
                     "BADGER", "COW")
  }
  
  return(species)
}

getMonaghanTipInformation <- 

checkIfWicklowOrMonaghan <- function(tipLabel, wicklowVCFs, monaghanVCFs){
  
  # Create a variable to record whether tip is from wicklow of monaghan dataset
  dataset <- "UNKNOWN"
  
  # Check if from WICKLOW
  if(tipLabel %in% wicklowVCFs[, 1]){
    dataset <- "WICKLOW"
  }
  
  # Check if from MONAGHAN
  if(tipLabel %in% monaghanVCFs[, 1]){
    dataset <- "MONAGHAN"
  }
  
  return(dataset)
}

runRAXML <- function(fastaFile, date, path, nBootstraps=100, nThreads=10, outgroup=NULL, model="GTRCAT"){
  
  # Note the directory for the RAxML output files
  directory <- paste(path, "RAxML_", date, sep="")
  
  # Check if analyses already run
  alreadyRun <- dir.exists(directory)
  
  # If not already run, create output directory for RAxML
  if(alreadyRun == FALSE){
    suppressWarnings(dir.create(directory))
  }
  
  # Set the Working directory - this will be where the output files are dumped
  setwd(directory)
  
  # Build analysis name
  analysisName <- paste("RaxML-R_", date, sep="")
  
  # Check if already Run and just want to retrieve tree
  if(alreadyRun == FALSE){
    
    # Build the command
    seeds <- sample(1:100000000, size=2, replace=FALSE) # For parsimony tree and boostrapping
    
    if(is.null(outgroup)){
      command <- paste("raxmlHPC", 
                       " -f a", # Algorithm: Rapid boostrap inference
                       " -N ", nBootstraps,
                       " -T ", nThreads,
                       " -m ", model, " -V", # -V means no rate heterogenity
                       " -p ", seeds[1], " -x ", seeds[2], # Parsimony and boostrapping seeds
                       " -n ", analysisName,
                       " -s ", fastaFile, sep="")
    }else{
      command <- paste("raxmlHPC", 
                       " -f a", # Algorithm: Rapid boostrap inference
                       " -N ", nBootstraps,
                       " -T ", nThreads,
                       " -m ", model, " -V", # -V means no rate heterogenity
                       " -p ", seeds[1], " -x ", seeds[2], # Parsimony and boostrapping seeds
                       " -n ", analysisName,
                       " -s ", fastaFile, 
                       " -o ", outgroup, sep="")
    }
    
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
