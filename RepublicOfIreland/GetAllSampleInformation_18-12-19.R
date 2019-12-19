#### Preparation ####

# Load libraries
library(phyloHelpeR) # Running RAXML
library(ape) # Writing tree file

# Get today's date
today <- format(Sys.Date(), "%d-%m-%y")

# Note the path
path <- file.path("~", "storage", "Research", "RepublicOfIreland", "Mbovis")

#### Note the FASTQ and VCF info ####

# Get the FASTQ files
fastqs <- getFASTQs(path, directory="Monaghan", date="29-07-19")
fastqs <- rbind(fastqs, getFASTQs(path, directory="Monaghan", date="24-09-19"))
fastqs <- rbind(fastqs, getFASTQs(path, directory="Monaghan", date="16-12-19"))
fastqs <- rbind(fastqs, getFASTQs(path, directory="Wicklow", date="07-01-18"))
fastqs <- rbind(fastqs, getFASTQs(path, directory="Wicklow", date="09-07-18", ignore="Mbovis"))
fastqs <- rbind(fastqs, getFASTQs(path, directory="Wicklow", date="15-03-19", ignore="Bovis76"))

# Get the VCF files
vcfs <- getVCFs(path, "Monaghan")
vcfs <- rbind(vcfs, getVCFs(path, "Wicklow"))

# Combine the FASTQ and VCF information into a single table
sampleInfo <- merge(x=fastqs, y=vcfs, by="SeqID", all=TRUE) # Outer join

# Note FASTQ files that can't find VCF - duplicates or poor coverage
noVCFs <- sampleInfo[is.na(sampleInfo$VCF), ]
noVCFs[noVCFs$NumberMapped > 100000 & noVCFs$MappingProp > 0.1, 
       c("SeqID", "Forward", "NumberMapped", "MappingProp", "County", "VCF", "Coverage", "SequencingDate")]

# Remove FASTQs if no VCF present
sampleInfo <- sampleInfo[is.na(sampleInfo$VCF) == FALSE, ]

#### Create a new directory structure for ROI data ####

# Copy all the files into a structured set of directories
copyFiles(sampleInfo, date="18-12-19")

#### Construct a phylogeny for each set of VCFs ####

# Wicklow phylogeny
fasta <- file.path(path , "ROI_18-12-19", "Wicklow", "vcfFiles", "sequences_Prox-10_18-12-2019.fasta")
tree <- runRAXML(fastaFile=fasta, date="18-12-19", path=file.path(path , "ROI_18-12-19", "Wicklow/"),
                 nThreads=10, nBootstraps=100, outgroup="\\>Ref-1997")
write.tree(tree, file.path(path , "ROI_18-12-19", "Wicklow", "mlTree_WICKLOW_18-12-2019.tree"))

# Monaghan phylogeny
fasta <- file.path(path , "ROI_18-12-19", "Monaghan", "vcfFiles", "sequences_Prox-10_18-12-2019.fasta")
tree <- runRAXML(fastaFile=fasta, date="18-12-19", path=file.path(path , "ROI_18-12-19", "Monaghan/"),
                 nThreads=10, nBootstraps=100, outgroup="\\>Ref-1997")
write.tree(tree, file.path(path , "ROI_18-12-19", "Monaghan", "mlTree_MONAGHAN_18-12-2019.tree"))

# All ROI phylogeny
fasta <- file.path(path , "ROI_18-12-19", "vcfFiles", "sequences_Prox-10_18-12-2019.fasta")
tree <- runRAXML(fastaFile=fasta, date="18-12-19", path=file.path(path , "ROI_18-12-19/"),
                 nThreads=10, nBootstraps=100, outgroup="\\>Ref-1997")
write.tree(tree, file.path(path , "ROI_18-12-19", "mlTree_ROI_18-12-2019.tree"))

#### Construct an aliquot for each sequence ID ####

# Note the two genomes from Northern Ireland
sampleInfo[sampleInfo$SeqID %in% c("161-MBovis", "182-MBovis"), "County"] <- "SomewhereInNorthernIreland"

# Note the various sample info files for Wicklow data
file <- file.path(path, "Wicklow", "Mbovis_SamplingInfo_17-07-18.tsv")
wicklowLinkTable <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep="\t")
file <- file.path(path, "Wicklow", "IsolateSpeciesAndYear_26-04-19.csv")
wicklowMetadata <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep=",")

# Note the various sample info files for Monaghan data
file <- file.path(path, "ROI_18-12-19", "EartagsAndSettIDs_REMOVE_17-12-19.csv")
monaghanAnimalTags <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)
file <- file.path(path, "ROI_18-12-19", "Animal_HerdIDs_REMOVE_17-12-19.csv")
monaghanHerdIDs <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Create aliquot and get an animal ID and species for each sequence
sampleInfo <- createAliquotsAndLinkToMetadata(sampleInfo, wicklowLinkTable, wicklowMetadata, monaghanAnimalTags)

# Check the samples I can't find species for
noSpecies <- sampleInfo[is.na(sampleInfo$Species), ]

# Check the samples I can't find animal ID
noAnimalID <- sampleInfo[is.na(sampleInfo$AnimalID), ]

#### Store the sample information ####

# Remove the path columns
sampleInfo <- sampleInfo[, grepl(colnames(sampleInfo), pattern="PathTo") == FALSE]

# Store the sample information table with and without animal IDs
write.table(sampleInfo, file=file.path(path, "ROI_18-12-19", "sampleInformation_18-12-19.csv"), quote=FALSE, sep=",", row.names=FALSE)
write.table(sampleInfo[, -ncol(sampleInfo)], file=file.path(path, "ROI_18-12-19", "sampleInformation_NoAnimalIDs_18-12-19.csv"),
            quote=FALSE, sep=",", row.names=FALSE)

#### FUNCTIONS ####

createAliquotsAndLinkToMetadata <- function(sampleInfo, wicklowLinkTable, wicklowMetadata, 
                                            monaghanAnimalTags){
  
  # Create additional columns in sample information to store Aliquot, species, and animal ID
  sampleInfo$Aliquot <- NA
  sampleInfo$Species <- NA
  sampleInfo$Date <- NA
  sampleInfo$AnimalID <- NA
  
  # Examine the information for each sequence
  for(row in 1:nrow(sampleInfo)){
    
    # Check if Wicklow
    if(sampleInfo[row, "County"] == "Wicklow"){
      
      # Check if require link table
      if(grepl(sampleInfo[row, "SeqID"], pattern="-MBovis")){
        
        # Get the sequence number from the curren tip label
        sequenceNumber <- strsplit(sampleInfo[row, "SeqID"], split="-")[[1]][1]
        
        # Find the row in the link table
        linkTableRow <- which(wicklowLinkTable$Isolate.Code == sequenceNumber)
        
        # Get the current tips aliquot code
        if(length(linkTableRow) != 0){
          sampleInfo[row, "Aliquot"] <- wicklowLinkTable[linkTableRow, "Aliquot"]
        }else{
          warning(paste("Couldn't find sequence number: ", sequenceNumber, " ", sampleInfo[row, "SeqID"], "\n"))
        }
        
      # Otherwise get from metadata
      }else{
        
        # Get the second part of the tip label - looks like an aliuot label without 00s
        aliquotCodePart <- strsplit(sampleInfo[row, "SeqID"], split="-")[[1]][2]
        
        # Find row that matches above part
        metadataRow <- which(grepl(wicklowMetadata$Aliquot, pattern=aliquotCodePart))
        
        # Get the full aliquot code
        if(length(row) == 1){
          sampleInfo[row, "Aliquot"] <- wicklowMetadata[metadataRow, "Aliquot"]
        }else{
          warning(paste("Couldn't find aliquot part: ", aliquotCodePart, " (found ", length(metadataRow), " matches)\n\n"))
        }
      }
      
      # Check if found aliquot
      if(is.na(sampleInfo[row, "Aliquot"]) == FALSE){
        
        # Get the row in the metadata table for the current aliquot code
        metadataRow <- which(wicklowMetadata$Aliquot == sampleInfo[row, "Aliquot"])
        
        # Note the species and convert multiple "Cow labels to single "Cow label
        sampleInfo[row, "Species"] <- wicklowMetadata[metadataRow, "Species"]
        if(sampleInfo[row, "Species"] %in% c("Heifer", "Steer", "Calf", "Bull", "Bovine")){
          sampleInfo[row, "Species"] <- "Cow"
        }
        
        # Note the sampling date
        if(wicklowMetadata[metadataRow, "Received.Test.date"] != ""){
          sampleInfo[row, "Date"] <- as.character(as.Date(wicklowMetadata[metadataRow, "Received.Test.date"], format="%d/%m/%Y"))
        }
      }
    
    # Check if Monaghan
    }else if(sampleInfo[row, "County"] == "Monaghan"){
      
      # Split the tip label if dash
      aliquot <- gsub("p", "", sampleInfo[row, "SeqID"])
      if(grepl(aliquot, pattern="-")){
        aliquot <- strsplit(aliquot, split="-")[[1]][2]
      }
      
      # Remove trailing p if present
      
      # Build an aliquot code for the current isolate
      aliquot <- paste0("TB19-", paste(rep(0, 6-nchar(aliquot)), collapse=""), aliquot)
      sampleInfo[row, "Aliquot"] <- aliquot
      
      # Find the row in the sample information table for the current tip
      tagRow <- which(monaghanAnimalTags$Aliquot == aliquot)
      
      # Check that row was found
      if(length(tagRow) == 0){
        warning("Unable to find sampling information for: ", sampleInfo[row, "SeqID"], "\tAliquot: ", aliquot)
      }else if(length(tagRow) > 1){
        warning("Multiple entries in sampling information for: ", sampleInfo[row, "SeqID"], "\tAliquot: ", aliquot)
      }else{
        # Store the animal ID
        sampleInfo[row, "AnimalID"] <- monaghanAnimalTags[tagRow, "Animal.ID"]
        
        # Store the species of the current tip
        sampleInfo[row, "Species"] <- ifelse(grepl(sampleInfo[row, "AnimalID"], pattern="^RR"), "Badger", "Cow")
      }
    }
  }
  
  return(sampleInfo)
}

copyFiles <- function(sampleInfo, date){
  
  # Create a directory to store all data
  home <- file.path(path, paste0("ROI_", date))
  createDirectory(home)
  
  # Examine each of the file information lines
  for(row in 1:nrow(sampleInfo)){
    
    # Create county directory, if doesn't already exist
    countyDirectory <- file.path(home, sampleInfo[row, "County"])
    createDirectory(countyDirectory)
    
    # Create folder for FASTQ files, if doesn't already exist
    fastqDirectory <- file.path(countyDirectory, paste0("Fastqs_", sampleInfo[row, "SequencingDate"]))
    createDirectory(fastqDirectory)
    
    # Create folder for vcf file, if doesn't already exist
    vcfDirectory <- file.path(countyDirectory, "vcfFiles")
    createDirectory(vcfDirectory)
    
    # Copy the FASTQ files
    forward <- file.path(sampleInfo[row, "PathToFASTQ"], sampleInfo[row, "Forward"])
    reverse <- file.path(sampleInfo[row, "PathToFASTQ"], sampleInfo[row, "Reverse"])
    file.copy(forward, fastqDirectory)
    file.copy(reverse, fastqDirectory)
    
    # Copy the VCF file into county folder
    vcf <- file.path(sampleInfo[row, "PathToVCF"], sampleInfo[row, "VCF"])
    file.copy(vcf, vcfDirectory)
    
    # Copy the VCF file into All folder
    vcfDirectory <- file.path(home, "vcfFiles")
    createDirectory(vcfDirectory)
    file.copy(vcf, vcfDirectory)
  }
}

getVCFs <- function(path, directory){
  
  # Get the VCF files
  folder <- file.path(path, directory, "vcfFiles")
  vcfFiles <- getFilesInDirectory(folder, matchToKeep=".vcf.gz")
  
  # Get the coverage information
  file <- getFilesInDirectory(folder, matchToKeep="isolateCoverageSummary_DP-20", matchToRemove="pdf", keepPath=TRUE)
  coverage <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  
  # Create a data.frame to store VCF information
  vcfInfo <- data.frame("SeqID"=rep(NA, length(vcfFiles)),
                        "VCF"=vcfFiles,
                        "Coverage"=rep(NA, length(vcfFiles)),
                        "AverageDepth"=rep(NA, length(vcfFiles)),
                        "PathToVCF"=rep(folder, length(vcfFiles)),
                        stringsAsFactors=FALSE)
  
  # Examine each VCF file
  for(index in seq_along(vcfFiles)){
    
    # Get the sequence ID
    vcfInfo[index, "SeqID"] <- strsplit(vcfFiles[index], split="_")[[1]][1]
    
    # Store the coverage and average read depth for mapping to M. bovis
    coverageRow <- which(coverage$IsolateID == vcfFiles[index])
    vcfInfo[index, "Coverage"] <- coverage[coverageRow, "PercentageCoverage"]
    vcfInfo[index, "AverageDepth"] <- coverage[coverageRow, "MeanDepth"]
  }
  
  return(vcfInfo)
}

getFASTQs <- function(path, directory, date, ignore=NULL){
  
  # Get the FORWARD FASTQ files
  folder <- paste0(path, "/", directory, "/", "Fastqs_", date, "/")
  fastqs <- getFilesInDirectory(folder, matchToKeep="R1_001.fastq.gz", matchToRemove=ignore, keepPath=FALSE)
  
  # Get the mapping information
  file <- getFilesInDirectory(folder, matchToKeep="isolateMappingSummary", matchToRemove="pdf", keepPath=TRUE)
  mapping <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  
  # Create a data.frame to store FASTQ information
  fastqInfo <- data.frame("SeqID"=rep(NA, length(fastqs)),
                          "Forward"=fastqs,
                          "Reverse"=rep(NA, length(fastqs)),
                          "NumberMapped"=rep(NA, length(fastqs)),
                          "MappingProp"=rep(NA, length(fastqs)),
                          "SequencingDate"=rep(date, length(fastqs)),
                          "County"=rep(directory, length(fastqs)),
                          "PathToFASTQ"=rep(folder, length(fastqs)),
                          stringsAsFactors=FALSE)
  
  # Examine each FORWARD FASTQ file
  for(index in seq_along(fastqs)){
    
    # Create a REVERSE fastq name
    fastqInfo[index, "Reverse"] <- gsub(pattern="R1", replacement="R2", x=fastqs[index])
    
    # Get the sequence ID
    fastqInfo[index, "SeqID"] <- strsplit(fastqs[index], split="_")[[1]][1]
    
    # Store the proportion reads mapped to M. bovis
    mappingRow <- which(mapping$Isolate == fastqInfo[index, "SeqID"])
    propMapped <- mapping[mappingRow, "NumberMappedReads"] / 
      (mapping[mappingRow, "NumberMappedReads"] + mapping[mappingRow, "NumberUnmappedReads"])
    fastqInfo[index, "MappingProp"] <- propMapped
    fastqInfo[index, "NumberMapped"] <- mapping[mappingRow, "NumberMappedReads"]
  }
  
  return(fastqInfo)
}

createDirectory <- function(directory){

  # If doesn't already exist, create directory
  if(dir.exists(directory) == FALSE){
    dir.create(directory)
  }
}

getFilesInDirectory <- function(directory, matchToKeep=NULL, matchToRemove=NULL, keepPath=FALSE){
  
  # List the files in the directory
  files <- list.files(path=directory)
  
  # Select files with file match pattern if provided
  if(is.null(matchToKeep) == FALSE){
    
    files <- files[grepl(files, pattern=matchToKeep)]
  }
  
  # Remove if match pattern for files to remove provided
  if(is.null(matchToRemove) == FALSE){
    
    files <- files[grepl(files, pattern=matchToRemove) == FALSE]
  }
  
  # Check if want to add path
  if(keepPath){
    files <- ifelse(endsWith(directory, "/"), paste0(directory, files), paste0(directory, "/", files))
  }
  
  return(files)
}
