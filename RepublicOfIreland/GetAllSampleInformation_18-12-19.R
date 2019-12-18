#### Preparation ####

# Load libraries

# Get today's date
today <- format(Sys.Date(), "%d-%m-%y")

# Note the path
path <- file.path("~", "storage", "Research", "RepublicOfIreland", "Mbovis")

#### Note the FASTQ nad VCF info ####

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
noVCFs[noVCFs$NumberMapped > 100000 & noVCFs$MappingProp > 0.1, c("SeqID", "Forward", "NumberMapped", "MappingProp", "County", "VCF", "Coverage", "Date")]

# Remove FASTQs if no VCF present
sampleInfo <- sampleInfo[is.na(sampleInfo$VCF) == FALSE, ]

#### Create a new directory structure for ROI data ####

# Copy all the files into a structured set of directories
copyFiles(sampleInfo, date="18-12-19")

# Copy each of the FASTQ and vcfFiles

#### Construct an aliquot for each sequence ID ####

#### Get an animal ID and species for each sequence ####

#### Create sample information file ####

#### FUNCTIONS ####

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
    fastqDirectory <- file.path(countyDirectory, paste0("Fastqs_", sampleInfo[row, "Date"]))
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
                          "Date"=rep(date, length(fastqs)),
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
