#### Preparation ####

# Load libraries

# Get today's date
today <- format(Sys.Date(), "%d-%m-%y")

# Note the path
path <- file.path("~", "storage", "Research", "RepublicOfIreland", "Mbovis")

# Create a directory to store all data
createDirectory(file.path(path, "ROI_18-12-19"))

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


#### Note the VCF file names ####

#### Construct single sample information table ####

#### FUNCTIONS ####

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
                        "County"=rep(directory, length(vcfFiles)),
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
