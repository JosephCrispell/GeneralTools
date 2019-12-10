#### Read in the tip information tables ####

# Note today's date
date <- format(Sys.Date(), "%d-%m-%y")

# Read in the Monaghan tip information
file <- file.path("Desktop", "MONAGHAN_09-12-19", "TipInformation_09-12-19.csv")
monaghanTipInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Read in the Wicklow tip information
file <- file.path("Desktop", "WICKLOW_09-12-19", "TipInformation_09-12-19.csv")
wicklowTipInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

#### Find the FASTQ files ####

# Set the path to the ROI research folder
path <- file.path("storage", "Research", "RepublicOfIreland", "Mbovis")
monaghanPath <- file.path(path, "Monaghan")
wicklowPath <- file.path(path, "Wicklow")

# Get the Monaghan FASTQ file names
fastqDirectories <- getFilesInDirectory(monaghanPath, match="Fastqs", addPath=TRUE)
fastqs <- getFASTQsInDirectories(fastqDirectories)
monaghanFastqs <- keepFASTQsFromPhylogeny(fastqs, monaghanTipInfo)

# Copy the Monaghan FASTQ files
file.copy(monaghanFastqs, file.path("Desktop", "MONAGHAN_09-12-19", "FASTQs"))

# Get the Wicklow FASTQ file names
fastqDirectories <- getFilesInDirectory(wicklowPath, match="Fastqs", addPath=TRUE)
fastqs <- getFASTQsInDirectories(fastqDirectories)
wicklowFastqs <- keepFASTQsFromPhylogeny(fastqs, wicklowTipInfo)

# Copy the Wicklow FASTQ files
file.copy(wicklowFastqs, file.path("Desktop", "WICKLOW_09-12-19", "FASTQs"))

#### Find the FASTA files ####

# Find and copy the Monaghan FASTA and positions files
monaghanFastaFiles <- getFilesInDirectory(file.path(monaghanPath, "vcfFiles"), "fasta", addPath=TRUE)
file.copy(monaghanFastaFiles, file.path("Desktop", "MONAGHAN_09-12-19"))

# Find and copy the Wicklow FASTA and positions files
wicklowFastaFiles <- getFilesInDirectory(file.path(wicklowPath, "vcfFiles"), "fasta", addPath=TRUE)
file.copy(wicklowFastaFiles, file.path("Desktop", "WICKLOW_09-12-19"))

#### Find the phylogeny files ####

# Find, copy, and rename Monaghan phylogeny file
monaghanPhylogeny <- getFilesInDirectory(file.path(monaghanPath, "RAxML_24-09-19"), "RAxML_bipartitions[.]", addPath=TRUE)
file.copy(monaghanPhylogeny, file.path("Desktop", "MONAGHAN_09-12-19"))
copiedPhylogeny <- getFilesInDirectory(file.path("Desktop", "MONAGHAN_09-12-19"), "RAxML_bipartitions[.]", addPath=TRUE)
file.rename(copiedPhylogeny, file.path("Desktop", "MONAGHAN_09-12-19", paste0("mlTree_MONAGHAN_", date, ".tree")))

# Find, copy, and rename Wicklow phylogeny file
wicklowPhylogeny <- getFilesInDirectory(file.path(wicklowPath, "RAxML_21-10-19"), "RAxML_bipartitions[.]", addPath=TRUE)
file.copy(wicklowPhylogeny, file.path("Desktop", "WICKLOW_09-12-19"))
copiedPhylogeny <- getFilesInDirectory(file.path("Desktop", "WICKLOW_09-12-19"), "RAxML_bipartitions[.]", addPath=TRUE)
file.rename(copiedPhylogeny, file.path("Desktop", "WICKLOW_09-12-19", paste0("mlTree_WICKLOW_", date, ".tree")))

#### Check FASTQ file present for each tip ####

# Get all the Monaghan FASTQ files copied and check if all tips present
monaghanFastqs <- getFilesInDirectory(file.path("Desktop", "MONAGHAN_09-12-19", "FASTQs"), ".fastq.gz", addPath=TRUE)
found <- reportWhichTipsAreAbsent(monaghanFastqs, monaghanTipInfo)

# Get all the Wicklow FASTQ files copied and check if all tips present
wicklowFastqs <- getFilesInDirectory(file.path("Desktop", "WICKLOW_09-12-19", "FASTQs"), ".fastq.gz", addPath=TRUE)
found <- reportWhichTipsAreAbsent(wicklowFastqs, wicklowTipInfo)

#### FUNCTIONS ####

reportWhichTipsAreAbsent <- function(fastqs, tipInfo){
  
  # Initialise an array reporting which tips were found
  found <- rep(FALSE, nrow(tipInfo))
  
  # Examine each tip
  for(row in 1:nrow(tipInfo)){
    
    # Find FASTQs
    fastqIndices <- which(grepl(fastqs, pattern=paste0("/", tipInfo[row, 1], "_")))
    
    # Check if need trailing p
    if(length(fastqIndices) == 0){
      fastqIndices <- which(grepl(fastqs, pattern=paste0("/", tipInfo[row, 1], "p_")))
    }
    
    # Check if two FASTQs found
    if(length(fastqIndices) == 2){
      
      found[row] <- TRUE
      
    # Throw message if found more than two
    }else if(length(fastqIndices) > 2){
      warning(paste0("Found ", length(fastqIndices), " associated with the current ID: ", tipInfo[row, 1]))
      print(fastqs[fastqIndices])
      
    # Throw message if not found
    }else{
      warning(paste0("Current ID: ", tipInfo[row, 1], " not found!"))
    }
  }
  
  return(found)
}

keepFASTQsFromPhylogeny <- function(fastqs, tipInfo){
  
  # Initialise an array to store the fastqs to keep
  keep <- c()
  
  # Examine each of the FASTQs
  for(fastq in fastqs){
    
    # Extract the ID from the current fastq
    parts <- strsplit(fastq, split="/")[[1]]
    id <- strsplit(parts[length(parts)], split="_")[[1]][1]
    
    # Remove trailing p is present
    if(grepl(id, pattern="p$")){
      id <- substr(id, 1, nchar(id)-1)
    }
    
    # Check if ID in tipInfo table
    if(id %in% tipInfo[, 1]){
      keep[length(keep) + 1] <- fastq
    }
  }
  
  return(keep)
}

getFASTQsInDirectories <- function(directories){
  
  # Initialise a vector to store the FASTQs
  fastqs <- c()
  
  # Examine each directory
  for(directory in directories){
    
    # Get the fastqs in the current directory
    fastqs <- c(fastqs, file.path(directory, getFilesInDirectory(directory, match="fastq.gz")))
  }
  
  return(fastqs)
}

getFilesInDirectory <- function(directory, match=NULL, addPath=FALSE){
  
  # List the files in the directory
  files <- list.files(path=directory)
  
  # Select files with file ending if provided
  if(is.null(match) == FALSE){
    
    files <- files[grepl(files, pattern=match)]
  }

  # Add path if requested
  if(addPath){
    files <- file.path(directory, files)
  }
  
  return(files)
}