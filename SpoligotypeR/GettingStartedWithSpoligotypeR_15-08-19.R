#### libraries ####

library(R.utils)

#### Preparation ####

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

# Create a path variable
path <- "/home/josephcrispell/Desktop/Research/spoligotypeR/"

#### Build local blastn database ####

# Note the path to blast
pathToBlast="/home/josephcrispell/Desktop/Research/ncbi-blast-2.6.0+/bin/"

# Build the local database
createLocalBlastDatabase(paste0(path, "Spacers/spacer.fasta"), pathToBlast)

#### Dealing with VCF file ####

# Open the VCF file
vcfFile <- paste0(path, "TestingData/vcfFiles/14-MBovis_1.vcf.gz")
vcf <- readVCF(vcfFile)

# Create a FASTA file from vcf
fastaFile <- paste0(vcfFile, ".fasta")
parts <- strsplit(vcfFile, split="/")[[1]]
sequenceName <- parts[length(parts)]
convertVCFToFASTA(vcf, fastaFile, seqName=sequenceName, depthThreshold=20, supportThreshold=0.95)

#### Dealing with FASTQ file ####

#### Dealing with FASTA file ####

#### Issues so far ####

# Unzipping is ubuntu specific

# Need full path to blastn

# Multithreading would be useful!

#### FUNCTIONS - general ####

gzip <- function(fileName, unzip=FALSE){
  
  # Check if zipping or unzipping
  if(unzip){
    
    # Unzip the file
    R.utils::gunzip(fileName)
    
    # Note the new file name
    fileName <- substr(fileName, 1, nchar(fileName)-3)
  }else{
    
    # Zip up the file
    R.utils::gzip(fileName)
    
    # Note the fileName
    fileName <- paste0(fileName, ".gz")
  }
  
  return(fileName)
}

#### FUNCTIONS - VCF ####

convertVCFToFASTA <- function(vcf, fastaFile, seqName, depthThreshold=20, supportThreshold=0.95){
  
  # Create a vector of nucleotides from the vcf
  nucleotides <- c()
  nucleotides[vcf$POS] <- apply(vcf, MARGIN=1, FUN=getNucleotideFromVCFRow, depthThreshold, supportThreshold)
  
  # Fill in the gaps
  nucleotides[is.na(nucleotides)] <- "-"
  
  # Print the nucleotide sequence to file
  fileConnection <- file(fastaFile)
  writeLines(c(paste0(">", seqName), paste(nucleotides, collapse="")), fileConnection)
  close(fileConnection)
}

getNucleotideFromVCFRow <- function(row, depthThreshold, supportThreshold){

  # Initialise a character for the allele called
  allele <- 'N'
  
  # Get values for the metrics of interest from the INFO column
  infoValues <- getINFO(row[8], metrics=c("DP", "DP4"))
  
  # Skip if didn't find metrics
  if(is.null(infoValues[["DP"]]) || is.null(infoValues[["DP4"]])){
    return(allele)
  }
  
  # Calculate number of high quality reads
  nHighQualityReads <- sum(infoValues$DP4)
  
  # Skip if no high quality reads available
  if(nHighQualityReads == 0){
    return(allele)
  }
  
  # Calculate allele support
  support <- c(
    (infoValues$DP4[1]+infoValues$DP4[2])/nHighQualityReads,
    (infoValues$DP4[3]+infoValues$DP4[4])/nHighQualityReads)
  
  # Check if passed filters
  if(infoValues$DP >= depthThreshold && max(support) >= supportThreshold){
    
    # Call the allele
    allele <- c(row[4], row[5])[which.max(support)]
  }
  
  return(allele)
}

getINFO <- function(infoColumn, metrics){
  
  # Create regular expression to quickly search for metrics
  metricMatch <- paste(metrics, collapse="|")
  
  # Split the INFO column into its columns
  columns <- strsplit(infoColumn, split=";")[[1]]
  
  # Initialise a list to store the values
  infoValues <- list()
  valuesFound <- 0
  
  # Examine each part fo the INFO column
  for(column in columns){
    
    # Check if found metric
    if(grepl(column, pattern=metricMatch)){
      
      # Split the current metric into name and value
      parts <- strsplit(column, split="=")[[1]]
      
      # Store the metric and its value(s)
      infoValues[[parts[1]]] <- as.numeric(strsplit(parts[2], split=",")[[1]])
      
      # Note that found value
      valuesFound <- valuesFound + 1
    }
    
    # Check if already found all values wanted
    if(valuesFound == length(metrics)){
      break
    }
  }
  
  return(infoValues)
}

readVCF <- function(fileName){
  
  # Create a variabel to remember whether zipped
  zipped <- FALSE
  
  # Check if zipped
  if(grepl(fileName, pattern=".gz$")){
    
    # Note that VCF file was zipped
    zipped <- TRUE
    
    # unzip the VCF file
    fileName <- gzip(fileName, unzip=TRUE)
  }
  
  # Read in the vcf file
  # Header automatically ignore as starts with comment
  vcf <- read.table(fileName, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                    col.names=c("CHROM","POS","ID", "REF", "ALT","QUAL", "FILTER", "INFO", "FORMAT", "FORMAT_values"))
  
  # If input file was zipped - zip it back up
  if(zipped){
    fileName <- gzip(fileName)
  }
  
  return(vcf)
}

#### FUNCTIONS - BLAST ####

createLocalBlastDatabase <- function(spacerFile, pathToBlast){
  
  # Check if local database already exists
  if(file.exists(paste0(spacerFile, ".nhr")) == FALSE){
    
    # Build the command to build the database
    command <- paste0(pathToBlast, "makeblastdb -in ", spacerFile, " -parse_seqids -dbtype \"nucl\" -out ", spacerFile)
    
    # Build the local database
    system(command, wait=TRUE)
  }
}