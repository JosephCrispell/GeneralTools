#### Read in the Biosample attributes table ####

# Set the path
path <- "/home/josephcrispell/Desktop/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/"

# Read in the biosample attributes table - maps sample IDs to biosample accession numbers
attributesFile <- paste0(path, "SRA_submission_19-02-19/Biosample_Attributes_19-02-19.tsv")
attributes <- read.table(attributesFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)

#### Note the Bioproject number ####

bioProject <- "PRJNA523164"

#### Create the SRA metadata file ####

# Read in the template file
templateFile <- paste0(path, "SRA_submission_19-02-19/SRA_SampleMetadata_template_19-02-19.csv")
template <- read.table(templateFile, header=TRUE, check.names=FALSE, sep=",", stringsAsFactors=FALSE)

# Create a table to hold the sample info
output <- as.data.frame(matrix(NA, nrow=nrow(attributes), ncol=ncol(template)), stringsAsFactors=FALSE)
colnames(output) <- colnames(template)

# Fill in the table with the sample (FASTQ) information
output$bioproject_accession <- rep(bioProject, nrow(output))
output$biosample_accession <- attributes$accession
output$library_ID <- attributes$sample_name
output$title <- rep("WGS for M. bovis from England", nrow(output))
output$library_strategy <- rep("WGS", nrow(output))
output$library_source <- rep("GENOMIC", nrow(output))
output$library_selection <- rep("RANDOM", nrow(output))
output$library_layout <- rep("paired", nrow(output))
output$platform <- rep("ILLUMINA", nrow(output))
output$instrument_model <- ifelse(grepl(attributes$sample_name, pattern="^AF"), "NextSeq 500", "Illumina MiSeq")
output$design_description <- rep("Library preparation available on request", nrow(output))
output$filetype <- rep("fastq", nrow(output))

#### Add in the FASTQ file names ####

# Get an array of the badger FASTQ files and select those used in analyses
badgerFastqs <- c(getFilesInDirectory(paste0(path, "Badger_Batch1/"), "fastq.gz"),
                      getFilesInDirectory(paste0(path, "Badger_Batch2/"), "fastq.gz"))

# Get an array of the cattle FASTQ files and select those used in analyses
cattleFastqs <- c(getFilesInDirectory(paste0(path, "Cattle_FASTQ_1stRound/"), "fastq.gz"),
                      getFilesInDirectory(paste0(path, "Cattle_FASTQ_2ndRound/"), "fastq.gz"),
                      getFilesInDirectory(paste0(path, "NewCattle_16-03-18/fastqs/"), "fastq.gz"))

# Note the FASTQ files associated with each of the IDs that we are interested in
fastqsForEachID <- noteFastqFilesForID(badgerFastqs, cattleFastqs, attributes)

# Add in the FASTQ file names into the output table
for(row in seq_len(nrow(output))){
  output[row, c("filename", "filename2")] <- fastqsForEachID[[output[row, "library_ID"]]]
}

#### Create the output file ####

# Write the output table to file
outputFile <- paste0(path, "SRA_submission_19-02-19/SampleInfo_NCBI_19-02-19.tsv")
write.table(output, file=outputFile, quote=FALSE, sep="\t", row.names=FALSE, na="")

#### Copy the FASTQ files into a single directory ####

# Create a directory for the files
directory <- paste0(path, "SRA_submission_19-02-19/FASTQs/")
dir.create(directory)

# Copy the badger FASTQ files into the folder - only those that were used in analysis
copyFASTQsInDirectoryToAnother(path, from="Badger_Batch1/", to=directory, attributes)
copyFASTQsInDirectoryToAnother(path, from="Badger_Batch2/", to=directory, attributes)

# Copy the cattle FASTQ files into the folder - only those that were used in analysis
copyFASTQsInDirectoryToAnother(path, from="Cattle_FASTQ_1stRound/", to=directory, attributes)
copyFASTQsInDirectoryToAnother(path, from="Cattle_FASTQ_2ndRound/", to=directory, attributes)
copyFASTQsInDirectoryToAnother(path, from="NewCattle_16-03-18/fastqs/", to=directory, attributes)

#### FUNCTIONS ####

copyFASTQsInDirectoryToAnother <- function(path, from, to, attributes){
  
  # Get a list of the FASTQ files in the from directory
  fastqs <- getFilesInDirectory(paste0(path, from), "fastq.gz")
  
  # Get the IDs from these FASTQ files
  ids <- getIDsFromFiles(fastqs)
  
  # Copy only the used FASTQs into to directory
  file.copy(paste0(path, from, fastqs[ids %in% attributes$sample_name]), to)
}

noteFastqFilesForID <- function(badgerFastqs, cattleFastqs, attributes){
  
  # Initialise a list to store the pair of FASTQs associated with each ID
  fastqsForEachID <- list()
  
  # Examine the FASTQ files
  for(fastq in c(badgerFastqs, cattleFastqs)){
    
    # Get the ID from the current FASTQ
    id <- strsplit(fastq, split="_")[[1]][1]
    
    # Check that ID exists in attributes table
    if(id %in% attributes$sample_name == FALSE){
      next
    }
    
    # Check if we have encountered the current ID
    if(is.null(fastqsForEachID[[id]])){
      fastqsForEachID[[id]] <- c(fastq)
    }else{
      fastqsForEachID[[id]] <- c(fastqsForEachID[[id]], fastq)
    }
  }
  
  return(fastqsForEachID)
}

getIDsFromFiles <- function(files){
  
  # Create an array to store the IDs
  ids <- c()
  
  # Examine each of the forward fastq files
  for(file in files){
    
    # Get the ID from the current file name
    ids[length(ids) + 1] <- strsplit(file, split="_")[[1]][1]
  }
  
  return(ids)
}

getFilesInDirectory <- function(directory, fileEnding=NULL){
  
  # List the files in the directory
  files <- list.files(path=directory)
  
  # Select files with file ending if provided
  if(is.null(fileEnding) == FALSE){
    
    files <- files[grepl(files, pattern=paste0(fileEnding, "$"))]
  }
  
  return(files)
}