#### Read in the data ####

# Set the path
path <- "/home/josephcrispell/Desktop/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/"

# Read in the isolate data for the cattle
cattleFile <- paste0(path, "IsolateData/CattleIsolateInfo_AddedNew_TB1484-TB1490_22-03-18.csv")
cattleInfo <- read.table(cattleFile, header=TRUE, stringsAsFactors=FALSE, sep=",")
cattleInfo <- cattleInfo[is.na(cattleInfo$StrainId) == FALSE, ]
cattleInfo$DateCultured <- as.Date(cattleInfo$DateCultured, format="%d/%m/%Y")

# Read in the isolate data for the cattle
badgerFile <- paste0(path, "IsolateData/BadgerInfo_08-04-15_LatLongs_XY_Centroids.csv")
badgerInfo <- read.table(badgerFile, header=TRUE, stringsAsFactors=FALSE, sep=",")
badgerInfo$date <- as.Date(badgerInfo$date, format="%d/%m/%Y")

#### Get the IDs of the badger and cattle isolates used for analyses (193 badgers and 159 cattle) ####

# Get a list of the VCF files
directory <- paste0(path, "vcfFiles/")
vcfFiles <- getFilesInDirectory(directory, ".vcf.gz")

# Get the isolate IDs from the VCF files
ids <- getIDsFromFiles(vcfFiles)
table(grepl(ids, pattern="^WB"))

# Select badger and cattle information only for IDs from VCF files used in analyses
badgerInfo <- badgerInfo[badgerInfo$WB_id %in% ids, ]
cattleInfo <- cattleInfo[cattleInfo$StrainId %in% ids, ]

#### Create the biosample ####

# Read in the template file
templateFile <- paste0(path, "SRA_submission_19-02-19/Biosample_template_19-02-19.tsv")
template <- read.table(templateFile, header=TRUE, check.names=FALSE, skip=11, sep="\t", stringsAsFactors=FALSE)

# Create a table to hold the badger and cattle info
output <- as.data.frame(matrix(NA, nrow=nrow(badgerInfo)+nrow(cattleInfo), ncol=ncol(template)), stringsAsFactors=FALSE)
colnames(output) <- colnames(template)

# Fill the table with the cattle and badger sequence data
output$`*sample_name` <- c(badgerInfo$WB_id, cattleInfo$StrainId)
output$isolate <- c(badgerInfo$WB_id, cattleInfo$StrainId)
output$`*organism` <- rep("Mycobacterium tuberculosis variant bovis", nrow(badgerInfo) + nrow(cattleInfo))
output$`*collection_date` <- c(badgerInfo$date, cattleInfo$DateCultured)
output$`*geo_loc_name` <- rep("United Kingdom: England: Woodchester Park", nrow(badgerInfo) + nrow(cattleInfo))
output$`*sample_type` <- rep("Culture", nrow(badgerInfo) + nrow(cattleInfo))
output$host <- c(rep("BADGER", nrow(badgerInfo)), rep("BOVINE", nrow(cattleInfo)))

# Write the output table to file
outputFile <- paste0(path, "SRA_submission_19-02-19/Metadata_NCBI_19-02-19.tsv")
write.table(output, file=outputFile, quote=FALSE, sep="\t", row.names=FALSE, na="")

#### FUNCTIONS ####

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