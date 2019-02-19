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
output$sample_name <- attributes$accession
output$library_ID <- attributes$sample_name
output$title <- rep("WGS for M. bovis from England", nrow(output))
output$library_strategy <- rep("WGS", nrow(output))
output$library_source <- rep("GENOMIC", nrow(output))
output$library_layout <- rep("paired", nrow(output))
output$platform <- rep("ILLUMINA", nrow(output))
output$instrument_model <- ifelse(grepl(attributes$sample_name, pattern="^WB"), "Illumina MiSeq", "Illumina HiSeq")
output$filetype <- rep("fastq", nrow(output))

#### Add in the FASTQ file names and MD5 sums ####

# Get an array of the badger FASTQ files and select those used in analyses

# Get an array of the cattle FASTQ files and select those used in analyses

# md5sums?
