# Reading Spotyping output to get spoligotype codes for Wicklow genomes

#### Read in the data ####

# Set the path
path <- file.path("~", "storage", "Research", "RepublicOfIreland", "Mbovis", "Wicklow", "Spoligotyping")

# Read in the Wicklow FASTQ spoligotypes (created using "RunSpoTypingOnFASTQs_05-11-20.sh")
types <- read.table(file.path(path, "SpotypingOutput_11-05-20.txt"), header=TRUE, sep="\t", colClasses="character")

# Read in the conversion table
conversion <- read.csv(file.path(path, "Mbovis_SpoligotypeCodes_11-05-20.csv"), colClasses="character")
conversion$Binary <- gsub(pattern=" ", replacement="", x=conversion$Binary)

#### Retrieve spoligotypes ####

# Identify the codes for thw Wicklow genomes
types$Spoligotype <- sapply(types$Code, 
                            FUN=function(code, conversion){
                               # Find the index of the current code in the conversion table
                               index <- which(conversion$Binary == code)
                         
                               # Return the spoligotype - if available
                               return(ifelse(length(index) == 0, NA, conversion[index, "SB.Number"]))
                            }, conversion)

# Extract the ID from the FASTQ file names
types$ID <- sapply(types$Files, 
                   FUN=function(fileName){
                     return(strsplit(fileName, split="_")[[1]][1])
                   })

types[, c("ID", "Spoligotype")]
