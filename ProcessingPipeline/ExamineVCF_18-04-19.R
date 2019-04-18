#### Load packages ####

library(pegas)

#### Load the vcf file ####

# Set the path variable
path <- "/home/josephcrispell/Desktop/Research/RepublicOfIreland/Fastqs_ALL_15-03-19/Alicia/"

# Read in the vcf file
vcfFile <- paste0(path, "variants.vcf")
vcf <- read.table(vcfFile, header=TRUE, stringsAsFactors=FALSE, sep="\t", comment.char="@", check.names=FALSE, skip=28)
#vcf <- VCFloci(vcfFile)


#### Examine the coverage ####

# Get the read depth for each site in vcf
#depth <- getINFO(vcf, what="DP")
depth <- apply(vcf, 1, getReadDepth)

# Plot the depth distribution
hist(depth, las=1, xlab="Read Depth")
legend("topright", legend=paste("mean =", round(mean(depth, na.rm=TRUE), digits=2)), bty="n")

#### Identify variant sites ####

# Identify loci with alternate allele
variants <- vcf[vcf$ALT != ".", ]

# Get the INFO entries for each variant
variants <- splitInfoColumn(variants)

#### FUNCTIONS ####

splitInfoColumn <- function(variants){
  
  # Examine the INFO column for each variant
  for(row in seq_len(nrow(variants))){
    
    # Split the INFO column for the current loci into its parts
    entries <- strsplit(variants[row, "INFO"], split=";")[[1]]
    
    # Examine each of the INFO column entries
    for(entry in entries){
      
      # Split the current entry into its ID and value
      parts <- strsplit(entry, split="=")[[1]]
      
      # Check if INFO metric exists in table column names
      colIndex <- which(colnames(variants) == parts[1])
      if(length(colIndex) == 1){
        
        variants[row, colIndex] <- parts[2]
      }else{
        variants[, parts[1]] <- 1
        variants[row, parts[1]] <- parts[2]
      }
    }
  }
  
  return(variants)
}

getReadDepth <- function(vcf){
  
  # Split the current loci's INFO column into its individual elements and extract DEPTH
  return(as.numeric(strsplit(strsplit(vcf["INFO"], split=";")[[1]][1], split="=")[[1]][2]))
}
