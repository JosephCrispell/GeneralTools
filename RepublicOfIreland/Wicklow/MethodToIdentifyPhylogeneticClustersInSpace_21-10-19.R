#### Preparation ####

# Set the path
path <- "/home/josephcrispell/Desktop/Research/RepublicOfIreland/Mbovis/Wicklow/PhylogeneticAndSpatialClusters/"

# Get the current date
today <- format(Sys.Date(), "%d-%m-%y")

#### Read in genetic and spatial distance matrices ####

# Read in the genetic distance matrix
geneticDistancesFile <- file.path(path, "GeneticDistanceMatrix_21-10-19.tsv")
geneticDistances <- read.table(geneticDistancesFile, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)

# Read in the spatial distance matrix
spatialDistancesFile <- file.path(path, "SpatialDistanceMatrix_21-10-19.tsv")
spatialDistances <- read.table(geneticDistancesFile, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)

#### Identify clusters on phylogeny and in space ####

# Set the SNP distance to define cluster
snpThreshold <- 3


#### FUNCTIONS ####

getClustersFromMatrixWithThreshold <- function(matrix, threshold){
  
  # Create a list to note the cluster each sequence has been assigned
  clusters <- list()
  cluster <- 0
  
  # Get the row and column namaes
  rowNames <- rownames(matrix)
  colNames <- colnames(matrix)
  
  # Examine each cell in the matrix
  for(i in seq_along(matrix)){
    
    # Get the current row ID
    row <- rowNames[i]
    
    for(j in seq_along(matrix)){
      
      # Skip diagonal, lower part of matrix and any values aboce threshold
      if(i >= j || matrix[i, j] > threshold){
        next
      }
      
      # Get the current column ID
      col <- colNames[j]
      
      # Check if row has already been assigned to a cluster
      if(is.null(clusters[[row]]) == FALSE && is.null(clusters[[col]]) == TRUE){
        
        # Assign column to same cluster
        clusters[[col]] <- clusters[[row]]
        
      # Check if column has already been assigned to a cluster
      }else if(is.null(clusters[[row]]) == TRUE && is.null(clusters[[col]]) == FALSE){
        
        # Assign row to same cluster
        clusters[[col]] <- clusters[[row]]
        
      # Check if neither row or column have been assigned to a cluster
      }else if(is.null(clusters[[row]]) == TRUE && is.null(clusters[[col]]) == TRUE){
        
        # Create a new cluster
        cluster <- cluster + 1
        
        # Put both row and olumn into new cluster
        clusters[[row]] <- cluster
        clusters[[col]] <- cluster
        
      # Check if have been assigned to different clusters!
      }else if(clusters[[row]] != clusters[[col]]){
        stop("Clustered row (", row, ":", clusters[[row]], ") and column (", col, ":", clusters[[col]],
             ") have been assigned to different clusters!")
      }
    }
  }
  
  # Initialise an output data.frame recording each sequences cluster
  clusterTable <- data.frame("SequenceID"=names(matrix), "Cluster"=NA)
  
  # Retrieve the cluster assignments from the clusters list
  for(sequenceID in names(clusters)){
    
    # Find row in cluster table
    row <- which(clusterTable$SequenceID == sequenceID)
    
    # Assign the cluster
    clusterTable[row, "Cluster"] <- clusters[[sequenceID]]
  }
  
  return(clusterTable)
}