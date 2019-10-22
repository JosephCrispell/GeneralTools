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
spatialDistances <- read.table(spatialDistancesFile, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)

#### Identify clusters on phylogeny and in space ####

# Define genetic distance clusters based on threshold
geneticClusters <- getClustersFromMatrixWithThreshold(geneticDistances, threshold=2)

# Define sub-clusters in the genetic clusters based on spatial distances
clusters <- examineSpatialDistancesInEachCluster(spatialDistances, geneticClusters, threshold=2500)

# Highlight the sub-clusters on the linked phylogeny and spatial locations plot
# Check this link for solution: https://stackoverflow.com/questions/41050532/is-there-any-way-to-draw-a-boundary-around-a-group-of-points-in-r

#### FUNCTIONS ####

getSpatialDistancesForCluster <- function(isolates, distances){
  
  # Initialise a matrix to store the distances associated with the current cluster
  matrix <- matrix(NA, nrow=length(isolates), ncol=length(isolates))
  rownames(matrix) <- isolates
  colnames(matrix) <- isolates
  
  # Fill in the distances for the matrix
  for(i in seq_along(isolates)){
    for(j in seq_along(isolates)){
      
      # Skip upper triangle
      if(i >= j){
        next
      }
      
      # Store the current distance
      matrix[i, j] <- distances[isolates[i], isolates[j]]
      matrix[j, i] <- distances[isolates[i], isolates[j]]
    }
  }
  
  return(matrix)
}

examineSpatialDistancesInEachCluster <- function(distances, clusters, threshold){
  
  # Add a sub cluster column to the clusters table
  clusters$SubSpatialCluster <- NA
  
  # Examine each cluster
  for(cluster in unique(clusters$Cluster)){
    
    # Get the isolate IDs for the current cluster
    isolates <- clusters[clusters$Cluster == cluster, "SequenceID"]
    
    # Get a spatial distance matrix for the current cluster
    matrix <- getSpatialDistancesForCluster(isolates, distances)
    
    # Identify any sub-clusters based upon spatial distances
    subClusters <- getClustersFromMatrixWithThreshold(matrix, threshold)
    
    # Assign the subclusters
    for(row in seq_len(nrow(subClusters))){
      clusters$SubSpatialCluster[which(clusters$SequenceID == subClusters[row, "SequenceID"])] <- subClusters[row, "Cluster"]
    }
  }
  
  return(clusters)
}

mergeClusters <- function(clusters, a, b){
  
  # Examine every entry in the clusters list
  for(id in names(clusters)){
    
    # Check if current isolate was assigned cluster b
    if(clusters[[id]] == b){
      
      # Reassign current isolate to cluster a
      clusters[[id]] <- a
    }
  }
  
  return(clusters)
}

getClustersFromMatrixWithThreshold <- function(matrix, threshold){
  
  # Create a list to note the cluster each sequence has been assigned
  clusters <- list()
  cluster <- 0
  
  # Get the row and column namaes
  rowNames <- rownames(matrix)
  colNames <- colnames(matrix)
  
  # Examine each cell in the matrix
  for(i in seq_along(rowNames)){
    
    # Get the current row ID
    row <- rowNames[i]
    
    for(j in seq_along(colNames)){
      
      # Skip diagonal, lower part of matrix and any values aboce threshold
      if(i >= j || is.na(matrix[i, j]) || matrix[i, j] > threshold){
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
        clusters[[row]] <- clusters[[col]]

      # Check if neither row or column have been assigned to a cluster
      }else if(is.null(clusters[[row]]) == TRUE && is.null(clusters[[col]]) == TRUE){

        # Create a new cluster
        cluster <- cluster + 1

        # Put both row and olumn into new cluster
        clusters[[row]] <- cluster
        clusters[[col]] <- cluster

      # Check if have been assigned to different clusters!
      }else if(clusters[[row]] != clusters[[col]]){

        # Merge the individuals into one of the clusters
        clusters <- mergeClusters(clusters, clusters[[row]], clusters[[col]])
      }
    }
  }
  
  # Initialise an output data.frame recording each sequences cluster
  clusterTable <- data.frame("SequenceID"=rownames(matrix), "Cluster"=NA, stringsAsFactors=FALSE)
  
  # Retrieve the cluster assignments from the clusters list
  for(sequenceID in names(clusters)){
    
    # Find row in cluster table
    row <- which(clusterTable$SequenceID == sequenceID)
    
    # Assign the cluster
    clusterTable[row, "Cluster"] <- clusters[[sequenceID]]
  }
  
  # Remove NA rows
  clusterTable <- clusterTable[is.na(clusterTable$Cluster) == FALSE, ]
  
  return(clusterTable)
}
