#### Load libraries ####

library(ape) # Reading and phylogeny
library(rgdal) # Convert X and Y to lat longs and reading in shape files
library(basicPlotteR) # Add SNV scale
library(grid) # Used to plot lines between plot panels
library(raster) # Binding two shape files
library(OpenStreetMap) # Great tutorial here: https://www.r-bloggers.com/the-openstreetmap-package-opens-up/
library(geiger) # For the tips function
library(phyloHelpeR) # Building phylogeny with RAxML

#### Read in the sample data ####

# Set the path 
path <- "J:\\WGS_Monaghan\\"
#path <- "/home/josephcrispell/storage/Research/RepublicOfIreland/Mbovis/Monaghan/"

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

# Read in the animal IDs
file <- paste0(path, "EartagsAndSettIDs_17-12-19.csv")
animalTags <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Read in the cattle herd IDs
file <- paste0(path, "Animal_HerdIDs_17-12-19.csv")
herdIDs <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Read in the FASTA file
fastaFile <- paste0(path, "vcfFiles/sequences_Prox-10_17-12-2019.fasta")
sequences <- read.dna(fastaFile, as.character=TRUE, format="fasta")
nSites <- ncol(sequences)

# Read in the coverage information
coverageFile <- paste0(path, "vcfFiles/isolateCoverageSummary_DP-20_17-12-2019.txt")
coverage <- read.table(coverageFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Read in the test summaries tables
file <- paste0(path, "AHCS_data/Batch1_TestingInformation_20190731.csv")
cattleTestInfoBatch1 <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE,
                                   check.names=FALSE)
file <- paste0(path, "AHCS_data/Batch2_TestingInformation_20191119.csv")
cattleTestInfoBatch2 <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE,
                                   check.names=FALSE)

# Read in the badger capture information
file <- paste0(path, "Badger_data/TBL_captured_badgers_25_11_2019.txt")
badgerCaptureData <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Read in the sett capture event information
file <- paste0(path, "Badger_data/TBL_capture_events_25_11_2019.txt")
settCaptureEventData <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)


#### Build the phylogeny ####

# Build a phylogeny using RAxML
tree <- runRAXML(fastaFile, date="17-12-19", path, outgroup="\\>Ref-1997", alreadyRun=TRUE)

# Remove Reference
tree <- drop.tip(tree, tree$tip.label[grepl(tree$tip.label, pattern=">Ref-1997")])

# Edit the tip labels
tree$tip.label <- editTipLabels(tree$tip.label)

# Get the tip information (species and sampling date)
tipInfo <- getTipInfo(tree$tip.label, animalTags, coverage, herdIDs)

# Convert branch lengths to SNPs
tree$edge.length <- tree$edge.length * nSites

# Add in slaughter and capture dates
tipInfo <- addSlaughterOrCaptureDates(tipInfo, cattleTestInfoBatch1, cattleTestInfoBatch2,
                                      badgerCaptureData, settCaptureEventData)

# Note badgers not able to find capture data for
missing <- tipInfo[is.na(tipInfo$HerdOrSettID) & tipInfo$Species == "BADGER", ]
missing <- missing[is.na(missing$Tip) == FALSE, ]
write.table(missing[, c("Aliquot", "AnimalID")], sep=",", quote=FALSE, row.names=FALSE,
            file=paste0(path, "MissingBadgerTags_", date, ".csv"))

#### EXAMINE THE SEQUENCE QUALITY!!! (REMOVE UNINFORMATIVE SITES!) ####

tipInfo[tipInfo$Coverage < 0.9, c("Tip", "Aliquot", "Coverage", "Species")]


#### Plot the phylogeny ####

# Set the tip shapes and colours
tipShapesAndColours <- list("BADGER"=c(rgb(1,0,0,1), 19),
                            "COW"=c(rgb(0,0,1,1), 17),
                            "NA"=c("grey", 15))

# Open a PDF
pdf(paste0(path, "Figures/MbovisAnnotatedPhylogeny_", date, ".pdf"), height=14)

plotPhylogeny(tree, tipInfo, tipShapesAndColours=tipShapesAndColours,
              addBranchColours=FALSE)



# Close the pdf
dev.off()

pdf(paste0(path, "Figures/MbovisAnnotatedPhylogeny_FAN_", date, ".pdf"), width=14)

plotPhylogeny(tree, tipInfo, tipShapesAndColours=tipShapesAndColours,
              addBranchColours=FALSE, type="fan", scaleSize=50)

# Close the pdf
dev.off()



# Look at the node numbers
pdf("test.pdf", width=50, height=20)
plot.phylo(tree, show.tip.label=FALSE)
nodelabels(frame="none")
dev.off()


#### Plot the temporal sampling ####

# Open a PDF
pdf(paste0(path, "Figures/TemporalSampling_", date, ".pdf"), width=14)

plotTemporalSamplingRange(tipInfo)

# Close the pdf
dev.off()

#### Read in the spatial data ####

# Get the badger sett locations
shapeFile <- paste0(path, "Badger_data/badgers.shp")
badgerShapeFile <- readOGR(dsn=shapeFile)
settLocations <- data.frame("SettID"=badgerShapeFile@data$SETT_NO,
                            "X"=badgerShapeFile@coords[, 1],
                            "Y"=badgerShapeFile@coords[, 2])

# Read in the cattle land parcels
shapeFile <- paste0(path, "Cattle_data/herds_mn_july_2019.shp")
cattleShapeFileJuly <- readOGR(dsn = shapeFile)
shapeFile <- paste0(path, "Cattle_data/herds_mn_november_2019.shp")
cattleShapeFileNovember <- readOGR(dsn = shapeFile)
cattleShapeFile <- bind(cattleShapeFileJuly, cattleShapeFileNovember)

# Extract the polygon coords
landParcelCoordsJuly <- getPolygonCoords(cattleShapeFileJuly)
landParcelCoordsNovember <- getPolygonCoords(cattleShapeFileNovember)
landParcelCoords <- c(landParcelCoordsJuly, landParcelCoordsNovember)

# Calculate the herd centroids
landParcelCentroids <- calculateLandParcelCentroids(landParcelCoords)

# Score the herds based upon distances to land parcels
landParcelCentroids <- calculateDistancesToOverallCentroidForEachHerd(landParcelCentroids,
                                                                      tipInfo)

#### Plot the sampling locations ####

# Set the shapes and colours
tipShapesAndColours <- list("BADGER"=c("red", 19, 21),
                            "COW"=c("blue", 17, 24),
                            "NA"=c("grey", 15, 15))

# Add herd centroids and sett coordinates into tip information table
tipInfo <- addHerdCentroidsAndSettLocationsToTipInfo(tipInfo, settLocations,
                                                     landParcelCentroids)

# Open a PDF
outputPlotFile <- paste0(path, "Figures\\SamplingLocations_", date, ".pdf")
pdf(outputPlotFile)

# Plot the sampling locations
plotSamplingLocations(tipInfo, tipShapesAndColours)

# Plot the sampling locations with land parcels
plotSamplingLocations(tipInfo, tipShapesAndColours, plotPolygons=TRUE,
                      cattleShapeFile, landParcelCoords, landParcelCentroids)

# Close the PDF
dev.off()

#### Plot the phylogeny and sampling locations together ####

# Construct spatial and genetic distance matrices
geneticDistanceMatrix <- calculateGeneticDistances(sequences, tree)
spatialDistanceMatrix <- calculateSpatialDistances(tipInfo)

# Identify spatial clusters based upon genetic and spatial distance thresholds
geneticClusters <- getClustersFromMatrixWithThreshold(geneticDistanceMatrix, threshold=3)
clusters <- examineSpatialDistancesInEachCluster(spatialDistanceMatrix, geneticClusters,
                                                 threshold=2500)

# Set the rownames of the tip information dataframe
rownames(tipInfo) <- tipInfo$Tip

# Set the tip shapes and colours
tipShapesAndColours <- list("BADGER"=c(rgb(1,0,0,1), 19),
                            "COW"=c(rgb(0,0,1,1), 17),
                            "NA"=c("grey", 15))

# Open a pdf
file <- paste0(path, "Figures/", "PhylogenyAndLocations_", date, ".pdf")
pdf(file, width=20, height=14)

# Plot the phylogeny linked to sampling locations
plotPhylogenyAndLocations(tree, tipInfo, tipShapesAndColours,
                          tipShapeCexOnPhylogeny=3, scaleCex=3,
                          addTipIndices=TRUE,
                          layoutMatrix=matrix(c(1,1,2,2,2), nrow=1, ncol=5, byrow=TRUE),
                          tipLabelOffset=10, snpScaleSize=50,
                          spatialClusters=clusters)

# Close the PDF
dev.off()


#### DATA I NEED ####

# Update:
#   - Cattle herd land parcels - Guy
#   - Badger sett information - Guy
#   - Sampling information for new genomes - Kevin
#   - Sampled cattle testing histories - Dan


#### FUNCTIONS - plot sampling locations ####

plotSamplingLocations <- function(tipInfo, tipShapesAndColours,
                                  plotPolygons=FALSE,
                                  cattleShapeFile=NULL, landParcelCoords=NULL,
                                  herdCentroids=NULL){
  
  
  # Note the sampled herds
  sampledHerds <- getSampledHerds(tipInfo)
  
  # Create an empty plot
  par(mar=c(0,0,0,0))
  plot(tipInfo$X, tipInfo$Y, col="white", xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
  
  # Plot the sampled herd land parcels
  if(plotPolygons){
    plotHerdLandParcels(sampledHerds, cattleShapeFile, landParcelCoords,
                        herdCentroids)
  }
  
  # Plot the sampling locations
  points(tipInfo[, c("X", "Y")], 
         pch=getTipShapeOrColourBasedOnSpecies(tipInfo, tipShapesAndColours,
                                               which="shape", shapeIndex=3),
         bg=getTipShapeOrColourBasedOnSpecies(tipInfo, tipShapesAndColours,
                                              which="colour", alpha=1),
         cex=2, bty="n", xaxt="n", yaxt="n",
         xlab="", ylab="")
  
  # Add a scale
  addScale(size=10000)
  
  # Add a legend
  legend("topright",
         legend=c("Badger", "Cow"), pch=c(21, 24), text.col=c("red", "blue"),
         pt.bg=c("red", "blue"), pt.cex=3, xpd=TRUE, horiz=FALSE, bty="n", cex=2)
  
  # Add a legend for the cattle herd location certainty
  range <- range(tipInfo$MeanDistanceToLandParcels, na.rm=TRUE)
  legendValues <- c(500, 1000, 5000, 10000)
  legend("topleft", legend=legendValues/1000, pch=24, 
         pt.bg=c(rgb(0,0,1, 1 - legendValues[1]/range[2]),
                 rgb(0,0,1, 1 - legendValues[2]/range[2]),
                 rgb(0,0,1, 1 - legendValues[3]/range[2]),
                 rgb(0,0,1, 1 - 1)),
         bty="n", pt.cex=2, col="black",
         title="Mean distance to land parcels (km)", xpd=TRUE)
}

plotHerdLandParcels <- function(sampledHerds, cattleShapeFile, landParcelCoords,
                                herdCentroids){
  
  # Plot the land parcels
  for(herd in sampledHerds){
    
    # Get the polygon set ID for the current herd
    polygonSetID <- which(cattleShapeFile@data$HERD_NO == herd)[1]
    fieldIndex <- 0
    
    # Plot each of the current herds land parcels
    for(coords in landParcelCoords[[polygonSetID]]){
      
      # Increment the field index
      fieldIndex <- fieldIndex + 1
      
      # Plot the current field
      polygon(coords)
      
      # Plot line from current field to centroid of all fields
      points(x=c(herdCentroids[[herd]]$X, 
                 herdCentroids[[herd]]$CentroidXs[fieldIndex]),
             y=c(herdCentroids[[herd]]$Y, 
                 herdCentroids[[herd]]$CentroidYs[fieldIndex]),
             type="l", col=rgb(0,0,1, 0.5))
    }
  }
}

#### FUNCTIONS - Identify genetic and spatial clusters ####

plotShapeAroundSpatialClusters <- function(tipInfo, clusters,
                                           expand=0.1,
                                           lty=2, lwd=2, 
                                           border="blue",
                                           fill=rgb(0,0,1, 0.1),
                                           pch=23, cex=1, threshold=2500){
  
  # Taken from this answer: https://stackoverflow.com/questions/41050532/is-there-any-way-to-draw-a-boundary-around-a-group-of-points-in-r
  
  # Examine each of the clusters
  for(clusterID in names(clusters)){
    
    # Get coordinates of points in cluster
    clusterInfo <- tipInfo[tipInfo$Tip %in% clusters[[clusterID]], ]
    
    # Get the X and Y coordinates
    x <- clusterInfo$X
    y <- clusterInfo$Y
    
    # Keep only the unique coordinates
    keep <- c()
    uniqueCoords <- c()
    for(i in seq_along(x)){
      
      # Create a location key
      coords <- paste0(x[i], ":", y[i])
      
      # Check if already noted this location
      if(coords %in% uniqueCoords){
        
      # Otherwise store it and its index
      }else{
        uniqueCoords[length(uniqueCoords) + 1] <- coords
        keep[length(keep) + 1] <- i
      }
    }
    x <- x[keep]
    y <- y[keep]
    
    # Check if only two coordinates
    if(length(x) == 2){

      # Find the middle between the two coordinates
      x <- x[1] + (0.5*(x[2] - x[1]))
      y <- y[1] + (0.5*(y[2] - y[1]))
      
      # Create four new coordinates
      x <- c(x-(0.5*threshold), x, x+(0.5*threshold), x)
      y <- c(y, y+(0.5*threshold), y, y-(0.5*threshold))
      
      #points(x, y, type="l", lty=lty, lwd=lwd, col=border)
      #next
      
    # Check if only one coordinate
    }else if(length(x) == 1){

      # Create four new coordinates
      x <- c(x[1]-(0.5*threshold), x[1], x[1]+(0.5*threshold), x[1])
      y <- c(y[1], y[1]+(0.5*threshold), y[1], y[1]-(0.5*threshold))
      
      #points(x, y, pch=pch, cex=cex, col=border, bg=fill)
      #next
    }

    # Calculate the centre of random points
    meanX = mean(x)
    meanY = mean(y)
    
    # Create convex hull around the random points
    convexHullIndices = chull(data.frame(x, y))
    convexHullXCoords <- x[convexHullIndices]
    convexHullYCoords <- y[convexHullIndices]
    
    # Slightly expand convex hull away from centre
    expandedConvexHullXCoords <- convexHullXCoords + 
      (expand*(convexHullXCoords - meanX))
    expandedConvexHullYCoords <- convexHullYCoords + 
      (expand*(convexHullYCoords - meanY))

    # Plot a polygon for the expanded convex hull
    polygon(expandedConvexHullXCoords, expandedConvexHullYCoords,
            lty=lty, lwd=lwd, border=border, col=fill)
  }
}

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
  
  # Note the overall cluster assignment
  clusters$OverallCluster <- ifelse(is.na(clusters$SubSpatialCluster), NA,
                                    paste0(clusters$Cluster, "-", clusters$SubSpatialCluster))
  
  # Return a list containing the isolates assigned to each cluster
  output <- list()
  for(row in seq_len(nrow(clusters))){
    
    # Skip NA rows
    if(is.na(clusters[row, "OverallCluster"])){
      next
    }
    
    # Check if found cluster before
    if(is.null(clusters[row, "OverallCluster"]) == FALSE){
      output[[clusters[row, "OverallCluster"]]] <- c(output[[clusters[row, "OverallCluster"]]],
                                                     clusters[row, "SequenceID"])
    }else{
      output[[clusters[row, "OverallCluster"]]] <- c(clusters[row, "SequenceID"])
    }
  }
  
  return(output)
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

#### FUNCTIONS - distances ####

calculateGeneticDistances <- function(sequences, tree){
  
  # Calculate the genetic distances between all sequences
  distances <- as.matrix(dist.dna(as.DNAbin(sequences), model="N"))
  
  # Parse the sequence names
  colnames(distances) <- editTipLabels(colnames(distances))
  rownames(distances) <- colnames(distances)
  
  # Select the distances for the current tree
  return(distances[tree$tip.label, tree$tip.label])
}

calculateSpatialDistances <- function(tipInfo){
  
  # Initialise a matrix to store the distances
  distances <- matrix(NA, nrow=nrow(tipInfo), ncol=nrow(tipInfo))
  colnames(distances) <- tipInfo$Tip
  rownames(distances) <- tipInfo$Tip
  
  # Compare each of the sampling locations
  for(i in seq_len(nrow(tipInfo))){
    
    # Skip current i if no location available
    if(is.na(tipInfo[i, "X"]) || is.na(tipInfo[i, "Y"])){
      next
    }
    
    for(j in seq_len(nrow(tipInfo))){
      
      # Skip current j if no location available
      if(is.na(tipInfo[j, "X"]) || is.na(tipInfo[j, "Y"])){
        next
      }
      
      # Skip self comparisons and making same comparison twice
      if(i >= j){
        next
      }
      
      # Calculate the distance between the current pair of sampling dates
      distances[i, j] <- euclideanDistance(x1=tipInfo[i, "X"], y1=tipInfo[i, "Y"],
                                           x2=tipInfo[j, "X"], y2=tipInfo[j, "Y"])
      distances[j, i] <- distances[i, j]
    }
  }
  
  return(distances)
}

euclideanDistance <- function(x1, y1, x2, y2){
  return(sqrt((x1 - x2)^2 + (y1 - y2)^2))
}

#### FUNCTIONS - joint phylogeny and spatial data ####

plotPhylogenyAndLocations <- function(tree, tipInfo, tipShapesAndColours,
                                      tipShapeCexOnPhylogeny=2, scaleCex=1.5,
                                      scaleTextColour="black", addTipIndices=FALSE,
                                      connectingLinesWidth=1, connectingLinesAlpha=0.1,
                                      tipIndexBackground=rgb(0,0,0, 0.1),
                                      tipCexOnMap=2, 
                                      layoutMatrix=NULL, tipLabelOffset=1,
                                      tipLabelCexOnPhylogeny=2, showTipLabels=FALSE,
                                      spatialClusters=NULL, 
                                      clusterExpandFactor=0.01,
                                      clusterLty=2, clusterLwd=2, 
                                      clusterBorder=rgb(0.5,0.5,0.5,1), 
                                      clusterFill=rgb(0,0,0, 0.1), snpScaleSize=2,
                                      clusterPCH=23, clusterPointCex=10,
                                      nodesDefiningClades, cladeColours,
                                      addBranchColours=FALSE, branchColour="dimgrey"){

  # Get and set the margins
  currentMar <- par()$mar
  par(mar=c(0,0,0,0))

  # Set the number plots in window - 1 row of 2
  # Leave more space for map in second column
  if(is.null(layoutMatrix)){
    layout(matrix(c(1,2,2,2,2), nrow=1, ncol=5, byrow=TRUE))
  }else{
    layout(layoutMatrix)
  }
  
  ### Plot the phylogeny on the left
  plotPhylogeny(tree, tipInfo, tipCex=tipShapeCexOnPhylogeny, scaleCex=scaleCex,
                tipShapesAndColours=tipShapesAndColours,
                addTipIndices=addTipIndices, indexCex=tipLabelCexOnPhylogeny,
                labelOffset=tipLabelOffset, showTipLabels=showTipLabels, 
                scaleSize=snpScaleSize, nodesDefiningClades=nodesDefiningClades, 
                cladeColours=cladeColours, addBranchColours=addBranchColours,
                branchColour=branchColour)
  
  # Note the tip coordinates
  tipCoordsOnPhylogeny <- getTipCoordinatesOnPhylogeny(tree$tip.label)
  
  ### Plot the sampling locations map
  
  # Plot the map
  par(mar=c(4,0,0,0))
  plot(x=NULL, y=NULL, xlim=range(tipInfo$X, na.rm=TRUE),
       ylim=range(tipInfo$Y, na.rm=TRUE),
       xaxt="n", yaxt="n", xlab="", ylab="", bty="n", asp=1)

  # Add scale bar
  addScale(size=10000, labelCex=scaleCex)
  
  # Plot shape around spatial clusters
  if(is.null(spatialClusters) == FALSE){
    plotShapeAroundSpatialClusters(tipInfo, spatialClusters,
                                   expand=clusterExpandFactor,
                                   lty=clusterLty, lwd=clusterLwd, 
                                   border=clusterBorder,
                                   fill=clusterFill,
                                   pch=clusterPCH,
                                   cex=clusterPointCex)
  }
  
  # Plot the sampling locations
  labels <- c(1:nrow(tipInfo))[is.na(tipInfo$X) == FALSE]
  xCoords <- tipInfo$X[is.na(tipInfo$X) == FALSE]
  yCoords <- tipInfo$Y[is.na(tipInfo$X) == FALSE]
  colours <- getTipShapeOrColourBasedOnSpecies(tipInfo, tipShapesAndColours,
                                               which="colour", alpha=0.75, 
                                               scoreForHerds=FALSE)[
                                                 is.na(tipInfo$X) == FALSE]
  addTextLabels(xCoords, yCoords, labels=labels,
                col.label=colours, avoidPoints=FALSE, cex.label=tipCexOnMap,
                col.background=tipIndexBackground)
  
  # Get coordinates of plotted sampling locations in plotting window
  tipCoordsOnMap <- getTipCoordinatesOnMap(tipInfo)
  
  ### Connect tips on phylogeny to sampling locations on map
  plotConnectingLines(tree$tip.label, tipCoordsOnPhylogeny,
                      tipCoordsOnMap, tipInfo, lwd=connectingLinesWidth,
                      alpha=connectingLinesAlpha)
  
  # Add species legend
  legend("topright", legend=names(tipShapesAndColours), bty="n", cex=scaleCex, pt.cex=tipCexOnMap,
         pch=as.numeric(unlist(tipShapesAndColours)[seq(from=2, by=2, length.out=length(tipShapesAndColours))]),
         col=unlist(tipShapesAndColours)[seq(from=1, by=2, length.out=length(tipShapesAndColours))])

  
  # Reset the margins
  par(mar=currentMar)
}

setAlpha <- function(colour, alpha){
  
  # Convert the input colour into rgb values
  rgbValues <- col2rgb(colour)
  
  # Place rgb values within rgb function and insert alpha value
  # Note that col2rgb returns rgbvlues from 0 to 255
  rgbColour <- rgb(rgbValues["red", 1], rgbValues["green", 1], rgbValues["blue", 1],
                   alpha=alpha*255, maxColorValue=255)
  return(rgbColour)
}

plotConnectingLines <- function(tipLabels, tipCoordsOnPhylogeny,
                                tipCoordsOnMap, tipInfo, alpha=1, ...){
  
  # Prepare for adding lines across the plot panels - using grid package
  pushViewport(viewport())
  popViewport()
  
  # Examine each of the tips
  for(tipLabel in tipLabels){
    
    # Skip the tips without locations
    if(is.na(tipCoordsOnMap[[tipLabel]][1])){
      next
    }
    
    # Define the line colour based on species
    lineColour <- "red"
    if(tipInfo[which(tipInfo$Tip == tipLabel), "Species"] == "COW"){
      lineColour <- "blue"
    }
    
    # Set the line alpha
    lineColour <- setAlpha(lineColour, alpha)
    
    # Prepare to add a single line
    pushViewport(viewport())
    
    # Plot line from phylogeny to sampling location for current tip
    grid.lines(x = c(tipCoordsOnPhylogeny[[tipLabel]][1],
                     tipCoordsOnMap[[tipLabel]][1]), 
               y = c(tipCoordsOnPhylogeny[[tipLabel]][2],
                     tipCoordsOnMap[[tipLabel]][2]), 
               gp = gpar(col=lineColour, ...))
    
    # Add the changes to the plot (the line)
    popViewport()
  }
  
}

getSatelliteImage <- function(landparcelCoords, badgerShapeFile, deerShapeFile, 
                              cattleShapeFile, includeLandParcels=FALSE,
                              xExpand=NULL, yExpand=NULL){
  
  # Initialise a variable to store the map bounds
  bounds <- NULL
  
  # Check bounds should include land parcels
  if(includeLandParcels){
    
    # Calculate bounds based on badger and deer sampling locations and herd land parcels
    bounds <- calculateBoundsForAllSpatialData(badgerShapeFile@bbox, deerShapeFile@bbox,
                                               cattleShapeFile@bbox)
  }else{
    
    # Calculate the herd centroids based upon the Irish grid coordinates
    herdCentroidsAsXY <- calculateLandParcelCentroids(landParcelCoords, useGrid=TRUE)
    
    # Calculate bounds just using badger and deer locations and overal herd centroids
    bounds <- calculateBoundsForAllSpatialData(badgerShapeFile@bbox, deerShapeFile@bbox,
                                               getLimitsOfHerdCentroids(herdCentroidsAsXY))
  }
  
  # Expand the X bounds if requested
  if(is.null(xExpand) == FALSE){
    
    # Calculate the X length
    xLength <- bounds[1, 2] - bounds[1, 1]
    
    # Adjust the bounds
    bounds[1, 1] <- bounds[1, 1] - (0.5*xExpand*xLength)
    bounds[1, 2] <- bounds[1, 2] + (0.5*xExpand*xLength)
  }
  
  # Expand the Y bounds if requested
  if(is.null(yExpand) == FALSE){
    
    # Calculate the X length
    yLength <- bounds[2, 2] - bounds[2, 1]
    
    # Adjust the bounds
    bounds[2, 1] <- bounds[2, 1] - (0.5*yExpand*yLength)
    bounds[2, 2] <- bounds[2, 2] + (0.5*yExpand*yLength)
  }
  
  # Convert the bounds to longitudes and latitudes
  boundsLongsLats <- convertXYToLatLongs(x=bounds[1, ],y=bounds[2, ])
  
  # Get a satellite image of the area
  map <- openmap(upperLeft=c(boundsLongsLats[2, 1], boundsLongsLats[1, 2]), 
                 lowerRight=c(boundsLongsLats[1, 1], boundsLongsLats[2, 2]),
                 type="bing")
  
  return(map)
}

getLimitsOfHerdCentroids <- function(herdCentroids){
  
  # Initialise a table to store the coordinates bounds
  bbox <- matrix(c(Inf, -Inf), nrow=2, ncol=2, byrow=TRUE)
  colnames(bbox) <- c("min", "max")
  rownames(bbox) <- c("x", "y")
  
  # Examine each of the herds
  for(herdID in names(herdCentroids)){
    
    # Check if we need to update the X min and max
    if(herdCentroids[[herdID]]$X < bbox["x", "min"]){
      bbox["x", "min"] <- herdCentroids[[herdID]]$X
    }
    if(herdCentroids[[herdID]]$X > bbox["x", "max"]){
      bbox["x", "max"] <- herdCentroids[[herdID]]$X
    }
    
    # Check if we need to update the Y min and max
    if(herdCentroids[[herdID]]$Y < bbox["y", "min"]){
      bbox["y", "min"] <- herdCentroids[[herdID]]$Y
    }
    if(herdCentroids[[herdID]]$Y > bbox["y", "max"]){
      bbox["y", "max"] <- herdCentroids[[herdID]]$Y
    }
  }
  
  return(bbox)
}

getTipCoordinatesOnMap <- function(tipInfo){
  
  # Convert the available mercator coordinates into grid coordinates
  tipInfo$GridX <- grconvertX(tipInfo$X, from="user", to="ndc")
  tipInfo$GridY <- grconvertY(tipInfo$Y, from="user", to="ndc")
  
  # Initialise a list to store the grid coordinates
  coords <- list()
  
  # For the cattle - convert the herd centroids to grid coordinates
  for(row in seq_len(nrow(tipInfo))){
    
    # Store the coordinates for the current row
    coords[[tipInfo[row, "Tip"]]] <- c(tipInfo[row, "GridX"], tipInfo[row, "GridY"])
  }
  
  return(coords)
}

getTipCoordinatesOnPhylogeny <- function(tipLabels){
  
  # Get all the information about the last phylogenetic tree plotted
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  
  # Create an empty list to store the coordinates
  tips <- list()
  
  # Examine each of the tip labels - order must match tree$tip.labels of plotted tree
  for(i in 1:length(tipLabels)){
    
    # Get and store the coordinates for the current tip label
    # Note that you are converting them to actual coordinates within the plotting window
    tips[[as.character(tipLabels[i])]] <- c(grconvertX(lastPP$xx[i], "user", "ndc"), 
                                            grconvertY(lastPP$yy[i], "user", "ndc"))
  }
  
  return(tips)
}

plotPhylogeny <- function(tree, tipInfo, tipCex=1.25, 
                          scaleCex=1, tipShapesAndColours,
                          addTipIndices=FALSE, indexCex=1, labelOffset=1,
                          showTipLabels=FALSE, scaleSize=2, nodesDefiningClades=NULL, 
                          cladeColours=NULL, addBranchColours=TRUE, branchColour="dimgrey",
                          type="phylogram"){
  
  # Get and set the plotting margins
  currentMar <- par()$mar
  par(mar=c(3,0,0,0))
  
  # Define the branch colours - branches in clades are coloured by clade
  branchColours <- branchColour
  if(addBranchColours){
    branchColours <- defineBranchColoursOfClades(tree, nodesDefiningClades,
                                                 cladeColours, "lightgrey")
  }
  
  # Plot the phylogeny
  if(addTipIndices){
    
    # Re-define tip labels as indices
    tree$tip.label <- 1:length(tree$tip.label)
    
    # Plot tree with tip labels
    plot.phylo(tree, show.tip.label=TRUE, edge.color=branchColours, edge.width=4,
               cex=indexCex, label.offset=labelOffset, xpd=TRUE, type=type,
               tip.color=getTipShapeOrColourBasedOnSpecies(tipInfo,
                                                           tipShapesAndColours,
                                                           which="colour",
                                                           scoreForHerds=FALSE))
  }else if(showTipLabels){
    plot.phylo(tree, show.tip.label=TRUE, edge.color=branchColours, edge.width=4,
               cex=indexCex, label.offset=labelOffset, xpd=TRUE, type=type)
  }else{
    plot.phylo(tree, show.tip.label=FALSE, edge.color=branchColours, edge.width=4, type=type)
  }
  
  # Add tips coloured by species
  tiplabels(pch=getTipShapeOrColourBasedOnSpecies(tipInfo,
                                                  tipShapesAndColours,
                                                  which="shape"),
            col=getTipShapeOrColourBasedOnSpecies(tipInfo,
                                                  tipShapesAndColours,
                                                  which="colour",scoreForHerds=FALSE),
            cex=tipCex)
  
  # Add scale bar
  addSNPScale(position="bottom", size=scaleSize, cex=scaleCex, lineWidth=4)
  
  # Reset the plotting margins
  par(mar=currentMar)
}

defineBranchColoursOfClades <- function(tree, nodesDefiningClades,
                                        CladeColours, defaultColour){
  branchColours <- rep(defaultColour, dim(tree$edge)[1])
  for(i in 1:length(nodesDefiningClades)){
    clade <- tips(tree, node=nodesDefiningClades[i])
    branchesInClades <- which.edge(tree, clade)
    branchColours[branchesInClades] <- cladeColours[i]
  }
  return(branchColours)
}


#### FUNCTIONS - spatial data ####

getSampledHerds <- function(tipInfo){
  
  # Get a list of the sampled herds
  sampledHerds <- tipInfo[is.na(tipInfo$Species) == FALSE & 
                            tipInfo$Species == "COW" & 
                            is.na(tipInfo$HerdOrSettID) == FALSE, "HerdOrSettID"]
  
  return(sampledHerds)
}

addScale <- function(size, xPad=0.1, yPad=0.1, lwd=1, labelCex=2){
  
  # Get the axis limits
  axisLimits <- par()$usr
  
  # Calculate the lengths of the X and Y axes
  xLength <- axisLimits[2] - axisLimits[1]
  yLength <- axisLimits[4] - axisLimits[3]
  
  # Note the bottom left corner of scale box
  xLeft <- axisLimits[1] + (xPad*xLength)
  yBottom <- axisLimits[3] + (yPad*yLength)
  
  # Plot a scale box
  rect(xleft=xLeft, ybottom=yBottom, xright=xLeft+size, ytop=yBottom+size, lwd=lwd)
  
  # Add a label
  label <- bquote(paste(.((size/1000)*(size/1000)), "km")^2)
  text(x=xLeft + (0.5*size), y=yBottom, labels=label, pos=1, cex=labelCex)
}

addHerdCentroidsAndSettLocationsToTipInfo <- function(tipInfo, settLocations,
                                                      landParcelCentroids){
  
  # Initialise columns to store the coordinates
  tipInfo$X <- NA
  tipInfo$Y <- NA
  tipInfo$HerdCentroidScore <- NA
  tipInfo$MeanDistanceToLandParcels <- NA
  
  # Examine each row of the tip information
  for(row in 1:nrow(tipInfo)){
    
    # Skip animals with no herd ID or sett ID
    if(is.na(tipInfo[row, "HerdOrSettID"])){
      next
    }
    
    # Check if cow
    if(tipInfo[row, "Species"] == "COW"){
      
      # Check if herd coordinates are available
      if(tipInfo[row, "HerdOrSettID"] %in% names(landParcelCentroids) == FALSE){
        warning("No land parcels available for herd: ", 
                tipInfo[row, "HerdOrSettID"])
        next
      }
      
      # Store herd centroid
      tipInfo[row, "X"] <- landParcelCentroids[[tipInfo[row, "HerdOrSettID"]]]$X
      tipInfo[row, "Y"] <- landParcelCentroids[[tipInfo[row, "HerdOrSettID"]]]$Y
      
      # Store the herd score
      tipInfo[row, "HerdCentroidScore"] <- landParcelCentroids[[tipInfo[row, "HerdOrSettID"]]]$Score
      tipInfo[row, "MeanDistanceToLandParcels"] <- landParcelCentroids[[tipInfo[row, "HerdOrSettID"]]]$SummaryOfDistancesToCentroid$Mean
          
    # Otherwise if badger store sett coordinates
    }else{
      
      # Check if sett coordinates are available
      settRow <- which(settLocations$SettID == tipInfo[row, "HerdOrSettID"])
      if(length(settRow) == 0){
        warning("No coordinates available for sett: ", 
                tipInfo[row, "HerdOrSettID"])
        next
      }
      
      # Store the sett coordinates
      tipInfo[row, "X"] <- settLocations[settRow, "X"]
      tipInfo[row, "Y"] <- settLocations[settRow, "Y"]
    }
  }
  
  return(tipInfo)
}

calculateLandParcelCentroids <- function(landParcelCoords){
  
  # Initialise a list to store the centroids associated with each herd
  herdCentroids <- list()
  
  # Examine each of the herd
  for(herd in names(landParcelCoords)){
    
    # Initialise two vectors to store the X and Y coordinates of each field's centroid
    centroidsX <- c()
    centroidsY <- c()
    
    # Examine the field polygons for the current herd
    for(fieldIndex in seq_len(length(landParcelCoords[[herd]]))){
      
      # Get the polygon coordinates for the current field
      coords <- landParcelCoords[[herd]][[fieldIndex]]
      
      # Calculate the centroid for the current field
      centroidsX[fieldIndex] <- mean(coords[, "X"])
      centroidsY[fieldIndex] <- mean(coords[, "Y"])

    }

    # Store the calculated centroids
    herdCentroids[[herd]] <- list(
      "X"=mean(centroidsX),
      "Y"=mean(centroidsY),
      "CentroidXs"=centroidsX,
      "CentroidYs"=centroidsY)
  }
  
  return(herdCentroids)
}

getPolygonCoords <- function(spatialDataFrame){
  
  # Got some good information from: https://stackoverflow.com/questions/29803253/r-extracting-coordinates-from-spatialpolygonsdataframe
  # Also from my blog: https://josephcrispell.github.io/BlogPosts/GetPolygonsFromShapeFile_11-10-17/GetPolygonsFromShapeFile_11-10-17.html
  # WHICH NEEDS UPDATED!!!!
  
  
  # Note the number of polygon sets - one for each ID
  nSets <- length(spatialDataFrame@polygons)
  
  # Initialise a list to store the IDs and coordinates of each polygon
  output <- list()
  
  # Loop through all the sets of polygons
  for(setIndex in seq_len(nSets)){
    
    # Get the ID of the current polygon set
    id <- spatialDataFrame@polygons[[setIndex]]@ID
    
    # Convert the ID to a herd ID
    id <- as.character(spatialDataFrame@data[as.numeric(id) + 1, "HERD_NO"])
    
    # Note the number of polygons in the current set
    nPolygons <- length(spatialDataFrame@polygons[[setIndex]]@Polygons)
    
    # Initialise a list to store the coordinates of each polygon in current set
    polygons <- list()
    
    # Examine each polygon in current set
    for(polygonIndex in seq_len(nPolygons)){
      
      # Get the coordinates for the current polygon
      coords <- spatialDataFrame@polygons[[setIndex]]@Polygons[[polygonIndex]]@coords
      colnames(coords) <- c("X", "Y")

      # Store the latitude and longitudes and X and Y coords
      polygons[[polygonIndex]] <- coords
    }
    
    # Store the polygon coordinates
    output[[id]] <- polygons
  }
  
  return(output)
}

calculateDistancesToOverallCentroidForEachHerd <- function(herdCentroids, tipInfo){
  
  # Get a list of the sampled herds
  sampledHerds <- getSampledHerds(tipInfo)
  
  # Create vector to store each herd's mean distance of land parcels to overall centroid
  meanDistancesToCentroids <- c()
  
  # Examine each herd
  for(herdID in names(herdCentroids)){
    
    # Skip unsampled herds
    if(herdID %in% sampledHerds == FALSE){
      next
    }
    
    # Initialise vector to store distance of each land parcel centroid to overall centroid
    distances <- c()
    
    # Examine each land parcel centroid
    for(index in seq_along(herdCentroids[[herdID]]$CentroidXs)){
      
      # Calculate the distance of the current land parcel centroid to overall centroid
      distances[index] <- euclideanDistance(herdCentroids[[herdID]]$X, 
                                            herdCentroids[[herdID]]$Y,
                                            herdCentroids[[herdID]]$CentroidXs[index],
                                            herdCentroids[[herdID]]$CentroidYs[index])
    }
    
    # Store the distances and summarise the distribution
    herdCentroids[[herdID]]$DistancesToCentroid <- distances
    herdCentroids[[herdID]]$SummaryOfDistancesToCentroid <- list("Mean"=mean(distances), 
                                                                 "Variance"=var(distances),
                                                                 "Min"=min(distances),
                                                                 "Max"=max(distances))
    
    # Store the mean distance for the current herd
    meanDistancesToCentroids[length(meanDistancesToCentroids) + 1] <- 
      herdCentroids[[herdID]]$SummaryOfDistancesToCentroid$Mean
  }
  
  # Calculate the maximum mean distance of land parcels to a centroid
  maxMean <- max(meanDistancesToCentroids)
  
  # Score each herd based upon how close its land parcels are to its centroid
  for(herdID in names(herdCentroids)){
    herdCentroids[[herdID]]$Score <- 1 - (herdCentroids[[herdID]]$SummaryOfDistancesToCentroid$Mean
                                          / maxMean)
  }
  
  return(herdCentroids)
}

#### FUNCTIONS - plotting temporal sampling ####

plotTemporalSamplingRange <- function(tipInfo){
  
  # Convert the date format
  tipInfo$DateAtSlaughterOrCapture <- as.Date(tipInfo$DateAtSlaughterOrCapture)
  
  # Calculate the range in the sampling dates
  dateRange <- range(tipInfo$DateAtSlaughterOrCapture, na.rm=TRUE)
  
  # Get and set the plotting margins
  currentMar <- par()$mar
  par(mar=c(4.1, 5.1, 4.1, 4.4))
  
  # Define the locations of the sampling ranges on the Y axis
  yLocations <- c(0.15, 0.35)
  speciesShapes <- c(19, 19)
  
  # Create an empty plot
  plot(x=NULL, y=NULL, ylim=c(0, 0.5), xlim=dateRange,
       bty="n", ylab="", yaxt="n", xlab="", xaxt="n",
       main="Number of samples available in time",
       cex.main=2)
  
  # Add an X axis
  at <- as.Date(c("2018-05-15", "2018-08-15", "2018-11-15", "2019-02-15", "2019-05-15"),  format="%Y-%m-%d")
  axis(side=1, at=at, labels=FALSE, xpd=TRUE)
  axis(side=1, at=at, labels=format(at, "%d-%m-%Y"), tick=FALSE, cex.axis=1.5,
       line=1, xpd=TRUE)
  
  # Add y axis
  species <- c("BADGER", "COW")
  axis(side=2, at=yLocations, labels=species, las=1, tick=FALSE, cex.axis=1.5,
       line=-2)
  
  # Examine the sampling of each species
  for(speciesIndex in seq_along(species)){
    
    # Get the sampling dates for the current species
    dates <- tipInfo[tipInfo$Species == species[speciesIndex], "DateAtSlaughterOrCapture"]
    
    # Remove NAs
    dates <- dates[is.na(dates) == FALSE]
    
    # Order dates
    dates <- sort(dates)
    
    # Plot the dates
    points(x=dates, y=rep(yLocations[speciesIndex], length(dates)),
           type="o", pch=speciesShapes[speciesIndex], col=rgb(0,0,0, 0.25), cex=5)
    
    # Note the number of samples available
    text(x=dateRange[2]-(0.001*(dateRange[2] - dateRange[1])),
         y=yLocations[speciesIndex], adj=0, xpd=TRUE,
         labels=paste0("     ", length(dates), "/",
                       length(which(tipInfo$Species == species[speciesIndex]))),
         cex=1.5)
  }
  
  # Reset the plotting margins
  par(mar=currentMar)
}

#### FUNCTIONS - linking to metadata ####

addSlaughterOrCaptureDates <- function(tipInfo, cattleTestInfoBatch1, cattleTestInfoBatch2,
                                       badgerCaptureData, settCaptureEventData){
  
  # Add a date into the tip info table
  tipInfo$DateAtSlaughterOrCapture <- NA
  
  # Examine each row of tip information
  for(row in seq_len(nrow(tipInfo))){
    
    # Skip if no animal ID available
    if(is.na(tipInfo[row, "AnimalID"])){
      next
    }
    
    # Check if cow
    if(tipInfo[row, "Species"] == "COW"){
      
      # Try to find test information row for current cow
      testRowBatch1 <- which(cattleTestInfoBatch1[, "Animal No."] == tipInfo[row, "AnimalID"])
      testRowBatch2 <- which(cattleTestInfoBatch2[, "Animal ID"] == tipInfo[row, "AnimalID"])
      
      # Check if found in batch 1
      if(length(testRowBatch1) > 0){
        
        # Get the unique slaughter dates for the current animal
        slaughterDates <- unique(as.Date(cattleTestInfoBatch1[testRowBatch1, "PM (Date of SL)"],
                                         format="%d/%m/%Y"))
        
        # Check if more than unique date available
        if(length(slaughterDates) > 1){
          warning("Multiple slaughter dates available for cow: ", tipInfo[row, "AnimalID"])
        }
        
        # Store a slaughter date
        tipInfo[row, "DateAtSlaughterOrCapture"] <- as.character(slaughterDates[1])
        
      # Check if found in batch 2
      }else if(length(testRowBatch2) > 0){
        
        # Get the test date
        testDate <- as.Date(cattleTestInfoBatch2[testRowBatch2, "Test Date"], format="%d/%m/%Y")
        
        # Store a test date
        tipInfo[row, "DateAtSlaughterOrCapture"] <- as.character(testDate)
        
      # Otherwise throw an error
      }else{
        warning("Not able to find testing information for cattle eartag: ", tipInfo[row, "AnimalID"])
      }
    
    # Get the capture date for the badger
    }else{
      
      # Try and find badger in the capture event table
      badgersRow <- which(badgerCaptureData$BADGER_ID == tipInfo[row, "AnimalID"])
      
      # Skip badger if didn't find a record
      if(length(badgersRow) == 0){
        warning(paste0("Unable to find capture data for badger: ", tipInfo[row, "AnimalID"]))
        next
      }
      
      # Store the sett ID
      tipInfo[row, "HerdOrSettID"] <- badgerCaptureData[badgersRow, "SETT_ID"]
      
      # Find the sett capture event information
      eventRow <- which(settCaptureEventData$CAPTURE_BLOCK_EVENT == badgerCaptureData[badgersRow, "CAPTURE_BLOCK_EVENT"])
      
      # Skip if weren't able to find capture event row
      if(length(eventRow) == 0){
        warning(paste0("Unable to find information for sett capture event:: ", badgerCaptureData[badgersRow, "CAPTURE_BLOCK_EVENT"]))
        next
      }
      
      # Estimate the date of capture for the current badger
      tipInfo[row, "DateAtSlaughterOrCapture"] <- estimateCaptureDate(settCaptureEventData[eventRow, "DATE_COMMENCED"],
                                                                      settCaptureEventData[eventRow, "DATE_COMPLETED"])
    }
  }
  
  return(tipInfo)
}

estimateCaptureDate <- function(dateTrappingCommenced, dateTrappingCompleted){
  
  # Convert the date strings to date objects
  dateTrappingCommenced <- as.Date(strsplit(dateTrappingCommenced, split=" ")[[1]][1], format="%d-%b-%y")
  dateTrappingCompleted <- as.Date(strsplit(dateTrappingCompleted, split=" ")[[1]][1], format="%d-%b-%y")
  
  # Get the middle date
  middleDate <- dateTrappingCommenced + ((dateTrappingCompleted - dateTrappingCommenced)/2)
  
  return(as.character(middleDate))
}

getTipInfo <- function(tipLabels, animalTags, coverage, herdIDs){
  
  # Initialise a dataframe to store the tip information
  tipInfo <- data.frame("Tip"=tipLabels,"Aliquot"=NA, "AnimalID"=NA, "Species"=NA,
                        "Coverage"=NA, "HerdOrSettID"=NA, 
                        stringsAsFactors=FALSE)
  
  # Examine each of the tips
  for(index in seq_along(tipLabels)){
    
    # Split the tip label if dash
    aliquot <- tipLabels[index]
    if(grepl(tipLabels[index], pattern="-")){
      aliquot <- strsplit(tipLabels[index], split="-")[[1]][2]
    }
    
    # Build an aliquot code for the current isolate
    aliquot <- paste0("TB19-", paste(rep(0, 6-nchar(aliquot)), collapse=""), aliquot)
    tipInfo[index, "Aliquot"] <- aliquot
    
    # Find the row in the coverage information for the current tip
    coverageRow <- which(grepl(coverage$IsolateID, pattern=paste0("^", tipLabels[index], "_"))
                         | grepl(coverage$IsolateID, pattern=paste0("^", tipLabels[index], "p_")))
    
    # Check that row was found
    if(length(coverageRow) == 0){
      warning("Unable to find coverage information for: ", tipLabels[index],
              "\tAliquot: ", aliquot)
    }else if(length(coverageRow) > 1){
      warning("Multiple entries in coverage information for: ", tipLabels[index],
              "\tAliquot: ", aliquot)
    }else{
      # Store the coverage information
      tipInfo[index, "Coverage"] <- coverage[coverageRow, "PercentageCoverage"]
    }
    
    # Find the row in the sample information table for the current tip
    tagRow <- which(animalTags$Aliquot == aliquot)
    
    # Check that row was found
    if(length(tagRow) == 0){
      warning("Unable to find sampling information for: ", tipLabels[index],
              "\tAliquot: ", aliquot)
    }else if(length(tagRow) > 1){
      warning("Multiple entries in sampling information for: ", tipLabels[index],
              "\tAliquot: ", aliquot)
    }else{
      # Store the animal ID
      tipInfo[index, "AnimalID"] <- animalTags[tagRow, "Animal.ID"]
      
      # Store the species of the current tip
      tipInfo[index, "Species"] <- ifelse(grepl(tipInfo[index, "AnimalID"], pattern="^RR"),
                                          "BADGER", "COW")
    }
    
    # Skip badgers
    if(is.na(tipInfo[index, "Species"]) == FALSE && 
       tipInfo[index, "Species"] == "BADGER"){
      next
    }
    
    # Find the row in the herd IDs for the current tip (note we'll ignore badgers)
    herdRow <- which(herdIDs$Aliquot == aliquot)
    
    # Check that row was found
    if(length(herdRow) == 0){
      warning("Unable to find herd for: ", tipLabels[index],
              "\tAliquot: ", aliquot)
    }else if(length(herdRow) > 1){
      warning("Multiple entries in herd IDs table for: ", tipLabels[index],
              "\tAliquot: ", aliquot)
    }else{
      # Store the animal ID
      tipInfo[index, "HerdOrSettID"] <- herdIDs[herdRow, "Herd_ID"]
    }
  }
  
  return(tipInfo)
}

#### FUNCTIONS - plotting phylogeny ####

editTipLabels <- function(tipLabels){
  
  # Initialise a vector to store the new labels
  output <- c()
  
  # Examine each tip label
  for(label in tipLabels){
    
    # Remove the ">" prefix
    if(grepl(label, pattern="^>")){
      label <- substr(label, 2, nchar(label))
    }
    
    # Split the label and retain first part
    label <- strsplit(label, split="_")[[1]][1]
    
    # Remove a trailing "p"
    label <- gsub("p", "", label)
    
    # Store the new label
    output[length(output) + 1] <- label
  }
  
  return(output)
}

#### FUNCTIONS - Phylogeny ####

getTipShapeOrColourBasedOnSpecies <- function(tipInfo, tipShapesAndColours, which,
                                              alpha=1, scoreForHerds=TRUE, shapeIndex=2){
  
  # Initialise a vector to store the shapes or colours
  output <- c()
  
  # Examine each tip
  for(row in seq_len(nrow(tipInfo))){
    
    # Check if Species available
    if(is.na(tipInfo[row, "Species"])){
      
      # Check if wanting shape or colour
      if(which == "shape"){
        output[row] <- as.numeric(tipShapesAndColours[["NA"]][shapeIndex])
      }else if(which == "colour"){
        output[row] <- setAlpha(tipShapesAndColours[["NA"]][1], alpha)
      }else{
        cat(paste("Error! Option", which, "not recognised!"))
      }
      
      # If species available assign appropriate colour or shape
    }else{
      
      # Check if wanting shape or colour
      if(which == "shape"){
        output[row] <- as.numeric(tipShapesAndColours[[tipInfo[row, "Species"]]][shapeIndex])
      }else if(which == "colour"){
        
        # Get the species colour
        colour <- tipShapesAndColours[[tipInfo[row, "Species"]]][1]
        
        # Check if herd centroid score available
        if(is.na(tipInfo[row, "HerdCentroidScore"]) == FALSE && scoreForHerds == TRUE){
          colour <- setAlpha(tipShapesAndColours[[tipInfo[row, "Species"]]][1], 
                             tipInfo[row, "HerdCentroidScore"])
        }else{
          colour <- setAlpha(tipShapesAndColours[[tipInfo[row, "Species"]]][1], alpha)
        }
        
        # Set the colour
        output[row] <- colour
        
      }else{
        cat(paste("Error! Option", which, "not recognised!"))
      }
    }
  }
  
  return(output)
}

getNSitesInFASTA <- function(fastaFile){
  
  # Open a connection to a file to read (open="r")
  connection <- file(fastaFile, open="r")
  
  # Get first line of file
  firstLine <- readLines(connection, n=1)
  
  # Close file connection
  close(connection)
  
  # Get the number of sites used in the FASTA file from first line
  nSites <- as.numeric(strsplit(firstLine, " ")[[1]][2])
  
  return(nSites)
}

runRAXML <- function(fastaFile, date, path, nBootstraps=100, nThreads=6, alreadyRun=FALSE, outgroup=NULL, model="GTRCAT"){
  
  # Create a directory for the output file
  directory <- paste(path, "RAxML_", date, sep="")
  suppressWarnings(dir.create(directory))
  
  # Set the Working directory - this will be where the output files are dumped
  setwd(directory)
  
  # Build analysis name
  analysisName <- paste("RaxML-R_", date, sep="")
  
  # Check if already Run and just want to retrieve tree
  if(alreadyRun == FALSE){
    
    # Build the command
    seeds <- sample(1:100000000, size=2, replace=FALSE) # For parsimony tree and boostrapping
    
    if(is.null(outgroup)){
      command <- paste("raxmlHPC", 
                       " -f a", # Algorithm: Rapid boostrap inference
                       " -N ", nBootstraps,
                       " -T ", nThreads,
                       " -m ", model, " -V", # -V means no rate heterogenity
                       " -p ", seeds[1], " -x ", seeds[2], # Parsimony and boostrapping seeds
                       " -n ", analysisName,
                       " -s ", fastaFile, sep="")
    }else{
      command <- paste("raxmlHPC", 
                       " -f a", # Algorithm: Rapid boostrap inference
                       " -N ", nBootstraps,
                       " -T ", nThreads,
                       " -m ", model, " -V", # -V means no rate heterogenity
                       " -p ", seeds[1], " -x ", seeds[2], # Parsimony and boostrapping seeds
                       " -n ", analysisName,
                       " -s ", fastaFile, 
                       " -o ", outgroup, sep="")
    }
    
    system(command, intern=TRUE)
  }
  
  # Get the tree and read it in
  treeBS <- getTreeFileWithSupportValues(analysisName)
  
  return(treeBS)
}

getTreeFileWithSupportValues <- function(analysisName){
  
  # Get files in current working directory
  files <- list.files()
  
  # Select the tree file with BS support values
  treeBSFile <- files[grepl(files, pattern=paste("RAxML_bipartitions[.]", analysisName, sep="")) == TRUE]
  
  # Open the file
  treeBS <- read.tree(treeBSFile)
  
  return(treeBS)
}
