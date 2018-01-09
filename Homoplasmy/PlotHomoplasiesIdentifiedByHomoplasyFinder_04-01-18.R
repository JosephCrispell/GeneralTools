###############
# Preparation #
###############

# Packages
library(phangorn)
library(gplots)

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Homoplasmy/"

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

################################################
# Read in the java HomoplasyFinder tool output #
################################################

file <- paste(path, "homoplasyReport_", date, ".txt", sep="")
homoplasyInfo <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                            check.names=FALSE)

#########################
# Read in the phylogeny #
#########################

file <- paste(path, "example_", date, ".tree", sep="")
tree <- read.tree(file)

##############################################################
# Plot the phylogeny and annotate the homoplasies identified #
##############################################################

plotPhylogenyAndAnnotateHomoplasyEvents(tree, homoplasyInfo)


#############
# FUNCTIONS #
#############

plotPhylogenyAndAnnotateHomoplasyEvents <- function(tree, homoplasyInfo){
  
  par(mar=c(0, 0, 0, 0))
  
  # Create a list to note the homoplasy events already examined
  done <- list()
  
  # Examine each homoplasy event identified
  for(row in 1:nrow(homoplasyInfo)){
    
    # Define the tip colors based upon the homoplasies identified
    tipColours <- rep("black", length(tree$tip.label))
    
    # Create a key for the homoplasy
    idFoundIn <- homoplasyInfo[row, "Node/IsolateID"]
    alsoFoundIn <- homoplasyInfo[row, "IsolatesAlleleFoundIn"]
    
    # Check if already observed and skip
    if(is.null(done[[paste(alsoFoundIn, idFoundIn, sep=":")]]) == FALSE){
      next
    }
    
    # Get the isolates associated with the current homoplasy identified
    isolates <- c(strsplit(idFoundIn, split="-")[[1]],
                  strsplit(alsoFoundIn, split="-")[[1]])
    
    # Colour the isolates according to the current homoplasy event
    for(i in 1:length(tree$tip.label)){
      
      if(tree$tip.label[i] %in% isolates == FALSE || tipColours[i] != "black"){
        next
      }
      tipColours[i] <- "red"
    }

    # Plot the tree
    plot.phylo(tree, show.tip.label=TRUE, type="phylogram",
               edge.color="grey", edge.width=3,
               show.node.label=TRUE, label.offset=0.15,
               tip.color=tipColours, cex=1.5)
    
    # Get the axis limits
    axisLimits <- par("usr")
    
    # Add scale
    xPos <- axisLimits[2] - (0.1 * (axisLimits[2] - axisLimits[1]))
    yPos <- axisLimits[3] + (0.075 * (axisLimits[4] - axisLimits[3]))
    lines(x=c(xPos,xPos+1), y=c(yPos, yPos), lwd=4)
    text(x=xPos+0.5, y=axisLimits[3], labels="1 SNP", pos=3)
    
    # Note that finished with homoplasy event
    done[[paste(idFoundIn, alsoFoundIn, sep=":")]] <- 1
  }
  
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}
