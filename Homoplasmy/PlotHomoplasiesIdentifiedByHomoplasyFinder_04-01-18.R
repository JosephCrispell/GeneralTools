###############
# Preparation #
###############

# Packages
library(phangorn)
library(gplots)

# Set the path
#path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/NewZealand/NewAnalyses_12-05-16/MLTree/"
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_13-07-17/vcfFiles/";

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

################################################
# Read in the java HomoplasyFinder tool output #
################################################

file <- paste(path, "homoplasyReport_", date, ".txt", sep="")
homoplasyInfo <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                            check.names=FALSE, comment.char = "@")

#########################
# Read in the phylogeny #
#########################

# "#" in isolate IDs messes up read.tree - replace in tree and FASTA file for New Zealand
#file <- paste(path, "mlTree_withRef_ReplacedHashes_14-06-16.tree", sep="")
file <- paste(path, "mlTree_29-09-2017.tree", sep="")
tree <- read.tree(file)

##############################################################
# Plot the phylogeny and annotate the homoplasies identified #
##############################################################

pdf(paste(path, "test.pdf", sep=""))
plotPhylogenyAndAnnotateHomoplasyEvents(tree, homoplasyInfo, nNodesBack=0)
dev.off()

hist(homoplasyInfo$NNodesBackToRootAlleleFound)

#############
# FUNCTIONS #
#############

plotPhylogenyAndAnnotateHomoplasyEvents <- function(tree, homoplasyInfo, nNodesBack=0){
  
  par(mar=c(0, 0, 0, 0))
  
  # Create a list to note the homoplasy events already examined
  done <- list()
  
  # Examine each homoplasy event identified
  for(row in 1:nrow(homoplasyInfo)){
    
    # Skip homoplasies that are found only 1 node back
    if(homoplasyInfo[row, "NNodesBackToRootAlleleFound"] <= nNodesBack){
      next
    }
    
    # Define the tip colors based upon the homoplasies identified
    tipColours <- rep(rgb(0,0,0, 1), length(tree$tip.label))
    
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
      
      if(tree$tip.label[i] %in% isolates == FALSE){
        next
      }
      tipColours[i] <- "red"
    }

    # Plot the tree
    plot.phylo(tree, show.tip.label=TRUE, type="fan",
               edge.color="grey", edge.width=3,
               show.node.label=TRUE, label.offset=0,
               tip.color=tipColours, cex=0.25)
    
    # Get the axis limits
    axisLimits <- par("usr")
    
    # Add a scale
    xPos <- axisLimits[2] - (0.5 * (axisLimits[2] - axisLimits[1]))
    text(x=xPos, y=axisLimits[4], 
         labels=paste("Pos: ", homoplasyInfo[row, "Position"],
                      "    Allele: ", homoplasyInfo[row, "Allele"], sep=""), pos=1)
    
    # Note that finished with homoplasy event
    done[[paste(idFoundIn, alsoFoundIn, sep=":")]] <- 1
  }
  
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}
