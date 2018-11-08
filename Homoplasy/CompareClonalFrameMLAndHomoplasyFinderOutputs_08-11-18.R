#### Load libraries ####

library(ape)

#### Load the HomoplasyFinder output ####

# Set the path
path <- "/home/josephcrispell/Desktop/ClonalFrameML_SaureusData/"

# Read in the HomoplasyFinder report
homoplasyFinderResults <- read.table(paste0(path, "consistencyIndexReport_BEFORE_07-11-18.txt"), sep="\t", header=TRUE)

#### Get the ClonalFrameML output ####
# Majority of this code taken directly from: https://github.com/xavierdidelot/ClonalFrameML/blob/master/src/cfml_results.R

# Read options from command line
prefix <- paste0(path, "example.output")
coresites_list <- paste0(path, "Saureus_core-sites.txt")
treefile <- paste(prefix,".labelled_tree.newick",sep="")
xreffile <- paste(prefix,".position_cross_reference.txt",sep="")
ML_seqfile <- paste(prefix,".ML_sequence.fasta",sep="")

# Get a vector of the core sites
coresites <- scan(coresites_list)

# Load a list cross-referencing patterns in the original data to the output FASTA file
xref <- scan(xreffile,sep=",")

# Count the patterns in the output sequences on the phylogeny
n.mut <- countPatternsOnPhylogeny(ML_seqfile, treefile)

#### Compare ClonalFrameML and HomoplasyFinder outputs ####

# Some notes about the ClonalFrameML R code about
# n.mut               For each mutation, how many times did it occur on the tree? In output sequences??
# xref                Index is reference sequence position. Value is output sequences position. Zero means no mutation
#                     2,902,619 long. 106,480 aren't zeros. These refer to only 20,790 unique positions in output sequences
# output sequences    One position for each UNIQUE pattern in the input sequences
#                     Therefore multiple reference genome positions can be linked to a position in the output sequences

# Compare the ClonalFrameML and HomoplasyFinder outputs
output <- compareClonalFrameMLAndHomoplasyFinderOutputs(n.mut, homoplasyFinderResults)

# Summarise the comparison
table(output$HomoplasyFinderPositionsFound)
table(output$ClonalFrameMLPositionsFound, useNA="ifany")
table(n.mut > 1)

# All 13743 homoplasies identified by ClonalFrameML were identified by HomoplasyFinder
# 19810/34364 (~58%) of homoplasies identified by HomoplasyFinder were identified by ClonalFrameML
#       Possibly some of these homoplasies are actually recombination events that ClonalFrameML has discounted

#### FUNCTIONS ####

compareClonalFrameMLAndHomoplasyFinderOutputs <- function(n.mut, homoplasyFinderResults){
  
  # Create a vector to record which of the homoplasies identified by ClonalFrameML were found by HomoplasyFinder
  found <- rep(NA, length(n.mut))
  
  # Add a column to the HomoplasyFinder results to identify which homoplasies it found were identified by ClonalFrameML
  homoplasyFinderResults$FoundByClonalFrame <- FALSE
  
  # Examine the number of times the patterns in the output sequence positions were observed on the phylogeny
  for(outputSeqPosition in seq_along(n.mut)){
    
    # Print progress information
    cat(paste("\r", outputSeqPosition, "of", length(n.mut)))
    
    # Skip when number of times pattern observed is one or less - non-homoplasious
    if(n.mut[outputSeqPosition] < 2){
      next
    }
    
    # Assume current homoplasious position not found
    found[outputSeqPosition] <- FALSE
    
    # Identify positions on input sequence that had the current pattern
    refPositions <- which(xref == outputSeqPosition)
    
    # Check which were found by HomoplasyFinder
    presentInHomoplasyFinderOutput <- (refPositions %in% homoplasyFinderResults$Position)
    
    # Examine each of the positions in the input sequence that the pattern was found at
    for(positionInInputSeq in seq_along(refPositions)){
      
      # If found by HomoplasyFinder - record that in found vector and in HomoplasyFinder output table
      if(presentInHomoplasyFinderOutput[positionInInputSeq]){
        homoplasyFinderResults[which(homoplasyFinderResults$Position == refPositions[positionInInputSeq]), "FoundByClonalFrame"] <- TRUE
        found[outputSeqPosition] <- TRUE
      }
    }
  }
  
  # Create a list to store the outputs
  output <- list(
    "HomoplasyFinderPositionsFound"=homoplasyFinderResults$FoundByClonalFrame,
    "ClonalFrameMLPositionsFound"=found)
  
  return(output)
}

countPatternsOnPhylogeny <- function(ML_seqfile, treefile){
  
  # All code for this function taken directly from:
  #  https://github.com/xavierdidelot/ClonalFrameML/blob/master/src/cfml_results.R
  
  # Load the phyML tree estimated from all core variant and invariant sites
  tree <- read.tree(treefile)
  
  # Load the imputed and reconstructed ancestral sequences
  ML_seq<-scan(ML_seqfile,what=character(0))
  tp <- substr(ML_seq[seq(1,length(ML_seq),by=2)],2,1000)
  ML_seq <- ML_seq[seq(2,length(ML_seq),by=2)]; names(ML_seq) = tp
  # M is a matrix containing the FASTA file base calls
  M <- matrix("",length(ML_seq),nchar(ML_seq[1]))
  for(i in 1:length(ML_seq)) {
    v <- unlist(strsplit(ML_seq[i],""))
    M[i,] <- v
    gc()
  }
  rownames(M) <- names(ML_seq)
  
  # Combine the tip and node labels
  treelabels <- c(tree$tip.label,tree$node.label)
  
  # For each row of M, identify the node index
  M_node_index <- match(rownames(M),treelabels)
  
  # For each row of M, identify the node index of its ancestor
  # To do this, identify the node index in tree$edge[,2] and read tree$edge[,1]
  M_anc_node_index <- tree$edge[match(M_node_index,tree$edge[,2]),1]
  
  # Find, by name, the ancestor
  M_anc_node <- treelabels[M_anc_node_index]
  
  # Find its position in M
  M_anc_node_M_index <- match(M_anc_node,rownames(M))
  
  # Precompute the positions of mutations on branches of the tree
  # For each pattern, record the mutated nodes
  # wh.mut is a matrix, in the same order as M, recording whether the base represents a mutation
  wh.mut <- apply(M,2,function(m) 1*(m!=m[M_anc_node_M_index])); wh.mut[nrow(wh.mut),] = 0
  
  # A homoplasy is a mutation that occurs on multiple branches. Count the number of homoplasic mutations per branch
  # Exclude reference sequences from the count
  gd <- !is.na(as.numeric(rownames(wh.mut))) | substr(rownames(wh.mut),1,4)=="NODE"
  n.mut <- apply(wh.mut[gd,],2,sum)
  
  return(n.mut)
}
