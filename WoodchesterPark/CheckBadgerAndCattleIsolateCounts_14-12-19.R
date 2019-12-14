
# Set the path
path <- "/home/josephcrispell/storage/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/"

file <- paste(path, "vcfFiles/",  "PreviousMergings/1_22-03-18/isolateCoverageSummary_DP-20_22-03-2018.txt", sep="")
table <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

file  <- cattleInfoFile <- paste(path, "IsolateData/",
                                 "CattleIsolateInfo_AddedNew_TB1484-TB1490_22-03-18.csv",
                                 sep="")
cattleInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

count <- 0

for(id in cattleInfo$StrainId){
  
  if(is.na(id)){
    next
  }
  
  if(length(which(grepl(table$IsolateID, pattern=id))) == 1){
    
    count <- count + 1
  }
}

badgerInfoFile <- paste(path, "IsolateData/",
                        "BadgerInfo_08-04-15_LatLongs_XY_Centroids.csv", sep="")
badgerInfo <- read.table(badgerInfoFile, header=TRUE, sep=",", stringsAsFactors=FALSE)

tattoos <- c()
for(row in 1:nrow(table)){
  
  id <- strsplit(table[row, "IsolateID"], split="_")[[1]][1]
  
  if(grepl(id, pattern="WB") == FALSE){
    next
  }
  
  tattoos[length(tattoos) + 1] <- badgerInfo[which(badgerInfo$WB_id == id), "tattoo"]
}