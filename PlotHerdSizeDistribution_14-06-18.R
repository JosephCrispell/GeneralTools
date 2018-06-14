path <- "/home/josephcrispell/Desktop/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/CattleMovementData-Post2001/"

table <- read.table(paste(path, "herdSizeDistribution_14-06-2018.txt", sep=""), header=TRUE, fill=TRUE, stringsAsFactors=FALSE)

dates <- as.Date(c("15/06/2003", "15/06/2004", "15/06/2005", "15/06/2006", "15/06/2007", "15/06/2008", "15/06/2009",
           "15/06/2010", "15/06/2011", "15/06/2012", "15/06/2013"), format="%d/%m/%Y")

remove <- c()
for(row in 1:nrow(table)){
  
  # Skip none BEEF, DAIRY or MIXED herds
  if(table[row, "HerdType"] %in% c("BEEF", "DAIRY", "MIXED") == FALSE){
    remove[length(remove) + 1] <- row
  }
  
  # Skip herds that were always empty
  uniqueHerdSizes <- unique(as.numeric(table[row, 3:(length(dates)+2)]))
  if(length(uniqueHerdSizes) == 1 && uniqueHerdSizes[1] == 0){
    remove[length(remove) + 1] <- row
  }
}
table <- table[-remove, ]

plot(dates, table[1, 3:(length(dates)+2)], 
     ylim=range(table[, 3:(length(dates)+2)]),
     xlim=range(dates),
     las=1, bty="n", type="l", pch=19,
     col="white",
     ylab="Herd Size", xlab="Year", main="Herd size distribution through time")
for(row in 1:nrow(table)){
  
  points(dates, table[row, 3:(length(dates)+2)], type="l", pch=19,
         col=rgb(0,0,0, 0.25))
}
