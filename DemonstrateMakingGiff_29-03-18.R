#########################################
# Script to demonstrate generating GIFF #
#########################################

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/ExampleGiff/"

##################################
# Create some example PNG images #
##################################

# Open a PNG file creator - creates and names each plot created as PNG
pathToPNGFiles <- paste(path, "ExampleImage_%02d.png", sep="")
png(file=pathToPNGFiles)

# Set the plot margins
par(mar=c(0, 0, 0, 0))

# Create plots to store in PNG files
for(i in 1:10){
  
  # Simple plot - plotting current i as coordinate
  plot(x=i, y=i, 
       xlim=c(1,10), ylim=c(1,10),
       xlab="", ylab="", xaxt="n", yaxt="n",
       las=1, bty="n", pch=19, cex=2)
}

# Close the png
dev.off()

# Reset plotting margins
par(mar=c(5.1, 4.1, 4.1, 2.1))

#######################################
# Combine all the PNG files into GIFF #
#######################################

# Note the path to the PNG files and convert to windows format
pathToPNGFiles <- paste(path, "ExampleImage_*.png", sep="")
formatedPathToPNGFiles <- shQuote(pathToPNGFiles, type="cmd")

# Create a name for the GIFF and convert to windows format
pathToGiffFile <- paste(path, "Example.gif", sep="")
formatedPathToGiffFile <- shQuote(pathToGiffFile, type="cmd")

# Use Magick to knit PNG files into GIFF
system(paste("magick -delay 50 ", formatedPathToPNGFiles, " ", formatedPathToGiffFile, sep=""))

# Delete the PNG files
unlink(pathToPNGFiles)
