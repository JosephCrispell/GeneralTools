library(basicPlotteR)

# Create some random points
n <- 45
testLabels <- c("short", "mediummm", "looooonnnnnnngggggg", "0090292002", "9", "A Different label")
coords <- data.frame(X=runif(n), Y=runif(n, min=0, max=100), Name=sample(testLabels, size=n, replace=TRUE),
                     stringsAsFactors = FALSE)

# Plot them without labels
plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col="red", xlab="X", ylab="Y")

pdf("/Users/josephcrispell/Desktop/addTextLabels.pdf", width=14, height=14)

# Set the margins
par(mar=c(0.5,0.5,0.5,0.5))

# Set the background colour
par(bg = rgb(0,0,1, 0.25))

# With potentially overlapping labels
plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col=rgb(1,0,0, 0.5), xlab="", ylab="", cex=2)
text(coords$X, coords$Y, labels=coords$Name, xpd=TRUE, cex=2)

# Plot them with non-overlapping labels
plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col=rgb(1,0,0, 0.5), xlab="", ylab="", cex=2)
addTextLabels(coords$X, coords$Y, coords$Name, cex.label=2, col.label="black", lty=2, col.line=rgb(0,0,0, 0.5))

# Plot them with non-overlapping labels
plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col=rgb(1,0,0, 0.5), xlab="", ylab="", cex=2)
addTextLabels(coords$X, coords$Y, coords$Name, cex.label=2, col.background=rgb(0,0,0, 0.5), col.label="white", lty=2, col.line=rgb(0,0,0, 0.5))

dev.off()


x = c(0.8846, 1.1554, 0.9317, 0.9703, 0.9053, 0.9454, 1.0146, 0.9012, 
      0.9055, 1.3307)
y = c(0.9828, 1.0329, 0.931, 1.3794, 0.9273, 0.9605, 1.0259, 0.9542, 
      0.9717, 0.9357)
ShortSci = c("MotAlb", "PruMod", "EriRub", "LusMeg", "PhoOch", "PhoPho", 
             "SaxRub", "TurMer", "TurPil", "TurPhi")

# basic plot
plot(x, y, asp=1)
abline(h = 1, col = "green")
abline(v = 1, col = "green")

addTextLabels(x, y, ShortSci, cex=0.9, col.background=rgb(0,0,0, 0.75), col.label="white")
addTextLabels(x, y, ShortSci, cex=0.9, col.background=rgb(0,0,0, 0.75), col.label="white", keepLabelsInside = FALSE)

xx=10**(-4:-1)
yy=0.5*(xx**0.5)

#this doesn't show the label
plot(xx,yy,log="x",type="b")

addTextLabels(xx,yy,c("A", "B", "C", "D"))

#this works
plot(log10(xx),log10(yy),type="b")
addTextLabels(log10(min(xx)),log10(min(yy)),"(min,min)")

### Plotting with  multiple colours
# Create random points
x = runif(n=10)
y = runif(n=10)

# Create random labels
testLabels <- c("short", "mediummm", "looooonnnnnnngggggg", "0090292002", "9", "A Different label")
labels <- sample(testLabels, size=10, replace=TRUE)

# Set label colours
colours <- c("red", "red", "red", "red", "red", "blue", "blue", "blue", "blue", "blue")

# Plot points
plot(x, y, asp=1)

# Add labels
addTextLabels(x, y, labels, cex=0.9, col.background=rgb(0,0,0, 0.75), col.label=colours, lty=c(1,2,3))


