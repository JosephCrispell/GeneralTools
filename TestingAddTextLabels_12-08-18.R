library(addTextLabels)

x = c(0.8846, 1.1554, 0.9317, 0.9703, 0.9053, 0.9454, 1.0146, 0.9012, 
      0.9055, 1.3307)
y = c(0.9828, 1.0329, 0.931, 1.3794, 0.9273, 0.9605, 1.0259, 0.9542, 
      0.9717, 0.9357)
ShortSci = c("MotAlb", "PruMod", "EriRub", "LusMeg", "PhoOch", "PhoPhoppp", 
             "SaxRub", "TurMer", "TurPil", "A really long label")

# basic plot
plot(x, y, asp=1)
abline(h = 1, col = "green")
abline(v = 1, col = "green")

addTextLabels(xCoords=x, yCoords=y, labels=ShortSci, col.line="grey", cex=1,
              col.background=rgb(0,0,0, 0.5))

x <- 1.1
y <- 1.35
label <- "!!!!!!!!!!ggggg"
text(x=x, y=y, labels=label)
xHalf <- strwidth(label) * 0.51
yHalf <- strheight(label) * 0.75
polygon(x=c(x - xHalf,
            x - xHalf,
            x + xHalf,
            x + xHalf),
        y=c(y - yHalf,
            y + yHalf,
            y + yHalf,
            y - yHalf), 
        col=rgb(0,0,0, 0.2),
        border=NA, xpd=TRUE)