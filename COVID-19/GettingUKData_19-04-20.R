ukCasesFile <- "https://raw.githubusercontent.com/tomwhite/covid-19-uk-data/master/data/covid-19-cases-uk.csv"
ukCases <- read.csv(ukCasesFile, stringsAsFactors=FALSE, check.names=FALSE)

# Convert case counts to numbers
ukCases$TotalCases <- as.numeric(ukCases$TotalCases)

# Collapse to counts for each country on each date
ukCases$CountryDate <- paste0(ukCases$Country, "_", ukCases$Date)
countsPerCountry <- aggregate(ukCases$TotalCases, by=list(ukCases$CountryDate), FUN=sum, na.rm=TRUE)

# Extract the country names and dates from aggregated table
countsPerCountry$Date <- sapply(countsPerCountry$Group.1, FUN=function(countryWithDate){strsplit(countryWithDate, "_")[[1]][2]})
countsPerCountry$Country <- sapply(countsPerCountry$Group.1, FUN=function(countryWithDate){strsplit(countryWithDate, "_")[[1]][1]})

# Build a table reporting case count in each country on each date
casesPerCountry <- data.frame("Date"=sort(unique(countsPerCountry$Date)), "Scotland"=0, "England"=0, "Wales"=0, "Northern Ireland"=0, 
                              stringsAsFactors=FALSE, check.names=FALSE)
rownames(casesPerCountry) <- casesPerCountry$Date
for(row in 1:nrow(countsPerCountry)){
  casesPerCountry[countsPerCountry[row, "Date"], countsPerCountry[row, "Country"]] <- countsPerCountry[row, "x"]
}

# Add a total column
casesPerCountry[, "Grand Total"] <- rowSums(casesPerCountry[, 2:4])



# Read in the UK shape file
shapeFile <- file.path(path, "Countries_December_2019_Boundaries_UK_BFE", 
                       "Countries_December_2019_Boundaries_UK_BFE.shp")
uk <- rgdal::readOGR(shapeFile)
