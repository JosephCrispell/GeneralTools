# Resources:
# - http://www.dartistics.com/api-google-analytics.html
# - https://www.rubenvezzoli.online/googleanalyticsr-beginners-guide/
# - view results by page: https://www.dataquest.io/blog/tutorial-blog-post-analysis-with-r-googleanalyticsr/
# - using dimension! https://rstudio-pubs-static.s3.amazonaws.com/276165_33dcdf0b0ed04ac99b5ffb04d6863868.html

#### Preparation ####

# Load libraries
library(googleAnalyticsR)

#### Authorisation ####

# Authorise google analytics
ga_auth() # Note that I had to confirm in email before permissions were granted

#### Download data ####

# Pull a full list of the views that you have access to
myAccounts <- ga_account_list()
id <- myAccounts[1, "viewId"] # Note this is the same ID as shown in url when logging into Google Analytics after the "p"

# Define the date range
start <- format(seq(Sys.Date(), length=2, by="-6 months")[2], "%Y-%m-%d")
today <- format(Sys.Date(), "%Y-%m-%d")

# Pull in the google analytics data
webData <- google_analytics(id, 
                            date_range = c(start, today), # When the pull data from
                            metrics = c("users", "timeOnPage"), # Which metrics to report for chosen dimentions
                            dimensions = c("date", "channelGrouping", "country", "latitude", "longitude", "pagePath"), # Single row given to each unique value in these
                            anti_sample = TRUE) # Split up the call to avoid sampling

#### Process data ####

#### Visualise ####                          