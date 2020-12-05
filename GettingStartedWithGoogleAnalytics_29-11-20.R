# Resources:
# - http://www.dartistics.com/api-google-analytics.html
# - https://www.rubenvezzoli.online/googleanalyticsr-beginners-guide/
# - view results by page: https://www.dataquest.io/blog/tutorial-blog-post-analysis-with-r-googleanalyticsr/


# Load libraries
library(googleAnalyticsR)

# Authorise google analytics
ga_auth() # Note that I had to confirm in email before permissions were granted

# Pull a full list of the views that you have access to
myAccounts <- ga_account_list()
id <- myAccounts[1, "viewId"] # Note this is the same ID as shown in url when logging into Google Analytics after the "p"

# Pull in the google analytics data
webData <- google_analytics(id, 
                            date_range = c("2020-01-01", "2020-11-01"),
                            metrics = c("users", "timeOnPage"), # meta data.frame (already in environment) provides all metrics available via API
                            dimensions = c("date", "channelGrouping", "country", "latitude", "longitude", "pagePath"),
                            anti_sample = TRUE) # Split up the call to avoid sampling

                            