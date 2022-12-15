### Install packages
#install.packages(c("poppr", "mmod", "magrittr", "treemap"), 
#                 repos = "http://cran.rstudio.com", dependencies = TRUE)
### Call libraries
library("poppr")
monpop <- read.genalex("monpop.csv")
### geneclone is a format of a document for population genetics
monpop
