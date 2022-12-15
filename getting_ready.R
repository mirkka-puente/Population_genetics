### Information about the tutorial (all the tutorial)
# https://grunwaldlab.github.io/Population_Genetics_in_R/Data_Preparation.html

### Install packages -------
#install.packages(c("poppr", "mmod", "magrittr", "treemap"), 
#                 repos = "http://cran.rstudio.com", dependencies = TRUE)

### Call libraries --------
library("poppr")
monpop <- read.genalex("monpop.csv")
### geneclone is a format of a document for population genetics
#monpop

#Data pinf
# Load the data
data("Pinf")
Pinf
gac <- genotype_curve(Pinf, sample = 1000, quiet = TRUE)
#MLG = multilocus genotypes

#Data microbov
data("microbov")
microbov
gac <- genotype_curve(microbov, sample = 1000, quiet = TRUE)

#Data nancycats
data("nancycats")
nancycats
gac <- genotype_curve(nancycats, sample = 1000, quiet = TRUE)

#Data Aeut
data("Aeut")
Aeut
gac <- genotype_curve(Aeut, sample = 1000, quiet = TRUE)


#Calculating allele frequencies gene diversity --------
#Look at the missing data and D is Simpson, H index
pinflt <- locus_table(Pinf)


info_table(Pinf, type = "missing", plot = TRUE)
tail(genind2df(Pinf, sep = "/"))
Pinf.ploidy <- info_table(Pinf, type = "ploidy",
                          plot = TRUE, low = "black", high = "orange")

#Calculating genotypic diversity per population ---------
poppr(Pinf)

P.tab <- mlg.table(Pinf)


#clone correction
#We need to delete the isolates that have repetitive alleles
splitStrata(monpop) <- ~Tree/Year/Symptom
monpop # After (Three distinct levels)

mcc_TY <- clonecorrect(monpop, strata = ~Tree/Year, keep = 1:2)
mcc_TY

# Set the original data to the same strata
setPop(monpop) <- ~Tree/Year

#Looking at the difference of cloned corrected data and uncorrected data.
cc <- locus_table(mcc_TY, info = FALSE)
mp <- locus_table(monpop, info = FALSE)
mp - cc

#Graph
locus_diff <- mp - cc

# Note that I need to select the column containing Simpson's Index. That's
# labeled as "1-D".
barplot(locus_diff[, "1-D"], ylab = "Change in Simpson's Index", xlab = "Locus",
        main = "Comparison of clone-corrected vs. uncorrected data")
#the clone corrected have higher gene diversity values (the differences were more negative)
#We are removing clones so the data is more accuarate

#Genotypic richness
setPop(monpop) <- ~Symptom
monpop
(monpop_diversity <- poppr(monpop))
#BB has a higher diversity, sexual and principal cause infection

library("vegan")
mon.tab <- mlg.table(monpop, plot = FALSE)
min_sample <- min(rowSums(mon.tab))
rarecurve(mon.tab, sample = min_sample, xlab = "Sample Size", ylab = "Expected MLGs")
title("Rarefaction of Fruit Rot and Blossom Blight")
#we can do this graph or the table with monpop_diversity

#LINKAGE DESEQUILIBRIUM
# P >= 0.01 --> accept Ho ==> sexual reproduction
# P < 0.01 --> Reject Ho ==> asexual 

library("magrittr")
MX <- popsub(Pinf, "North America")
ia(MX, sample = 999)

MX %>% clonecorrect(strata= ~Continent/Country) %>% ia(sample = 999)
MX
#p-vaue increased, now it is sexual reproduction


SA <- popsub(Pinf, "South America")
ia(SA, sample = 999)

#Let's see if it is not a result from clonality
SA %>% clonecorrect(strata= ~Continent/Country) %>% ia(sample=999)
# It's not a product of clonality because the analysis kept similar p-values
#since p-value < 0.01 then it is asexual reproduction for sure. 
