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
