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


#Population structure (compare values between populations)
library("mmod")

#Now we will use Hendrick’s standardized GST to assess population 
#structure among these populations

Gst_Hedrick(nancycats)

#Genetic distance calculation
library("ape") # To visualize the tree using the "nj" function
library("magrittr")
data(microbov)
set.seed(10)
ten_samples <- sample(nInd(microbov), 10)
mic10       <- microbov[ten_samples]
(micdist    <- provesti.dist(mic10))

#create a neighbor-joining tree.

theTree <- micdist %>%
  nj() %>%    # calculate neighbor-joining tree
  ladderize() # organize branches by clade
plot(theTree)
add.scale.bar(length = 0.05) # add a scale bar showing 5% difference.

set.seed(999)
aboot(mic10, dist = provesti.dist, sample = 200, tree = "nj", 
      cutoff = 50, quiet = TRUE)

strata(microbov) <- data.frame(other(microbov))
microbov

nameStrata(microbov) <- ~Country/Breed/Species

# Analysis
set.seed(999)
microbov %>%
  genind2genpop(pop = ~Country/Breed) %>%
  aboot(cutoff = 50, quiet = TRUE, sample = 1000, distance = nei.dist)


#K-means hierarchical clustering
#Software figures it out if the groups are close related or not
MX <- popsub(Pinf, "North America")
MXclust <- find.clusters(MX)
MXclust


SA <- popsub(Pinf, "South America")
SAclust <- find.clusters(SA)
SAclust

#TREE ----
pinfreps <- c(2, 2, 6, 2, 2, 2, 2, 2, 3, 3, 2)
#11 loci - need to include for the software analysis
MXtree <- bruvo.boot(MX, replen = pinfreps, cutoff = 50, quiet = TRUE)
SAtree <- bruvo.boot(SA, replen = pinfreps, cutoff = 50, quiet = TRUE)


#Let’s see how the groups we found with the clustering algorithm match up:
library("ape")
cols <- rainbow(4)
plot.phylo(MXtree, cex = 0.8, font = 2, adj = 0, tip.color = cols[MXclust$grp],
           label.offset = 0.0125)
nodelabels(MXtree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,
           font = 3, xpd = TRUE)
axisPhylo(3)
#more distance = more different between groups

#SA
plot.phylo(SAtree, cex = 0.8, font = 2, adj = 0, tip.color = cols[SAclust$grp],
           label.offset = 0.0125)
nodelabels(SAtree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,
           font = 3, xpd = TRUE)
axisPhylo(3)
#Everything clusters together nicely, further supporting a non-panmictic population.

#Interactive viewer 
imsn()
library("poppr")
library("magrittr")
data(monpop)
splitStrata(monpop) <- ~Tree/Year/Symptom
summary(monpop)

t26 <- monpop %>% setPop(~Tree) %>% popsub("26") %>% setPop(~Year/Symptom)
t26
# Set up our repeat lengths and populations to analyze
reps <- c(CHMFc4 = 7, CHMFc5 = 2, CHMFc12 = 4,
          SEA = 4, SED = 4, SEE = 2, SEG = 6,
          SEI = 3, SEL = 4, SEN = 2,
          SEP = 4, SEQ = 2, SER = 4)

sub9 <- c("9_BB", "9_FR")

# Calculate the MSN
t26.9msn <- bruvo.msn(t26, replen = reps, sublist = sub9, showplot = FALSE)

# Visualize the network
set.seed(120)
plot_poppr_msn(t26, t26.9msn, inds = "none", palette = cm.colors, nodebase = 1.25)
args(plot_poppr_msn)

