# DATA LOADING
library(ape)
library(phangorn)
library(seqinr)
library(readr)
primates <- read_csv("http://evolution.gs.washington.edu/book/primates.dna")
View(primates)

# CONVERT DATA TO phyDat class
primates <- read.phyDat("primates.dna", format = "phylip", type = "DNA")

# MODEL SELECTION with modeltest
primatesMT <- modelTest(primates)
# see results of comparison between models 
View(primatesMT)
# results of modelTest saved as call which can be evaluated to be a pml object
env <- attr(primatesMT, "env")
# list all possible models 
ls(envir = env)
# get the optimal fit for the substitution model chosen (GTR+G in this case)
# how the model is selected depends on what we're looking for in data
(fit <- eval(get("GTR+G", env), env))

# DISTANCE BASED to get base trees
# dist.ml to create a distance matrix
primatesDM <- dist.ml(primates)
# reconstruct a rooted tree
primatesUPGMA <- upgma(primatesDM)
# reconstruct an unrooted tree
primatesNJ <- NJ(primatesDM)
# parameters to plot both UPGMA and NJ on same plot
layout(matrix(c(1,2), 2, 1), height=c(1,2))
par(mar = c(0,0,2,0)+ 0.1)
# plot distance trees
plot(primatesUPGMA, main="UPGMA")
plot(primatesNJ, "unrooted", main="NJ")

# PARSIMONY
# returns parsimony score, number of changes which are
# at least needed to describe data for the input tree
parsimony(primatesUPGMA, primates)
parsimony(primatesNJ, primates)
# optimize (find lowest) parsimony score by rearranging the trees
UPGMApars <- optim.parsimony(primatesUPGMA, primates)
# using parsimony ratchet scores
primatesRatchet <-  pratchet(primates, trace = 0)
#compare the two scores 
parsimony(c(UPGMApars, primatesRatchet), primates)

# MAXIMUM LIKELIHOOD
# find the likelihood for a tree given the data
# returns object of class pml
# additional parameters available if not using NJ 
fit <- pml(primatesNJ, data = primates)
# view parameters of the model 
fit
# plot the tree for the object that estimated the likelihood
# of the tree suplied 
plot(fit)
# then optimize the branch lengths for the Jukes-Cantor model
fitJC <- optim.pml(fit, TRUE)
# gets log-likelihood value from maximum likelihood branch lengths
# degrees of freedom (df) gives number of estimated model parameters
logLik(fitJC)
# default pml values estimated JC model
fitGTR <- update(fit, k=4)
# update will change the model to the GTR+G and optimize all of its parameters
fitGTR <- optim.pml(fitGTR, model = "GTR", optGamma = TRUE, rearrangement = "ratchet", control = pml.control(trace = 0))
# the control parameters allow the thresholds for the fitting process to be changed
fitGTR

# MODEL COMPARISON
#compare nexted models for the JC and GTR+G with likelihood ratio stat
anova(fitJC, fitGTR)
# with the Shimodaira-Hasegawa test
SH.test(fitGTR, fitJC)

# with the AIC - amount of info lost when we use a particular model to
# approximate the real process of nucleotide substition, smaller is better
AIC(fitJC)
AIC(fitGTR)
# AIC bias correction is for regression and autoregression time series models
# good when sample size is small or when amount of fitted parameters is moderate
# to large fraction of the sample size
AICc(fitGTR)
#Bayesian information criterion, given equal priors for competing models
# model with smallest BIC is equivalent with max posterior probabilities 
BIC(fitGTR)

# BOOTSTRAP
# test how well the edges of tree are supported
bs = bootstrap.pml(fitJC, bs=100, optNni=TRUE, control = pml.control(trace = 0))

# Plot Bootstrap via phangorn
# placing two trees in one plot
par(mfrow=c(2,1))
par(mar=c(1,1,3,1))
# plot the unrooted tree (midpoint rooted) with bs values
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")
title("A")
#consensusnet from bootstrap to identify any potential conflicts
cnet <- consensusNet(bs, p=0.2)
plot(cnet, "2D", show.edge.label=TRUE)
title("B")


