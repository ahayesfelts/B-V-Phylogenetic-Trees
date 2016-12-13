# DATA LOADING of Pre-Aligned Sequences 
# load necessary packages
library(ape)
library(phangorn)
library(seqinr)

# import dataset with headers and column begining with row 1
# datast is in phylip and sequential format
# make sure working directory is set to the directory with dataset in it
ADAM7_phylip_sequential <- read.csv("~/Desktop/R/ADAM7_phylip_sequential.dna", sep="")
View(ADAM7_phylip_sequential)
adam7DNA  <-  ADAM7_phylip_sequential
#view dataset
View(adam7DNA)
# CONVERT DATA TO phyDat or pml class
# read.phyDat converts different file types to the pml or phyDat object 
# specify path to file in question with double quotes, type of data, and format of data
adam7DNA <- read.phyDat("~/Desktop/R/ADAM7_phylip_sequential.dna", type = "DNA", format = "phylip")

# MODEL SELECTION with modeltest
# compares different nucleotide/AA models with AIC, AICc, or BIC
# will only work if you've converted object to phyDat class
# default model for modelTest is for DNA
adam7mt <- modelTest(adam7DNA)
env <- attr(adam7mt, "env")
# list models in order 
ls(env = env)
# View model scores for adam7mt
View(adam7mt)
# choose model based on lowest loglikelihood score
# print out table of statistical scores, rate matrix, base frequencies
(fit <- eval(get("GTR+G+I", env), env))
# returns the loglikelihood, estimated parameters, AIC, AICc, and BIC tested models
# "env" dataframe has all the trees, data, and calls to get estimated models for starting point
# returns necessary values for distance modeling 
(GTRGI <- get("GTR+G+I", env))
# end of model selection, go to Distance Based Methods

# DISTANCE BASED METHODS

# create distance matrix from alignment data, using parameters from output of distance modeling
adam7dm  <- dist.ml(adam7DNA, bf = c(0.312409668248493, 0.216599733955101, 0.211665043116414, 0.259325554679992), Q = c(1.53447482678452, 3.99970934720379, 0.80489435056284, 1.29064215061731, 5.01487852044787, 1), inv = 0.121568347559429, k = 4, shape = 1.76956618505029)
# create UPGMA rooted tree with dm object from previous step
adamtreeUPGMA <- upgma(adam7dm)
# create NJ unrooted tree with dm object, will have to root later
adamtreeNJ <- NJ(adam7dm)
plot(adamtreeUPGMA, type = "phylogram", main="UPGMA")             



