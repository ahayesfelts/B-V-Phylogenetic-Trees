 # DATA LOADING	

# to read in datasets:
  library(phangorn)
  data("yeast")
  yeast
  baseFreq(yeast, all = TRUE)
    # all = TRUE will give counts for bases, ambiguous codes, missing data, and alignment gaps
          
                         
# TREE BUILDING METHODS
                         
## Distance Based Methods 
                         
### dist.ml function offers common sub models for AA
### use object from reading in alignments for distances
  dm <- dist.ml(yeast)
  
  # if you want to check if sequences are too different to define evolutionary distances
  yeastdm <- dist.dna(yeast, model = "raw")
  # if using substitution model = "JC69" or model = "F81":
  yeastdmJC <- dist.ml(yeast, model = "JC69")
  yeastdmF81 <- dist.ml(yeast, model = "F81")
  # distance values are very close showing that empirical base frequencies are very close to equal base frequencies
### sample code below uses UPGMA model for rooted tree
### + and NJ for unrooted tree (comparative)
    treeUPGMA <- upgma(dm)
    treeNJ <– NJ(dm)
                         
### PLOT TREES CODE
###	creating starting trees for MP and ML analyses
    layout(matrix(c(1,2), 2, 1), height=c(1,2))
    par(mar = c(0,0,2,0)+ 0.1)
    # parameters for plotting allow both to be in same plot
    plot(treeUPGMA, main="UPGMA")             
    plot(treeNJ, "unrooted", main="NJ")

    
    
## Parsimony Scores
## parsimony function returns parsimony score i.e. number of changes which are at 
## + least necessary to describe data for given tree
### useful for comparing parsimony scores between trees
    parsimony(treeUPGMA, yeast)
    parsimony(treeNJ, yeast) 


## optim.parsimony function does tree arrangements to find trees with 
## + lower parsimony score
    treeNNI  <-  optim.parsimony(treeUPGMA, yeast)
    ### prathchet function above is parsimony ratchet algorithm
    treeRatchet  <- pratchet(yeast, trace = 0)
    parsimony(c(treePars, treeRatchet), yeast)
    # parsimony compares the parsimony scores for both algorithms



## Maximum Likelihood (ML)

## pml function shows loglikelihood, unconstrained loglikelihood, rate matrix, base frequences
    fitUPGMA  <- pml(treeUPGMA, yeast)
    fitUPGMA


## check to see which functions we can use with pml
  methods(class="pml")
  [1] AICc BIC anova logLik plot print simSeq update vcov


## use function optim.pml
## object fit code above estimated likelihood for tree supplied 
## + and optimze branch lengths for Jukes-Cantor model with default pml parameters
  fitJC <- optim.pml(fitUPGMAh, TRUE)
  # TRUE refers to optimizing edge weights
  # other options to optimize other parts
  logLik(fitJC)
  # extracts loglikihood and returns degrees of freedom with estimated
  # + number of parameters for the model
## use function update.pml 
## + to change parameters for pml 
### code below changes model to F81 model and optimizes parameters
  fitF81 <- update(fit, k=4, inv=0.2)
  fitF81 <- optim.pml(fitF81, model="F81", optInv=TRUE, optGamma=TRUE,
                      + rearrangement = "NNI", control = pml.control(trace = 0))
  # takes a moment, be patient
  fitF81
  # pml.control function in optim.pml code above controlls fitting process

## NNI rearrangements can get stuck in local maxima with larger trees
## better trees but slower execution
## to improve topology search:
### 1) set rearrangement = "stochastic" for stochastic rearrangements
### + by random NNI permutation to tree that gets optimized to escape local optima
### 2) set rearrangement = "ratchet" that performs the likelihood ratchet 
      fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                            + rearrangement = "stochastic", control = pml.control(trace = 0))
      fitGTR
    

# MODEL SELECTION

## can compare nested models for JC and F81 model 
## + using likelihood ratio statistic 
### with Shimodaira-Hasegawa test in code below
  anova(fitJC, fitF81)
  
  SH.test(fitF81, fitJC)

### with AIC
  AIC(fitJC)
  AIC(fitF81)
  AICc(fitF81)
  BIC(fitF81)

## can use modelTest function to compare DNA or AA models 
## + with the AIC, AICc, or BIC 
## + to compute and optimize many models
### similar to the programs ModelTest and ProtTest
  mt <- modelTest(yeast)
## modelTest thresholds for optimization not as strict as optim.pml
## the call can be evaluated to get a pml object for further analysis
  env <- attr(mt, "env")
  
  ls(envir = env)
  # returns table with models in order of score
  
  fitF81 <- eval(get("F81", env), env)
  fitJC <- eval(get("JC", env), env)
  
  # returns loglikelihoods, etc for that model's values for data
   

## finally apply bootstrap to test support for tree edges
  bs <- bootstrap.pml(fitF81, bs = 100, optNni = TRUE, control = pml.control(trace = 0))
  # use fit object that best fits dataset from previous analyses 
### to plot tree with bootstrap support 
  par(mfrow = c(2,1))
  par(mar = c(1,1,3,1))
  ### par function controls graphical parameters
  plotBS(midpoint(fitF81$tree), bs, p = 50, type = "p", main = "a")


  
  


# TREE VISUALIZATION with ggtree
  
  # with phangorn data
  fitGTRtree <-  phyPML(fitGTR, type="ml")
  ggtree(fitGTRtree) + geom_text(aes(x=branch, vjust=-.5))

## plot basic phylogram
  ggtree(fitGTR) + ggtitle("(Phylogram) rectangular layout")
## change look of tree
  ggtree(fitGTR, color="firebrick", size=1, linetype="dotted")
## slanted phylogram
  ggtree(fitGTR, layout="slanted") + ggtitle("(Phylogram) slanted layout")
## circular phylogram
  ggtree(fitGTR, layout="circular") + ggtitle("(Phylogram) circular layout")
  # with display labels
  ggtree(fitGTR, layout="circular") + geom_tiplab(aes(angle=angle), color='blue')
## fan phylogram
  ggtree(fitGTR, layout="fan", open.angle=180) + ggtitle("(Phylogram) circular layout")
## set branch.length = "none"within ggtree function

## display nodes/tips
  ggtree(fitGTR)+geom_point(aes(shape=isTip, color=isTip), size=3)
  # other option
  p <- ggtree(fitGTR) + geom_nodepoint(color="#b5e521", alpha=1/4, size=10)
  p + geom_tippoint(color="#FDAC4F", shape=8, size=3)
  # display labels
  p + geom_tiplab(size=3, color="purple")

  
# visualize list of trees
  # bootstrapped trees in list
  btrees <- read.tree(system.file("extdata/RAxML", "RAxML_bootstrap.H3", package="ggtree"))
  ggtree(btrees) + facet_wrap(~.id, ncol=10)
  btrees
  # bootstrapped trees merged for density tree, with best tree on top
  p <- ggtree(btrees, layout="rectangular",   color="lightblue", alpha=.3)
  
  best_tree <- read.tree(system.file("extdata/RAxML", "RAxML_bipartitionsBranchLabels.H3", package="ggtree"))
  df <- fortify(best_tree, branch.length = 'none')
  p+geom_tree(data=df, color='firebrick')

  
### Visualize tree with multiple sequence alignment 
### + with msaplot function
  fasta <- system.file(yeast, package="ggtree")
  msaplot(ggtree(beast_tree), fasta) 
### show specific slice of alignment with window parameter
  msaplot(ggtree(beast_tree), fasta, window=c(150, 200)) + coord_polar(theta='y')
  
