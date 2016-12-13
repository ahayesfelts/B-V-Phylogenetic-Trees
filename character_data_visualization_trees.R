#############################################################################
#Tutorial on visualization of character data on phylogenetic trees in R.
#The data for this tutorial are derived from Price et al. 2012. Elevated
#rates of morphological and functional diversification in reef-dwelling
#Haemulid fishes. Evolution 67:417-428.
#These data are available via the Dryad data repository: http://datadryad.org/resource/doi:10.5061/dryad.s049s
#Tutorial by R. Glor, March 2014
#############################################################################

##############################
#Preliminaries
##############################
rm(list = ls()) #clear your work environment
setwd("~/Dropbox/Projects/Bodega_Bay/Dryad_Data/") #set working directory. You'll obviously need to change this line to reflect the location of the relevant files on your computer

#load relevant libraries
library(ape)
library(geiger)
library(diversitree)

##############################
#Load tree
##############################
read.nexus("Trees4dryad.nex") -> haemulidTrees #Load set of 500 ultrametric trees obtained from analyses of a multilocus data in BEAST
plot(haemulidTrees[[1]], cex=0.5) #Plot one of our trees to ensure that they loaded properly and are ultrametric. The double square brackets designate the first tree in a set of 500 trees 
lapply(haemulidTrees, ladderize) -> haemulidTreesLadderized #Use the lapply and ladderize functions to ladderize all 500 of our trees. Ladderization is a largely aesthetic change to how the tree is plotted
plot(haemulidTreesLadderized[[1]], cex=0.5) #Plot the first ladderized tree
add.scale.bar() #Add a simple scale bar indicating the scale for the branches in your tree
lapply(haemulidTreesLadderized, rescale, "depth", 1) -> haemulidTreesLadderized1 #Use the lapply and rescale functions to rescale all 500 of our trees to have a root to tip distance of 1. This rescaling will make subsequent plotting functions somewhat easier. Even more importantly, it will often improve the performance of likelihood functions.
plot(haemulidTreesLadderized1[[1]], cex=0.5) #Plot the first ladderized and rescaled tree
add.scale.bar() #Confirm that the tree has been rescaled

#In order to use trees in R we need to learn how trees are stored in R
help(read.tree) #Obtain information on how trees are read and stored in R. Notice the various elements of a tree object included under the Value heading
haemulidTreesLadderized1[[1]]$tip.label #Show list of tip labels
haemulidTreesLadderized1[[1]]$edge #Show list of branches designated by their start and end points
haemulidTreesLadderized1[[1]]$edge.length #Show list of branch lengths
plot(haemulidTreesLadderized1[[1]], label.offset=0.05, cex=0.5) #plot our tree
nodelabels(cex=0.5) #Plot node labels. Confirm that rows in the $edge component correspond with branches in this labeled tree
tiplabels(cex=0.5) #Plot tip labels. Confirm that numbers in the plot correspond with row numbers in the $tip.label component of the tree

##############################
#Load data
##############################
read.csv("haemulidDryadData.csv") -> haemulidDataInput 
# Load character data for haemulids. 
# This file includes only a subset of the full dataset available via Dryad. 
# Saved in comma delimited text format, this input file should include six columns with taxonomic, habitat, morphological and functional data: genus, species, habitat (reef[0] or non-reef[1]), standard length, raker length, and suction index.
haemulidDataInput 
#Check out your data
data.frame(haemulidDataInput[,3:6]) -> haemulidData 
#Convert the habitat, phenotypic, and functional data into a dataframe in R. 
# Dataframes are helpful in R because they store each column of a matrix as a separate object, which will make them easier to call later.
paste(haemulidDataInput[,1], haemulidDataInput[,2], sep="_", collapse=NULL) -> rownames(haemulidData) 
#Use the paste function to combine the genus and species names in columns 1 and 2 of our data table. 
#The purpose of doing this is to obtain row names for our table that correspond with the names of taxa in our phylogenetic tree, where OTUs are designated as genus_species.
name.check(haemulidTreesLadderized1[[1]], haemulidData) 
#Use the name.check function of geiger to test if the same taxa are found in the phylogenetic tree and the dataset. 
#You should see that many taxa differ between the two datasets, many apparently due to typographical errors.

read.csv("haemulidDryadData_corrected.csv") -> haemulidDataInput #Load in a new dataset where names in the data matrix have been changed to match those in the tree
data.frame(haemulidDataInput[,3:6]) -> haemulidData
paste(haemulidDataInput[,1], haemulidDataInput[,2], sep="_", collapse=NULL) -> rownames(haemulidData)
name.check(haemulidTreesLadderized1[[1]], haemulidData)
haemulidData[! rownames(haemulidData) %in% name.check(haemulidTreesLadderized1[[1]], haemulidData)$data_not_tree,] -> haemulidData #Delete taxa from the habitat dataset that are not in tree (i.e., all those taxa that are in the $data_not_tree component of the name.check output)
name.check(haemulidTreesLadderized1[[1]], haemulidData) #Confirm that tree and data are now in agreement

##############################
#Plot the data
##############################
plot(haemulidTreesLadderized1[[1]], cex=0.5, label.offset=0.05) #Plot a tree that leaves some room between the tree tips and taxon labels so that we can plot habitat use in this space
haemulidHabitatLabel  haemulidTreesLadderized1 #Reread the new tree file
plot(haemulidTreesLadderized1[[1]], cex=0.5, label.offset=0.05)
tiplabels(cex=0.5) #You should now see that the sequence of taxon labels in your tree object corresponds with the sequence of taxon labels in your plot.

plot(haemulidTreesLadderized1[[1]], cex=0.5, y.lim=c(0,55), x.lim=c(0,3))
points(rep(2, length(haemulidTreesLadderized1[[1]]$tip.label)), 1:length(haemulidTreesLadderized1[[1]]$tip.label), pch=21, bg=haemulidHabitatLabel[match(haemulidTreesLadderized1[[1]]$tip.label, names(haemulidHabitatLabel))], cex=1) #Check to see if the tips are now lined up accordingly
text(2, 52, "Habitat", cex=0.4)
#Add some quantitative data
abline(v=2.1, col="gray")
segments(2.1, 1:nrow(haemulidData), 2.1 + (haemulidData[,2]/1200), 1:nrow(haemulidData))
text(2.02, 0, "Std. Length", pos=4, cex=0.4)
abline(v=2.35, col="gray")
segments(2.35, 1:nrow(haemulidData), 2.35 + (haemulidData[,3]/50), 1:nrow(haemulidData))
text(2.3, 0, "Raker L.", pos=4, cex=0.4)
abline(v=2.6, col="gray")
points(rep(2.6, nrow(haemulidData)), 1:nrow(haemulidData), pch=21, bg="blue", col="white", lwd=0.25, cex=haemulidData[,4]*2.6)
text(2.6, 52, "Suction Index", cex=0.4)