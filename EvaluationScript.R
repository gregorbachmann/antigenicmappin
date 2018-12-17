################################################################################
# Script for evaluating algos
################################################################################

# Evaluate all different methods using all different evaluation methods

# ---------------------------------- Load Libraries ------------------------------------ #
library(ggplot2)
library(reshape)
library(clusterCrit)
library(phytools)
library(gridExtra)

#---------------------------------- Initialize needed stuff --------------------------------#
# Evaluate using these dimensions for true shape space
dimensions <- seq(2,20,by=1)
# Readin the given HIV tree
tree <- read.tree(file="/Users/gregorbachmann/Desktop/ETH/Semester9/CompBioProject/Data/HIV_Newick_tree")
# Fix cluster assignment for given tree by eye and assign 10 antibodies
clusterIndexTree <- c(1,1,1,1,1,1,1,1,2,2,2,3,3,3,3,3)
numClusters <- 3
numA <- 9
clusterIndexA <- c(1,1,1,2,2,2,3,3,3)


#--------------------------------------Random Clusters--------------------------------------#

# Initialize storage arrays
iso_box <- matrix(0, length(dimensions))        # Isomap
int_box <- matrix(0, length(dimensions))        # Interval
spline_box <- matrix(0, length(dimensions))     # Spline
ord_box <- matrix(0, length(dimensions))        # Ordinal
index <- 1

for (dim in dimensions){
  #iso_box[index] <- mean(test_algos(dim=dim, numtests=500, numpoints=60,method="Isomap"))
  #int_box[index] <- mean(test_algos(dim=dim, numtests=500, numpoints=60,method="interval"))
  #ord_box[index] <- mean(test_algos(dim=dim, numtests=500, numpoints=60,method="ordinal"))
  spline_box[index] <- mean(test_algos(dim=dim, numtests=500, numpoints=60,method="mspline"))
  index <- index + 1
}

# --------------------------- Transform data into dataframes ---------------------------- #
isobox <- as.data.frame(cbind(1:length(dimensions),iso_box))
colnames(intbox) <- c("Dimension","Score")
intbox <- as.data.frame(cbind(1:length(dimensions),int_box))
colnames(isobox) <- c("Dimension","Score")
ordbox <- as.data.frame(cbind(1:length(dimensions),ord_box))
colnames(ordbox) <- c("Dimension","Score")
splinebox <- as.data.frame(cbind(1:length(dimensions),spline_box))
colnames(splinebox) <- c("Dimension","Score")

# ------------------------------------ Plotting -------------------------------------------#
box = merge(isobox, intbox, by="Dimension")
colnames(box) <- c("Dimension","Isomap","Interval")
prescriptionMelted <- reshape2::melt(box, id.var='Dimension')
box = merge(ordbox, box, by="Dimension")
colnames(box) <- c("Dimension","Ordinal","Isomap","Interval")
prescriptionMelted <- reshape2::melt(box, id.var='Dimension')
box = merge(splinebox, box, by="Dimension")
colnames(box) <- c("Dimension","Spline","Ordinal","Isomap","Interval")
prescriptionMelted <- reshape2::melt(box, id.var='Dimension')
ggplot(prescriptionMelted, aes(x=Dimension, y=value, col=variable)) + geom_line()


#--------------------------------------Fixed Tree Framework---------------------------------#


# Evaluate using the data from the given tree
iso_fixed_tree     <- (test_algos_given_tree(tree,200,numA,"Isomap",numClusters,clusterIndexTree,clusterIndexA))
int_fixed_tree     <- (test_algos_given_tree(tree,200,numA,"interval",numClusters,clusterIndexTree,clusterIndexA))
ord_fixed_tree     <- (test_algos_given_tree(tree,200,numA,"ordinal",numClusters,clusterIndexTree,clusterIndexA))
spline_fixed_tree  <- (test_algos_given_tree(tree,200,numA,"mspline",numClusters,clusterIndexTree,clusterIndexA))

iso_fixed_score    <- mean(iso_fixed_tree)
int_fixed_score    <- mean(int_fixed_tree)
ord_fixed_score    <- mean(ord_fixed_tree)
spline_fixed_score <- mean(spline_fixed_tree)


# ------------------------------ Evaluate on random trees ------------------------------ #

# Try bigger trees with more antibodies and more centers the bigger the trees
treeSizes <- c(10,20,30,40,50,60,70,80,90,100)
numAs <- c(8,16,32,64,128)
maxClusters <- 8
iso_rand_scores <- matrix(0,length(treeSizes),maxClusters)
int_rand_scores <- matrix(0,length(treeSizes),maxClusters)
ord_rand_scores <- matrix(0,length(treeSizes),maxClusters)
spline_rand_scores <- matrix(0,length(treeSizes),maxClusters)

# Run for the defined vectors and use the evaluation scheme for each situation
# Also iterate over the number of clusters
for (clus in 1:maxClusters){
  for (size in 1:length(treeSizes)){
    m <- round(treeSizes[size]/clus)
    r <- treeSizes[size]%%m
    clusterIndexA <- rep(1:clus,m)
    clusterIndexA <- c(clusterIndexA,1:r)
    iso_rand_scores[size,clus] <- mean(test_algos_tree(100,numA,"Isomap",clus,treeSize,clusterIndexA))
    int_rand_scores[size,clus] <- mean(test_algos_tree(100,numA,"interval",clus,treeSize,clusterIndexA))
    ord_rand_scores[size,clus] <- mean(test_algos_tree(100,numA,"ordinal",clus,treeSize,clusterIndexA))
    spline_rand_scores[size,clus] <- mean(test_algos_tree(100,numA,"mspline",clus,treeSize,clusterIndexA))
    print(clus)
  }
}
# ----------------------------------------- Plotting ----------------------------------- #
# Two clusters:
myplots <- list()
for(i in 3:maxClusters){
  twoClusters <- cbind(treeSizes,iso_rand_scores[,i],int_rand_scores[,i],ord_rand_scores[,i],spline_rand_scores[,i])
  twoClusters <- as.data.frame(twoClusters)
  colnames(twoClusters) <- c("TreeSize", "Isomap", "Interval", "Ordinal", "Spline")
  twoClustersMelted <- reshape2::melt(twoClusters, id.var='TreeSize')
  myplots[[i-2]] <- ggplot(twoClustersMelted, aes(x=TreeSize, y=value, col=variable)) + geom_line()+ggtitle(paste('#Clusters',i))
}
grid.arrange(myplots[[1]],myplots[[2]],myplots[[3]],myplots[[4]],myplots[[5]],myplots[[6]],nrow=2)

#--------------------------------------------------------------------------------------------#

