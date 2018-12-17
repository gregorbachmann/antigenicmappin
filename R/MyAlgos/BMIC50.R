#############################################################################
# More natural way to simulate distances/clusters
#############################################################################

library(phytools)
# Read in Newick Tree coming from real viral data
tree <- read.newick(file="/Users/gregorbachmann/Desktop/ETH/Semester9/CompBioProject/Data/HIV_Newick_tree")
plot(tree)


create_OU_clusters <- function(tree,dim,numA,numClusters,clusterIndexTree,clusterIndexA){
  ############################################################################
  # Args:
  #     tree: Newick tree
  #      dim: Shape space dimension
  #     numA: Number of Antibodies (Number of Viruses is fixed in tree)
  # Returns:
  #           Non-symmetric distance matrix
  ############################################################################
  
  numV <- length(tree$tip.label)
  dist_matrix <- matrix(0,numA,numV)
  
  # To obtain IC-50 distances, run for each antibody a OU on the given tree
  # Try to create distances that correspond to the clusters in the tree
  # To this end we have to adapt the resulting OU trait values
  # Divide the number of antibodies by the number of clusters in the tree and
  # assign each group of antibodies to one cluster
  
  groupSize <- round(numA/numClusters)
  
  for (a in 1:numA){
    for(num in 1:numClusters){
        if (clusterIndexA[a]==num){
          # Save indices which belong to cluster num
          activeIndex <- which(clusterIndexTree==num)
          # Calculate trait values using OU
          trait <- multiOU(tree,alpha=2,1)
          minVal <- min(trait)
          # Take biggest distance in cluster num
          maxVal <- max(trait[activeIndex])
          # Make sure they have smallest distance to cluster num
          dist <- maxVal-minVal+abs(rnorm(1,mean(trait),1))
          trait[activeIndex] <- trait[activeIndex] - dist
          dist_matrix[a,] <- trait
        }
    }
  }
  # Make sure distance matrix only has positive entries
  dist_matrix <- dist_matrix + 2*abs(min(dist_matrix))  #- 2*min(dist_matrix)
  
  return(dist_matrix)
}


