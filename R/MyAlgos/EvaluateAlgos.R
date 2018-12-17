############################################################
## Several Methods to evaluate Performance
########################################################

# I will test the derived/used algorithms on their capability to detect clusters
# given clusters in high-dimensional spaces. For simplicity I will use k-means on the 
# embedded data and see how many groups are correctly identified

library(clusterCrit)

test_algos_artifical <- function(dim,numTests,numPoints,method){
  #####################################################################################
  # Args:
  #         dim: Dimension of true shape space 
  #    numTests: Number of tests you want to run
  #   numPoints: Number of virus and antibodies (summed up) you want to generate
  #      method: String containing the name of the method you want to test
  #              At the moment one of 
  #             c("mspline","ordinal","interval","ratio","Isomap","Sammon")
  # Returns:
  #       Vector of length numtests containing for each run the Czekanowski_Dice score of 
  #       each corresponding clustering
  #####################################################################################
  
  # Initialize random parameters for create_test_clusters
  # Initialize different number of clusters (from 1 to 20)
  numClusters <- sample(2:7,numTests,replace=TRUE)
  clusterSize <- matrix(0,numTests,20)
  n.v <- as.integer(numPoints/2)
  
  for (index in 1:numTests){
    # On average equally distributed clusters
    n = as.integer(numPoints/numClusters[index])
    current <- sample(1:n,numClusters[index]-1, replace=TRUE)
    # Enforce total number of numpoints
    current <- c(current,numPoints-sum(current))
    clusterSize[index,1:numClusters[index]] <- current
  }
  
  fit <- matrix(0,numTests)
  for (run in 1:numTests){
    # Generate random data
    clusters <- create_test_clusters(numClusters[run],dim, clusterSize[run,1:numClusters[run]],n.v)
    # The true clustering
    true_assignments <- as.integer(clusters$membership)
    matrices <- preprocess(clusters)
    # Use estimator from Loewe book
    if (is.element(method,c("mspline","ordinal","interval","ratio"))){
      embedding <- mdsvis(matrices,2,method,"torgerson",FALSE,FALSE)$embedding
    }
    if (method=="Isomap"){
      embedding <- isomapvis(matrices, 2, 25, plot=FALSE,label=FALSE)$embedding
    }
    if (method=="Sammon"){
      embedding <- sammon_mapping(matrices$symmetric,2,FALSE,FALSE,n.v,numPoints-n.v)
      print(embedding$stress)
    }
    # Estimate the clustering in the low dimensional space using kmeans
    assignments <- as.matrix(cluster(embedding,numclusters[run],FALSE,method)$cluster)
    # Evaluate clustering using the Czekanowski Criterion
    fit[run] <- extCriteria(as.integer(assignments),true_assignments,"Czekanowski_Dice")$czekanowski_dice
    #print(run)
  }
  print("done")
  return(fit)
}

# To assess how much error is introduced by kmeans I run the same simulation but without projecting
# the data down but just clustering immediately on the cluster data generated

assess_kmeans_error <- function(dim,numTests,numPoints){
  # Args:
  #         dim: Dimension of true shape space 
  #    numTests: Number of tests you want to run
  #   numPoints: Number of virus and antibodies (summed up) you want to generate
  # Returns:
  #       Vector of length numtests containing for each run the Czekanowski_Dice score of 
  #       each corresponding kmeans clustering
  # Initialize random parameters for create_test_clusters
  # Initialize different number of clusters (from 1 to 20)
  numClusters <- sample(2:7,numTests,replace=TRUE)
  clusterSize <- matrix(0,numTests,20)
  n.v <- as.integer(numPoints/2)
  
  for (index in 1:numTests){
    # On average equally distributed clusters
    n = as.integer(numPoints/numClusters[index])
    current <- sample(1:n,numClusters[index]-1, replace=TRUE)
    # Enforce total number of numpoints
    current <- c(current,numPoints-sum(current))
    clusterSize[index,1:numClusters[index]] <- current
  }
  
  fit <- matrix(0,numTests)
  for (run in 1:numTests){
    # Generate random data
    clusters <- create_test_clusters(numClusters[run],dim, clusterSize[run,1:numClusters[run]],n.v)
    data <- rbind(clusters$virus,clusters$antibody)
    # The true clustering
    true_assignments <- as.integer(clusters$membership)
    # Estimate the clustering in the low dimensional space 
    assignments <- as.matrix(kmeans(data,numClusters[index])$cluster)
    # Evaluate clustering using the Czekanowski Criterion
    fit[run] <- extCriteria(as.integer(assignments),true_assignments,"Czekanowski_Dice")$czekanowski_dice
    print(run)
  }
  return(fit)
}


test_algos_given_tree <- function(tree,numTests,numA,method,numClusters,clusterIndexTree,clusterIndexA){
  ######################################################################################
  # Args:
  #      tree: Phylotree in right format (from phylotools)
  #   numTest: Number of tests to evaluate
  #      numA: Number of antibodies
  #    method: Method you want to evaluate
  # Returns:
  #   Score of the method
  #######################################################################################
  
  fit <- matrix(0,numTests)
  for (run in 1:numTests){
    dist_matrix <- create_OU_clusters_given_tree(tree,numA,numClusters,clusterIndexTree,clusterIndexA)
    sym_dist <- sym_matrix(dist_matrix)
    matrices <- list(stand_dist=dist_matrix,symmetric=sym_dist)
    
    if (is.element(method,c("mspline","ordinal","interval","ratio"))){
      embedding <- mdsvis(matrices,2,method,"torgerson",FALSE,FALSE)$embedding
    }
    if (method=="Isomap"){
      embedding <- isomapvis(matrices, 2, 25, plot=FALSE,label=FALSE)$embedding
    }
    if (method=="Sammon"){
      embedding <- sammon_mapping(matrices$symmetric,2,FALSE,FALSE,n.v,numpoints-n.v)
      print(embedding$stress)
    }
    
    true_assignments <- as.integer(clusterIndexTree)
    assignments <- as.matrix(cluster(embedding,numClusters,FALSE,method)$cluster)
    assignments <- assignments[1:(length(assignments)-numA)]
    # Evaluate clustering using the Czekanowski Criterion
    fit[run] <- extCriteria(as.integer(assignments),true_assignments,"Czekanowski_Dice")$czekanowski_dice
    
  }
  return(fit)
}


###########################################################################################

test_algos_tree <- function(numTests,numA,method,numClusters,treeSize,clusterIndexA){
  ######################################################################################
  # Args:
  #          numTest: Number of tests to evaluate
  #             numA: Number of antibodies
  #           method: Method you want to evaluate
  #      numClusters:
  #         treeSize:
  #    clusterIndexA:
  # Returns:
  #   Score of the method
  #######################################################################################
  
  fit <- matrix(0,numTests)
  for (run in 1:numTests){
    # Sample random tree and find clusters in it
    stuff <- create_OU_clusters(numClusters,treeSize,numA,clusterIndexA)
    dist_matrix <- stuff$dist_matrix
    clusterIndexTree <- stuff$membership
    # Symmetrize distance matrix
    sym_dist <- sym_matrix(dist_matrix)
    matrices <- list(stand_dist=dist_matrix,symmetric=sym_dist)
    
    # Calculate resulting embedding 
    if (is.element(method,c("mspline","ordinal","interval","ratio"))){
      embedding <- mdsvis(matrices,2,method,"torgerson",FALSE,FALSE)$embedding
    }
    if (method=="Isomap"){
      embedding <- isomapvis(matrices, 2, 25, plot=FALSE,label=FALSE)$embedding
    }
    if (method=="Sammon"){
      embedding <- sammon_mapping(matrices$symmetric,2,FALSE,FALSE,n.v,numpoints-n.v)
      print(embedding$stress)
    }
    
    true_assignments <- as.integer(clusterIndexTree)
    assignments <- as.matrix(cluster(embedding,numClusters,FALSE,method)$cluster)
    assignments <- assignments[1:(length(assignments)-numA)]
    # Evaluate clustering using the Czekanowski Criterion
    fit[run] <- extCriteria(as.integer(assignments),true_assignments,"Czekanowski_Dice")$czekanowski_dice
    #if (run%%1==0){
    #  print(run)
    #}
  }
  return(fit)
}