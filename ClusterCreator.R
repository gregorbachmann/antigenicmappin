#########################################
###Generate entire data clusters in 2D###
#########################################

# I'll create artificial data with numclusters in a box with each length equal to five
# This cluster creation method just uses raw randomness

create_test_clusters <- function(numclusters,dim,numpoints,n.v,plot=FALSE){
  ######################################################################################
  # Args:
  #     numclusters: Number of clusters desired
  #             dim: Dimension of data point
  #       numpoints: Vector of length numcluster containing the number of samples 
  #                  belonging to the clusters
  #            plot: TRUE if plot of clusters should be produced
  # Returns:
  #     List consisting of virus and antibody coordinates, cluster centers and memberships
  #     Plot of clusters if plot=TRUE
  ######################################################################################
  
  # Initialize random cluster centers in the box
  cluster_centers <- matrix(0,numclusters,dim)
  
  for (i in 1:numclusters){
    for (j in 1:dim){
      cluster_centers[i,j] <- runif(1,0,10) # /sqrt(dim)
    }
  }
  # Now scatter some datapoints randomly around the centers, each component 
  # being Gaussian with a small variance to guarantee closeness
  data = matrix(0,sum(numpoints),dim)
  membership <- matrix(0,sum(numpoints))
  k <- 1
  for(i in 1:numclusters){
    for (sample in 1:numpoints[i]){
      for(j in 1:dim){
        data[k,j] <- rnorm(1,cluster_centers[i,j],0.5)
        
      }
      membership[k] <- i
      k <- k+1
    }
  }
  
  # Now assign each data point randomly to virus or antibody
  virus_index <- sample(1:sum(numpoints),n.v,FALSE)
  memb_virus <- membership[virus_index]
  virus <- data[virus_index,]
  antibody <- data[-virus_index,]
  memb_anti <- membership[-virus_index]
  membership <- as.factor(c(memb_virus,memb_anti))
  
  # Plot the clusters with colors corresponding to assignments and form corresponding to
  # virus/antibody
  
  if (plot){
    all <- rbind(virus,antibody)
    shape <- as.factor(cbind(rep("Antibody",length(memb_anti)),rep("Viral Strain",length(memb_virus))))
    all <- as.data.frame(all)
    colnames(all) <- c("Axis1","Axis2")
    p <- ggplot(all, aes(x=Axis1, y=Axis2,color=membership,size=1.5,shape=shape)) + 
      geom_point()
    print(p)
  }
  return(list(virus=virus,antibody=antibody,clusters=cluster_centers,membership=membership))
}


# Calculates the two distance matrices out of the cluster data

preprocess <- function(cluster){
  # Args: 
  #     cluster: List of the form returned by create_test_clusters
  # Returns: 
  #     Two type of distance matrices in form of a list
  
  n.v <- dim(cluster$virus)[1]
  n.a <- dim(cluster$antibody)[1]
  coord <- rbind(cluster$antibody,cluster$virus)
  stand_dist <- calc_dist(coord,n.v,n.a)
  sym_dist <- sym_matrix(stand_dist)
  
  return(list(standard=stand_dist,symmetric=sym_dist))
}


###########################################################################################

# Use a cluster creation method that incorporates knowledge on genetic distance between
# the viral strains 

create_OU_clusters_given_tree <- function(tree,numA,numClusters,clusterIndexTree,clusterIndexA,plot=FALSE,w=numClusters,h=1){
  ############################################################################
  # Args:
  #                 tree: Newick tree
  #                 numA: Number of Antibodies (Number of Viruses is fixed in tree)
  #          numClusters: Number of clusters/clades in your tree
  #     clusterIndexTree: Cluster Assignment for each tip
  #        clusterIndexA: Cluster Assignment for each antibody
  #                 plot: Bolean indicating whether plot of trait values for clusters 
  #                       should be produced
  #                  h,w: Parameters passed to par(mfrow=c(w,h)) for plotting
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
  
  par(mfrow=c(w,h))
  colnames(dist_matrix) <- tree$tip.label
  cur <- 0
  if (plot){
    for(a in 1:numA){
      if (clusterIndexA[a]!=cur){
        dotTree(tree,dist_matrix[a,])
        cur = cur + 1
      }
    }
  }
  
  return(dist_matrix)
}


##########################################################################################

# Create clusters based on random trees and the Ornstein Uhlenbeck process

create_OU_clusters <- function(numClusters,treeSize,numA,clusterIndexA,plot=FALSE,w=numClusters,h=1){
  ########################################################################################
  # Args:
  #      numClusters: Number of clusters  
  #         treeSize: Number of tips of the random trees 
  #             numA: Number of antibodies
  #    clusterIndexA: Assignment of antibodies to clusters
  #             plot: Bolean indicating whether plot of trait values for clusters 
  #                   should be produced
  #              h,w: Parameters passed to par(mfrow=c(w,h)) for plotting
  # 
  # Returns:
  #         List of distance matrix and true cluster assignments
  ########################################################################################
  library(phytools)
  # Generate random tree
  tree <- rtree(treeSize)
  # Find "good" clusters in the tree 
  clusterIndexTree <- find_clusters_in_tree(tree,numClusters)
  
  numV <- length(tree$tip.label)
  dist_matrix <- matrix(0,numA,numV)
  # Run OU model on random tree for each antibody
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
  dist_matrix <- dist_matrix + 2*abs(min(dist_matrix))
  colnames(dist_matrix) <- tree$tip.label
  cur <- 0
  par(mfrow=c(w,h))
  if (plot){
    for(a in 1:numA){
      if (clusterIndexA[a]!=cur){
        dotTree(tree,dist_matrix[a,])
        cur = cur + 1
      }
    }
  }
  
  return(list(dist_matrix=dist_matrix,membership=clusterIndexTree))
}


find_clusters_in_tree <- function(tree,numClusters){
  tree <- force.ultrametric(tree)
  new <- as.hclust.phylo(tree)
  index <- cutree(new,numClusters)
  return(index)
}