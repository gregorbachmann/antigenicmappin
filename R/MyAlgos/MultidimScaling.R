###############################
## Own implementation of MDS ##
###############################

# Some (maybe) useful functions:

# Important: I will always have the virus in rows and antibodies in columns for distance
# First generate antibody/virus shape space vectors, with arbitrary "true" underlying
# dimension k in order to have random initial points.

generate_shape_space_coord <- function(n.v, n.a, k=2){
  ###############################################################################
  # Args:
  #     n.v: number of viruses (integer)
  #     n.a: number of antibodies (integer)
  #       k: Underlying dimension (number of coordinates of vectors, usually 2) (integer)
  # Returns
  #       A list with two matrixes containing the coordinates of virus/antibodies
  ##################################################################################
  
  # Sample the coordinates from a standard normal
  virus_matrix <- matrix(0,n.v,k)
  for (i in 1:n.v){
    for (j in 1:k){
      virus_matrix[i,j] <- runif(1)
    }
  }
  
  ab_matrix <- matrix(0,n.a,k)
  for (i in 1:n.a){
    for (j in 1:k){
      ab_matrix[i,j] <- runif(1)
    }
  }
  
  return(rbind(ab_matrix,virus_matrix))
}

# Calculate distances between virus and antibody given their coordinates

calc_dist <- function(coord,n.v,n.a){
  ############################################################################
  # Args:
  #     virus_matrix: coordinate matrix of viruses
  #        ab_matrix: coordinate matrix of antibodies
  # Returns: 
  #     A matrix with the current pairwise distances (only between antibody and virus)
  #############################################################################
  
  ab_matrix <- coord[1:n.a,]
  virus_matrix <- coord[(n.a+1):(n.a+n.v),]
  dist_matrix <- matrix(0,n.a,n.v)
  
  for (i in 1:length(ab_matrix[,1])){
    for (j in 1:length(virus_matrix[,1])){
      dist_matrix[i,j] <- sqrt(sum((virus_matrix[j,]-ab_matrix[i,])^2))
    }
  }
  return(dist_matrix) 
}

# Function to evaluate the error/stress under current configuration

cost_function <- function(embedded_coord, true_dist, n.v, n.a, cost="quadratic"){
  #############################################################
  # Args: 
  #       true_dist: Experimentally measured distances (matrix)
  #   embedded_dist: distances of our embedded vectors (matrix)
  #            cost: cost measure (char)
  # Returns: 
  #       Current stress/error
  #############################################################
  
  embedded_dist <- calc_dist(embedded_coord,n.v,n.a)
  if(cost=="quadratic")
    return(sum((true_dist-t(embedded_dist))^2))
}

# In order to apply the "usual" algorithms, convert the distance matrix into a symmetric one
# and insert NA whenever we need the pairwise distance between antibodyd-antibody or
# virus-virus

sym_matrix <- function(dist_matrix){
  ####################################################################
  # Args:
  #   dist_matrix: Distance matrix of the form as returned by calc_dist
  # Returns:
  #   Symmetrized distance matrix
  ####################################################################
  
  n.v <- length(dist_matrix[1,])
  n.a <- length(dist_matrix[,1])
  sym <- matrix(NA,n.v+n.a,n.v+n.a)
  
  # Fill the diagonal with zeros
  for (i in 1:(n.v+n.a)){
    sym[i,i] <- 0
  }
  
  for (i in 1:n.a){
    sym[n.v+i,1:n.v] <- dist_matrix[i,]
    sym[1:n.v,n.v+i] <- dist_matrix[i,]
  }
  return(sym)
}

# In order to evaluate the performance of our chosen embeddings, we need random
# test data.

generate_test_data <- function(n.v,n.a,k){
  ###################################################################
  # Args:
  #   n.v: Number of viruses
  #   n.a: Number of antibodies
  #   k: Underlying shape space dimension
  # Returns:
  #   Symmetric distance matrix and the true shape space coordinates
  ###################################################################
  
  shape_space_coord <- generate_shape_space_coord(n.v,n.a,k)
  distance <- calc_dist(shape_space_coord$virus,shape_space_coord$ab)
  sym <- sym_matrix(distance)
  
  return(list(coordinates=shape_space_coord,symmetric=sym))
}