########################
##Entirely own Methods##
########################


sammon_cost <- function(unwrapped_coord,true_dist,dim,n.v,n.a){
  coord <- matrix(0,n.a+n.v,dim)
  for (i in 1:(n.a+n.v)){
    coord[i,] <- unwrapped_coord[((i-1)*(dim)+1):(i*dim)]
  }
  ab_matrix <- coord[1:n.a,]
  virus_matrix <- coord[(n.a+1):(n.a+n.v),]
  dist_matrix <- matrix(0,n.v,n.a)
  
  for (i in 1:length(virus_matrix[,1])){
    for (j in 1:length(ab_matrix[,1])){
      dist_matrix[i,j] <- sqrt(sum((virus_matrix[i,]-ab_matrix[j,])^2))
    }
  }
  return(sum(((true_dist-dist_matrix)^2)/true_dist))
}

derivative <- function(embedded_coord,true_dist,n.v,n.a,dim){
  embedded_coord <- matrix(embedded_coord,n.a+n.v,dim,byrow=TRUE)
  embedded_dist <- calc_dist(embedded_coord,n.v,n.a)
  A <- ((true_dist)-t(embedded_dist))/t(embedded_dist)
  A <- A/true_dist
  k <- dim
  gradientab <- matrix(0,n.a,k)
  for(i in 1:n.a){
    for(j in 1:k){
      B <- embedded_coord[(n.a+1):(n.a+n.v),j]-embedded_coord[i,j]
      gradientab[i,j] <- -2*sum(A[,i]*B)
    }
  }
  flattenab <- as.vector(t(gradientab))
  
  gradientvir <- matrix(0,n.v,k)
  for(i in 1:n.v){
    for(j in 1:k){
      B <- embedded_coord[1:n.a,j]-embedded_coord[n.a+i,j]
      gradientvir[i,j] <- 2*sum(A[i,]*B)
    }
  }
  flattenvir <- as.vector(t(gradientvir))
  
  gradient <- c(flattenab,flattenvir)
  return(gradient)
}

optimize <- function(init_config, true_dist, n.v, n.a,dim,method){
  library('Rcgmin')
  flatten <- as.vector(t(init_config))
  #opt <- Rcgmin(par=init_config, fn=cost_function,gr=derivative,true_dist=true_dist, n.v=n.v, n.a=n.a)
  #opt <- Rcgmin(par=init_config, fn=cost_function,true_dist=true_dist, n.v=n.v, n.a=n.a)
  opt <- optim(par=flatten, fn=sammon_cost,gr=derivative,true_dist=true_dist, n.v=n.v, n.a=n.a,dim=dim,method=method)
  return(opt)
}

mysammon <- function(data, dim, n.v, n.a, seed, label, plot=TRUE){
  library(plot3Drgl)
  library(plot3D)
  # Sample random starting configuration
  init_config <- generate_shape_space_coord(n.v, n.a, k = dim, seed)
  flatten <- as.vector(t(init_config))
  # True distance 
  true_dist <- t(data$standard)
  
  #opt <- gradient_descent(embedded_coord,true_dist,n.v,n.a,10,0.0001)
  opt <- optimize(flatten, true_dist, n.v, n.a,dim=dim,method="CG")
  embedding <- matrix(opt$par,n.a+n.v,dim)
  if (dim == 3 & plot == TRUE){
    scatter3D(embedding[,1],embedding[,2],embedding[,3],labels=data$names,colvar=1:length(data$symmetric[,1]),col=data$time)
    if (label == TRUE){
      text3D(embedding[,1],embedding[,2],embedding[,3],labels=data$names,add=TRUE)
    }
    plotrgl()
  }
  if (dim == 2 & plot == TRUE){
    scatter2D(embedding[,1],embedding[,2],colvar=1:length(data$symmetric[,1]),col=data$time)
    if (label == TRUE){
      text2D(embedding[,1],embedding[,2],labels=data$names,add=TRUE)
    }
  }
  return(opt)
}
