##########################
##Own written algorithms##
##########################

optimize <- function(init_config, true_dist, n.v, n.a,dim){
  library('Rcgmin')
  flatten <- as.vector(t(init_config))
  #opt <- Rcgmin(par=init_config, fn=cost_function,gr=derivative,true_dist=true_dist, n.v=n.v, n.a=n.a)
  #opt <- Rcgmin(par=init_config, fn=cost_function,true_dist=true_dist, n.v=n.v, n.a=n.a)
  opt <- Rcgmin(par=flatten, fn=new_cost_function,gr=derivative,true_dist=true_dist, n.v=n.v, n.a=n.a,dim=dim)
  return(opt)
}

derivative <- function(embedded_coord,true_dist,n.v,n.a,dim){
  embedded_dist <- calc_dist(embedded_coord,n.v,n.a)
  A <- (t(true_dist)-embedded_dist)/embedded_dist
  k <- dim
  gradientab <- matrix(0,n.a,k)
  for(i in 1:n.a){
    for(j in 1:k){
      B <- embedded_coord[(n.a+1):(n.a+n.v),j]-embedded_coord[i,j]
      gradientab[i,j] <- -2*sum(A[,j]*B)
    }
  }
  gradientvir <- matrix(0,n.v,k)
  for(i in 1:n.v){
    for(j in 1:k){
      B <- embedded_coord[1:n.a,j]-embedded_coord[n.a+i,j]
      gradientvir[i,j] <- 2*sum(A[j,]*B)
    }
  }
  gradient <- c(gradientab,gradientvir)
  return(gradient)
}

gradient_descent <- function(embedded_coord,true_dist,n.v,n.a,iteration,step){
  conf <- embedded_coord
  for (i in 1:iteration){
    gradient <- derivative(conf,true_dist,n.v,n.a)
    conf = conf - step*gradient
    cost <- cost_function(conf,true_dist,n.v,n.a)
    print(cost)
  }
  return(conf)
}
myscaling <- function(data, dim, n.v, n.a, seed, label, plot=TRUE){
  library(plot3Drgl)
  library(plot3D)
  # Sample random starting configuration
  init_config <- generate_shape_space_coord(n.v, n.a, k = dim, seed)
  
  # True distance 
  true_dist <- data$standard
  
  #opt <- gradient_descent(embedded_coord,true_dist,n.v,n.a,10,0.0001)
  opt <- optimize(init_config, true_dist, n.v, n.a,dim=dim)
  if (dim == 3 & plot == TRUE){
    scatter3D(opt$par[,1],opt$par[,2],opt$par[,3],labels=data$names,colvar=1:length(data$symmetric[,1]),col=data$time)
    if (label == TRUE){
      text3D(opt$par[,1],opt$par[,2],opt$par[,3],labels=data$names,add=TRUE)
    }
    plotrgl()
  }
  if (dim == 2 & plot == TRUE){
    scatter2D(opt$par[,1],opt$par[,2],colvar=1:length(data$symmetric[,1]),col=data$time)
    if (label == TRUE){
      text2D(opt$par[,1],opt$par[,2],labels=data$names,add=TRUE)
    }
  }
  return(opt)
}