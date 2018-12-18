#########################################
##Embedding functions and visualization##
#########################################


isomapvis <- function(data, dim, nn, plot=TRUE,label=FALSE){
  ##############################################################################
  # Args:
  #         data:   A list containing the symmetric matrix with NA's and the
  #                 names of viruses and antibodies (Result of the function ReadIn for instance)
  #          dim:    Desired dimension of the embedding, if 2 or 3 a plot can be made
  #           nn:     Number of Nearest Neighbours used in Isomap
  #         plot:   Either TRUE or FALSE
  #        label:  Either TRUE or FALSE
  # Returns:
  #         Embedding found via Isomap
  ################################################################################

  library('vegan')
  library('plot3D')
  library("plot3Drgl")
  embedding <- isomap((data$symmetric),ndim=dim,k=nn)
  Label <- as.factor(c(rep("Antibody",length(data$standard[,1])),rep("Strain",length(data$standard[1,]))))
  IsolationWeek <- as.factor(data$time)
  # 2d version
  if (dim == 2){
    emb <- as.data.frame(embedding$points)
    colnames(emb) <- c("Axis1","Axis2")
    p <- ggplot(emb, aes(x=Axis1, y=Axis2,label=data$names,shape=Label,col=IsolationWeek)) + geom_point(size=2)+ggtitle("Isomap Embedding")+
         theme(plot.title = element_text(hjust = 0.5))
    if (label == TRUE){
      p <- p + geom_label(vjust = 0) + geom_point(size=2)
    }
    if (plot==TRUE){
      print(p)
    }
  }
  # 3d Version
  if (dim == 3 & plot == TRUE){
    scatter3D(embedding$points[,1],embedding$points[,2],embedding$points[,3],labels=data$names,colvar=1:length(data$symmetric[,1]),col=data$col)
    if(label==TRUE){
      text3D(embedding$points[,1],embedding$points[,2],embedding$points[,3],labels=data$names,add=TRUE)
    }
    plotrgl()
  }
  if (dim==2){
    return(list(embedding=embedding,plot=p))
  }
  return(embedding)
}


mdsvis <- function(data, dim, method, init="random", plot=TRUE,label=TRUE,real){
  ################################################################################
  # Args:
  #     data:   A list containing the symmetric matrix with NA's and the names of viruses and antibodies
  #             (Result of the function ReadIn for instance)
  #      dim:    Desired dimension of the embedding, if 2 or 3 a plot is made
  #   method: Either "interval", "ratio", "ordinal" or "mspline"
  #     init:   Starting configuration
  #     plot:   Either TRUE or FALSE
  #     label:  Either TRUE or FALSE
  #     real:   Either TRUE or FALSE, indicating whether data has names and time attribut
  # Returns:
  #     Embedding found via SmaCof
  ###############################################################################

  library('smacof')
  library('plot3D')
  library("plot3Drgl")

  n <- length(data$symmetric[1,])
  # Need to specifiy a weight matrix to ignore the unknown distances
  weights <- matrix(1,n,n)
  for (i in 1:n){
    for (j in 1:n){
      if (is.na(data$symmetric[i,j])){
      #if( data$symmetric[i,j]==100){
        weights[i,j] <- 0
      }
      else weights[i,j] <- 1
    }
  }

  embedding <- smacofSym(data$symmetric,dim,method,weightmat=weights,init=init)
  Label <- as.factor(c(rep("Antibody",length(data$standard[,1])),rep("Strain",length(data$standard[1,]))))
  IsolationWeek <- as.factor(data$time)

  # 2d version
  if (dim == 2 && real==TRUE){
    emb <- as.data.frame(embedding$conf)
    colnames(emb) <- c("Axis1","Axis2")
    p <- ggplot(emb, aes(x=Axis1, y=Axis2,label=data$names,shape=Label,col=IsolationWeek)) + geom_point(size=2)+ggtitle(paste(method,"Embedding"))+
      theme(plot.title = element_text(hjust = 0.5))
    if (label == TRUE){
      p <- p + geom_label(vjust = 0) + geom_point(size=2)
    }
    if (plot == TRUE){
      print(p)
    }
  }
  if (dim==2 && real==FALSE){
    emb <- as.data.frame(embedding$conf)
    colnames(emb) <- c("Axis1","Axis2")
    p <- ggplot(emb, aes(x=Axis1, y=Axis2,shape=Label)) + geom_point(size=2)+ggtitle(paste(method,"Embedding"))+
      theme(plot.title = element_text(hjust = 0.5))
  }

  # 3d version
  if (dim == 3 & plot == TRUE){
    scatter3D(embedding$conf[,1],embedding$conf[,2],embedding$conf[,3],labels=data$names,colvar=1:length(data$symmetric[,1]),col=data$col)
    if (label == TRUE){
      text3D(embedding$conf[,1],embedding$conf[,2],embedding$conf[,3],labels=data$names,add=TRUE)

    }
    plotrgl()
  }
  if (dim==2){
    return(list(embedding=embedding,plot=p))
  }
  return(embedding)
}

# Use Sammon Mapping from package "MASS"

sammon_mapping <- function(data,dim,plot=TRUE,label=FALSE,n.v,n.a){
  # Args:
  #    data: List containing the symmetric and standard distance matrix and labels
  #     dim: Dimension you want to project down
  #    plot: Boolean indicating whether plot should be produced
  #   label: Boolean indicating whether labels should be added to plot
  #     n.v: Number of viruses in data
  #     n.a: Number of antibodies
  # Returns:
  #     Sammon Embedding
  ######################################################################################

  numruns <- 10
  best_stress <- Inf
  best_embedding <- 0

  # As algorithm gets stuck in local minimas very often, iterate 10 times:
  for (run in 1:numruns){
    # Get starting configuration (cmdscale doesn't accept NA's)
    init <-generate_shape_space_coord(n.v,n.a,2)
    #  Calculate embedding
    embedding <- sammon(matrices$symmetric,init , k = 2, niter = 300, trace = FALSE,magic=0.01,
                      tol = 1e-6)
    if (embedding$stress<best_stress){
      best_embedding <- embedding
      best_stress <- embedding$stress
    }
  }

   # 2D version
  if (dim==2 && plot==TRUE){
    plot(best_embedding$points)

    if (label == TRUE){
      text2D(best_embedding$points[,1],best_embedding$points[,2],labels=data$names,add=TRUE)
    }
  }

  # 3D version
  if (dim == 3 & plot == TRUE){
    scatter3D(best_embedding$points[,1],best_embedding$points[,2],best_embedding$points[,3])
    if (label == TRUE){
      text3D(best_embedding$points[,1],best_embedding$points[,2],best_embedding$points[,3],labels=data$names,add=TRUE)
    }
    plotrgl()
  }
  return(best_embedding)
}

cluster <- function(embedding,k,plot=TRUE,method){
  ########################################################################
  # Args:
  #     embedding: 2d coordinates for embedding
  #             k: number of clusters
  #        method: Method to cluster
  #                (one of c("mspline","ordinal","interval","ratio","Isomap","Sammon"))
  # Returns:
  #     kmeans object
  #     Plot with clusters visualized
  ########################################################################

  if (is.element(method,c("mspline","ordinal","interval","ratio"))){
    # Use kmeans with 10 initializations to avoid bad clusters
    clus <- kmeans(embedding$conf,k,nstart=10)
  }
  if (is.element(method,c("Isomap","Sammon"))){
    # Use kmeans with 10 initializations to avoid bad clusters
    clus <- kmeans(embedding$points,k,nstart=10)
  }
  if(plot){
    emb <- as.data.frame(embedding)
    rownames(emb) <- c("Axis1","Axis2")
    plot(embedding,col=clus$cluster)
    p <- ggplot(embedding, aes(x=Axis1, y=Axis2,col=as.factor(clus$cluster))) + geom_point()
  return(clus)
  }
}
