# File to readin the virus-antibody data and get them into the right format for the algorithms
# The dataset should contain a column with the name distance and should be ordered w.r.t
# the viruses, so Virus1, then all antibodies and the distances, then Virus2 and all distances
# etc.


read_in <- function(link,n.v,n.a,ordering){
  # Input: link: string indicating the file location on your computer
  #               e.g. "/Users/gregorbachmann/Desktop/ETH/Semester9/CompBioProject/Data/CAP256/CAP256_bnabdata.csv"
  #         n.v: number of viruses in dataset
  #         n.a: number of antibodies in dataset
  #         ordering: 0 if ordered w.r.t. virus and 1 if ordered w.r.t. antibody
  # Output: A list consisting of the symmetrized distance matrix, the standard distance matrix
  #         and a vector with the corresponding names
  # Note: Cleaning is applied to remove all viruses with maximal distance (here 50) to all
  #       antibodies.
  
  data <- read.csv(link,header=TRUE)
  distances <- data$distance
  virusnames <- unique(as.character(data$virusname))
  abnames <- unique(as.character(data$abname))
  names <- c(abnames,virusnames)
  time <- matrix(1,n.v+n.a)
  col <-  matrix(1,n.v+n.a)
  old <-  0
  count <- 1 
  if (ordering==0){
    for (i in 1:n.v){
      if (data$virus_week[(i-1)*n.a+1]!=old){
        count = count + 1
      }
      col[n.a+i] <- count
      time[n.a+i] <- data$virus_week[(i-1)*n.a+1]
      old <- time[n.a+i]
    }
  }
  if (ordering==1){
    for (i in 1:n.v){
      if (data$virus_week[i]!=old){
        count = count + 1
      }
      col[n.a+i] <- count
      time[n.a+i] <- data$virus_week[i]
      old <- time[n.a+i]
    }
    for (i in 1:n.a){
      time[i] <- data$ab_week[(i-1)*n.v+1]
    }
  }
  # Matrix to store the standard distances
  #if (ordering == 0){
  #  stand_dist <- matrix(0,n.v,n.a)
  #  for (i in 1:n.v){
  #    stand_dist[i,] <- distances[((i-1)*n.a+1):(n.a*i)]
  #  }
  #}
  #if (ordering == 1)
  stand_dist <- matrix(0,n.a,n.v)
  for (i in 1:n.a){
    stand_dist[i,] <- distances[((i-1)*n.v+1):(n.v*i)]
      
  }
  
  # Clean:
  a <- c()
  for (i in 1:n.v){
    if (sum(stand_dist[,i])==50*n.a){
      a <- c(a,i)
     
    }
  }
  if (length(a)!=0){
    stand_dist <- stand_dist[,-a]
  }
  
  # Symmetrize resulting matrix using NA's
  sym <- sym_matrix(stand_dist)
  if (length(a)!=0){
    names <- names[-(a+n.a)]
    time <- time[-(a+n.a)]
    col <- col[-(a+n.a)]
  }
  return(list(standard=stand_dist,symmetric=sym,names=names,time=time,col=col))
}



