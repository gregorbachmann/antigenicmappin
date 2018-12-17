# Example using the earth data
a =0
countries <- read.csv("http://egallic.fr/R/Blog/Cartes/countries_distances.csv",header=TRUE)
dist_matrix <-  matrix(0,186,186)

for (i in 1:186){
  for (j in 1:186){
    if (i == j){
      dist_matrix[i,j] <- 0
      a <- a+1
     
    }
    if (i != j){
      if(countries$dist[(i-1)*186+j-a]!=0){
        dist_matrix[i,j] <-  countries$dist[(i-1)*186+j-a]
      }
      if(countries$dist[(i-1)*186+j-a]==0)
        dist_matrix[i,j] <- 1
    }
  }  
}

dist_matrix$symmetric <- dist_matrix
dist_matrix$names <- matrix(0,186)
dist_matrix$names[1] <- "Canada"
dist_matrix$names[2:186] <- as.character(countries$pays2[1:185])
for (i in 1:186){
  if (i%%3==0){
    dist_matrix$names[i] <- ""
  }
}


dist_matrix$names <- matrix(0,186)
dist_matrix$names[1] <- "Canada"
dist_matrix$names[2:186] <- as.character(countries$pays2[1:185])
for (i in 1:186){
  if (i%%3==0){
    dist_matrix$names[i] <- ""
  }
}
fit_1 <- mdsvis(dist_matrix,2,init="torgerson","ratio",plot=FALSE)
gg <- ggplot(as.data.frame(fit_1$conf), aes(x = fit_1$conf[,1],y= fit_1$conf[,2],label = dist_matrix$names))
gg +geom_text(check_overlap = TRUE)  
fit <- (sammon(dist_matrix$symmetric,k=2))
plot(fit$points)

gg <- ggplot(as.data.frame(fit$points), aes(x = fit$points[,1],y= fit$points[,2],label = dist_matrix$names))
gg +geom_text(check_overlap = TRUE)  

library(tsne)
a <- tsne(data2$symmetric,k=2)
