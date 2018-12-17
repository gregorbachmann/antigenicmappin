# Generate the classic "Swiss Roll" to illustrate Isomap

library('plot3D')
library("plot3Drgl")

N = 600
t <- seq(from=0, to=10,length.out=N)
#t = 3*pi/2 * (1 + 2*rnorm(N,sd=1))
h = 11 * rnorm(N)
X = cbind(t*cos(t), h, t*sin(t)) #+noise*rnorm(N,sd=1)
scatter3D(X[,1],X[,2],X[,3],theta=0,phi=10,colvar=seq(from=0, to=10,length.out=N))
plotrgl()
# Unrolled version
scatter2D(t,h,colvar=seq(from=0, to=10,length.out=N))

# Project a smaller version of it with Isomap

distances <- dist(X)
embedding <- isomap(distances,ndim=2,k=6)
scatter2D(embedding$points[,1],embedding$points[,2],colvar=seq(from=0, to=10,length.out=N))
