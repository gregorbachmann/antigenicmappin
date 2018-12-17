####################################
######## Placing algorithms ########
####################################



######## placing algorithm using sammon()

placing.sammon<-function(dis,k=2,niter=1000,magic=0.2,tol=1e-4){
	
	# dis:			distance matrix of viruses and antibodies, dis<-read.csv("~/work/projects/mapping/testdatasets/virusantibodies/distances/sim11_nvi3_nab2_dist.csv")
	# k:			dimension of the resulting map, parameter for sammon()
	# niter:		maximum number of iterations, parameter for sammon()
	# magic:		initial value of the step size constant in diagonal Newton method, parameter for sammon()
	# tol:			tolerance for stopping, in units of stress, parameter for sammon()
	  
	 
	
	
	
	require(MASS)
	
	# number of virus strains
	nv<-length(unique(dis[,"virus"]))
	# number of antibodies
	na<-length(unique(dis[,"ab"]))
	
	# transform the distance matrix into a symmetric matrix
	sym.dis<-sym.mat.dis(dis)
		
	# find the coordinates, starting with a completely random matrix
	zw<-sammon(sym.dis,matrix(runif(k*(nv+na)),nrow=nv+na,ncol=k),niter=niter,magic=magic,tol=tol)
	
	coord<-data.frame(matrix(NA,nrow=nv+na,ncol=k+2,dimnames=list(rep("",nv+na),c("character","number",paste("x",1:k,sep="")))))
	coord[,"character"]<-c(rep("V",nv),rep("A",na))
	coord[,"number"]<-c(unique(dis[,"virus"]),unique(dis[,"ab"]))
	coord[,1:k+2]<-zw$points
	out<-list(stress=zw$stress,coordinates=coord)
	
	
}



# How does stress relate to random distance matrices?
# How does stress relate to the number of different distances, i.e. especially number of Ab and viruses?
# What does sammon() do?
# write a function to display points with circular shades around them.



#### displaying virus points with "uncertainty"

# Euclidean distance in L2
L2.dist<-function(x,y){sqrt(sum((x-y)^2))}


# function to determine a point that 1) lies on the line between two points a and v 2) and is dr from point a away 3) from the two possible solutions we want the one which is closer to v
farpoint<-function(a,v,dr){
	
	# a 	2d vector of coordinates of a
	# v  	2d vector of coordinates of v
	# dr 	real distance between a and the new point
	
	a1<-a[1]
	a2<-a[2]
	v1<-v[1]
	v2<-v[2]
	
	m<-(a2-v2)/(a1-v1)
	b<-(v2*a1 - v1*a2)/(a1-v1)
	c<-dr^2 - a1^2 - a2^2
	
	# caculating the two possible x values
	x.1<- - (b*m-a2*m-a1)/(1+m^2) + ( ( (b*m-a2*m-a1)/(1+m^2) )^2 -(b^2 - 2*a2*b -c)/(1+m^2) )^0.5
	x.2<- - (b*m-a2*m-a1)/(1+m^2) - ( ( (b*m-a2*m-a1)/(1+m^2) )^2 -(b^2 - 2*a2*b -c)/(1+m^2) )^0.5	
	# calculating the two possible y values
	
	y.1<- m*x.1 +b
	y.2<- m*x.2 +b
	
	nv.1 <- c(x.1,y.1)
	nv.2 <- c(x.2,y.2)
	
	# choose the point which is closer to v
	
	out<-nv.1
	if(L2.dist(nv.2,v)<L2.dist(nv.1,v)){
		out<-nv.2
	}
	
	out
}


# ordering points anti-clockwise in 2D around a point v
anti.clockwise<-function(v,xy){
	
	# v:		point around which the other points shall be ordered anti-clockwise
	# xy:		matrix with the other points
	
	# zw stores the angles between (v+(0,1)) and the vectors in the matrix rows
	
	ord<-c()
	ord2<-c()
	for(i in 1:(dim(xy)[1]) ){
			
		zw<-acos(sum( c(0,1)*(xy[i,]-v) ) / (sum(c(0,1)^2)^0.5 * sum((xy[i,]-v)^2)^0.5) ) 
		
		if( (xy[i,1]-v[1])>0 ){zw<-2*pi-zw}
		
		ord<-c(ord,zw) 	
	}
	
	xy[order(ord),]
}


plot.estcoord.uncertainty<-function(dis,k=2,niter=1000,magic=0.2,tol=1e-4,fol.out,fn.out,uncarea=TRUE){
	
	# dis:			distance matrix of viruses and antibodies, dis<-read.csv("~/work/projects/mapping/testdatasets/virusantibodies/distances/sim1_nvi10_nab5_dist.csv")
	# k:			dimension of the resulting map, parameter for sammon()
	# niter:		maximum number of iterations, parameter for sammon()
	# magic:		initial value of the step size constant in diagonal Newton method, parameter for sammon()
	# tol:			tolerance for stopping, in units of stress, parameter for sammon()
	# fol.out: 		the folder name where /coord/ and /figures/ is stored, eg. fol.out<-"~/work/projects/mapping/results/test_placing_unc/"
	# fn.out: 		file name eg fn.out<-"test_sim1"  
	# uncarea:		uncertainty areas drawn if uncarea==TRUE
	library(graphics)
	
	zw<-placing.sammon(dis,k,niter,magic,tol)
	coor<-zw$coordinates
	
	pdf(paste(fol.out,"figures/",fn.out,".pdf",sep=""),width=4,height=4)
	
		xlim<-range(coor[,"x1"])
		xlim[1]<-floor(xlim[1])
		xlim[2]<-ceiling(xlim[2])
		ylim<-range(coor[,"x2"])
		ylim[1]<-floor(ylim[1])
		ylim[2]<-ceiling(ylim[2])
		
		
		plot(c(0,1),c(0,1),xlim=xlim,ylim=ylim,type="n", xlab="antigenic coordinate 1",ylab="antigenic coordinate 2", main=paste("stress=",signif(zw$stress),sep=""))
		
		
		# adding stars around the virus points with polygon function
		
		# number of virus strains
		nv<-length(unique(dis[,"virus"]))
		# number of antibodies
		na<-length(unique(dis[,"ab"]))
		
		# vector with virus numbers
		v.v<-unique(coor[which(coor[,"character"]== "V" ) , "number"])
		# vector with antibody numbers
		v.a<-unique(coor[which(coor[,"character"]== "A" ) , "number"])
		
		# uncertainty areas are drawn firs
			if(uncarea){
		
				for(i in 1:nv){
			
				# determine the points of the stars 
				stars<-matrix(NA,nrow=0,ncol=2,dimnames=list(c(),c("x1","x2")))
				for(j in 1:na){
					stars<-rbind(stars, farpoint(as.numeric(coor[which( coor[,"character"]=="A" & coor[,"number"]==v.a[j] ), c("x1","x2")] ), as.numeric( coor[which( coor[,"character"]=="V" & coor[,"number"]==v.v[i] ), c("x1","x2")] ), dis[which(dis[,"virus"]==v.v[i] & dis[,"ab"]==v.a[j] ), "distance"]) )	
				}
				# order them anti-clockwise
				stars<-anti.clockwise(as.numeric( coor[which( coor[,"character"]=="V" & coor[,"number"]==v.v[i] ), c("x1","x2")] ),stars)
			
				# draw a polygon
				polygon(stars[,"x1"],stars[,"x2"],col=rgb(0,1,1,alpha=0.1),border=rgb(0,1,1),lwd=0.5)
			
			}
		}
		
		
		
		
		
		plot.coord(coor,xlim,ylim,add=TRUE)
		
		
	dev.off()
	
	write.table(coor,paste(fol.out,"coord/",fn.out,".csv",sep=""),quote=F,sep=",",row.names=F,col.names=T)
	
	return(zw)
	
}




