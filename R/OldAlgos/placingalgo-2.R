####################################
######## Placing algorithms ########
####################################




### definition of the error function



error.2D<-function(v.coord,n.v,n.a,dis){
	
	# v.coord:	vector with all coordinates of the virus and antibody points
	# n.v:		number of viruses
	# n.a: 		number of antibodies
	# dis:		distance matrix, e.g. dis<-read.csv("~/work/projects/mapping/testdatasets/virusantibodies/distances/sim11_nvi3_nab2_dist.csv")
	
	
	
	## first separate the coordinates vector into 
	# v.v1:		vector with the x-coordinates of the virus points
	# v.v2:		vector with the y-coordinates of the virus points
	# v.a1:		vector with the x-coordinates of the antibody points
	# v.a2: 	vector with the y-coordinates of the antibody points
	
	
	v.v1<-v.coord[1:n.v]
	v.v2<-v.coord[1:n.v + n.v]
	v.a1<-v.coord[2*n.v + 1:n.a]	
	v.a2<-v.coord[2*n.v + n.a + 1:n.a]	
	
	
	sum((dis[,"distance"]-((v.v1[dis[,"virus"]]-v.a1[dis[,"ab"]])^2 + (v.v2[dis[,"virus"]]-v.a2[dis[,"ab"]])^2 )^0.5)^2)
	
}


error.2D.1Ab1Virusfix<-function(v.cm1,c.v1,c.a1,n.v,n.a,dis){
	
	# v.cm1:	vector with all the coordinates of all viruses and antibodies except one virus and one antibody (virus 1 and ab 1 in the distance matrix)
	# c.v1:		coordinates of virus 1
	# c.a1: 	coordinates of antibody 1
	# n.v:		number of viruses
	# n.a: 		number of antibodies
	# dis:		distance matrix, e.g. dis<-read.csv("~/work/projects/mapping/testdatasets/virusantibodies/distances/sim11_nvi3_nab2_dist.csv")
	
	
	v.coord<-as.numeric(c(c.v1[1] , v.cm1[1:(n.v-1)], c.v1[2] ,v.cm1[(n.v-1)+1:(n.v-1)] , c.a1[1] , v.cm1[2*(n.v-1) + 1:(n.a-1) ] , c.a1[2],  v.cm1[2*(n.v-1) + (n.a-1) + 1:(n.a-1) ] ))
	
	out<-error.2D(v.coord,n.v,n.a,dis)
	out
}



### function to obtain the coordinate matrix out of the coordinates of virus 1 and antibody 1 and the remaining coordinate vector:

coordinate.matrix<-function(v.cm1,c.v1,c.a1,n.v,n.a){
	
	# v.cm1:	vector with all the coordinates of all viruses and antibodies except one virus and one antibody (virus 1 and ab 1 in the distance matrix)
	# c.v1:		coordinates of virus 1
	# c.a1: 	coordinates of antibody 1
	# n.v:		number of viruses
	# n.a: 		number of antibodies
	
	
	coo<-data.frame(matrix(NA,nrow=n.v+n.a,ncol=4,dimnames=list(rep("",n.v+n.a),c("character","number","x1","x2"))))
	
	coo[,"character"]<-c(rep("V",n.v),rep("A",n.a))
	coo[,"number"]<-c(1:n.v,1:n.a)
	coo[1,c("x1","x2")]<-as.numeric(c.v1)
	coo[n.v+1,c("x1","x2")]<-as.numeric(c.a1)
	
	coo[2:n.v,"x1"]<-v.cm1[1:(n.v-1)]
	coo[2:n.v,"x2"]<-v.cm1[(n.v-1)+1:(n.v-1)]
	coo[n.v+2:n.a,"x1"]<-v.cm1[2*(n.v-1)+1:(n.a-1)]
	coo[n.v+2:n.a,"x2"]<-v.cm1[2*(n.v-1) + n.a -1 + 1:(n.a-1)]
	
	coo
}


### function to calculate the difference between the data points in a coordinate.matrix as defined above
dist.v.Ab<-function(coord){
	
	# coord			coordinate matrix with type virus (V) and antibody (A)
	
	out<-matrix(NA,nrow=0,ncol=3,dimnames=list(c(),c("virus","ab","distance")))
	
	for(i in coord[which(coord[,"character"] =="V" ),"number"] ){
		for(j in coord[which(coord[,"character"] =="A" ),"number"] ){
			
			out<-rbind(out,c(i,j, sqrt( sum( (coord[which(coord[,"character"]=="V" & coord[,"number"]==i ),c("x1","x2")] - coord[which(coord[,"character"]=="A" & coord[,"number"]==j),c("x1","x2")] )^2 ) ) ))
			
		}
	}
	out
}



### function to plot a coordinate matrix

plot.coord<-function(coord,xlim,ylim,add=FALSE){
	
	if(!add){
		plot(xlim,ylim,type="n",xlab="x1",ylab="x2",xlim=xlim,ylim=ylim)
	}
	points(coord[which(coord[,"character"]=="V"),"x1"],coord[which(coord[,"character"]=="V"),"x2"],pch=16 ,cex=0.5)
	text(coord[which(coord[,"character"]=="V"),"x1"],coord[which(coord[,"character"]=="V"),"x2"], paste("V",coord[which(coord[,"character"]=="V"),"number"],sep=""), cex=0.7,pos=1)	
	
	points(coord[which(coord[,"character"]=="A"),"x1"],coord[which(coord[,"character"]=="A"),"x2"],pch=16,col=3 ,cex=0.5)
	text(coord[which(coord[,"character"]=="A"),"x1"],coord[which(coord[,"character"]=="A"),"x2"], paste("A",coord[which(coord[,"character"]=="A"),"number"],sep=""), cex=0.7,pos=1,col=3)
}



### placing function based on optim

placing.0<-function(dis, start, meth){
	
	# dis			matrix with distances, eg. dis<-read.csv("~/work/projects/mapping/testdatasets/virusantibodies/distances/sim11_nvi3_nab2_dist.csv")
	# start			matrix with starting values for viruses and antibodies: "character, number, x1, x2"; virus 1 and antibody 1 will be assumed to be known, or determined through their distance, eg start<-read.csv("~/work/projects/mapping/testdatasets/virusantibodies/coordinates/sim11_nvi3_nab2_coord.csv")
	# 
	
	
	## starting values for virus and antibody
	# virus 1
	c.v1<-as.numeric(start[which(start[,"character"]=="V" & start[,"number"]==1),c("x1","x2")])
	c.a1<-as.numeric(start[which(start[,"character"]=="A" & start[,"number"]==1),c("x1","x2")])

	sviab<-c(start[which(start[,"character"]=="V" & start[,"number"]>1),"x1"] ,start[which(start[,"character"]=="V" & start[,"number"]>1),"x2"], 	start[which(start[,"character"]=="A" & start[,"number"]>1),"x1"], start[which(start[,"character"]=="A" & start[,"number"]>1),"x2"] )
	
		
	n.v<- max(as.numeric(start[which(start[,"character"]=="V"),"number"]))
	n.a<- max(as.numeric(start[which(start[,"character"]=="A"),"number"]))
	
	# Nelder-Mead
	try(zw<-optim(sviab,error.2D.1Ab1Virusfix,c.v1=c.v1,c.a1=c.a1,n.v=n.v,n.a=n.a,dis=dis,control=list(maxit=1e6,reltol=1e-20),method=meth))
		
	
	
	if(zw$convergence == 0){
		print("convergence successful")
		print(zw)
		
		czw<-zw$par
		#matrix with new coordinates
		nc<-start
		nc[which(nc[,"character"]=="V" & nc[,"number"]>1),"x1"]<- czw[1:(n.v-1)]
		nc[which(nc[,"character"]=="V" & nc[,"number"]>1),"x2"]<- czw[(n.v-1)+1:(n.v-1)]
		nc[which(nc[,"character"]=="A" & nc[,"number"]>1),"x1"]<- czw[2*(n.v-1)+1:(n.a-1)]
		nc[which(nc[,"character"]=="A" & nc[,"number"]>1),"x2"]<- czw[2*(n.v-1) + (n.a-1) + 1:(n.a-1)]
		
		return(nc)
	}
	else{
		print("NO convergence")
		print(zw)
		czw<-zw$par
		#matrix with new coordinates
		nc<-start
		nc[which(nc[,"character"]=="V" & nc[,"number"]>1),"x1"]<- czw[1:(n.v-1)]
		nc[which(nc[,"character"]=="V" & nc[,"number"]>1),"x2"]<- czw[(n.v-1)+1:(n.v-1)]
		nc[which(nc[,"character"]=="A" & nc[,"number"]>1),"x1"]<- czw[2*(n.v-1)+1:(n.a-1)]
		nc[which(nc[,"character"]=="A" & nc[,"number"]>1),"x2"]<- czw[2*(n.v-1) + (n.a-1) + 1:(n.a-1)]
		
		return(nc)
	}
	
}



### placing function based on optimx

placing.1<-function(dis, start){
	
	# dis			matrix with distances, eg. dis<-read.csv("~/work/projects/mapping/testdatasets/virusantibodies/distances/sim11_nvi3_nab2_dist.csv")
	# start			matrix with starting values for viruses and antibodies: "character, number, x1, x2"; virus 1 and antibody 1 will be assumed to be known, or determined through their distance, eg start<-read.csv("~/work/projects/mapping/testdatasets/virusantibodies/coordinates/sim11_nvi3_nab2_coord.csv")
	# 
	
	
	## starting values for virus and antibody
	# virus 1
	c.v1<-as.numeric(start[which(start[,"character"]=="V" & start[,"number"]==1),c("x1","x2")])
	c.a1<-as.numeric(start[which(start[,"character"]=="A" & start[,"number"]==1),c("x1","x2")])

	sviab<-c(start[which(start[,"character"]=="V" & start[,"number"]>1),"x1"] ,start[which(start[,"character"]=="V" & start[,"number"]>1),"x2"], 	start[which(start[,"character"]=="A" & start[,"number"]>1),"x1"], start[which(start[,"character"]=="A" & start[,"number"]>1),"x2"] )
	
		
	n.v<- max(as.numeric(start[which(start[,"character"]=="V"),"number"]))
	n.a<- max(as.numeric(start[which(start[,"character"]=="A"),"number"]))
	

	try(zw<-optimx(sviab,error.2D.1Ab1Virusfix,c.v1=c.v1,c.a1=c.a1,n.v=n.v,n.a=n.a,dis=dis,control=list(maxit=1e6,reltol=1e-20)))
		
	
	
	if(zw$convergence == 0){
		print("convergence successful")
		print(zw)
		
		czw<-zw$par
		#matrix with new coordinates
		nc<-start
		nc[which(nc[,"character"]=="V" & nc[,"number"]>1),"x1"]<- czw[1:(n.v-1)]
		nc[which(nc[,"character"]=="V" & nc[,"number"]>1),"x2"]<- czw[(n.v-1)+1:(n.v-1)]
		nc[which(nc[,"character"]=="A" & nc[,"number"]>1),"x1"]<- czw[2*(n.v-1)+1:(n.a-1)]
		nc[which(nc[,"character"]=="A" & nc[,"number"]>1),"x2"]<- czw[2*(n.v-1) + (n.a-1) + 1:(n.a-1)]
		
		return(nc)
	}
	else{
		print("NO convergence")
		print(zw)
		czw<-zw$par
		#matrix with new coordinates
		nc<-start
		nc[which(nc[,"character"]=="V" & nc[,"number"]>1),"x1"]<- czw[1:(n.v-1)]
		nc[which(nc[,"character"]=="V" & nc[,"number"]>1),"x2"]<- czw[(n.v-1)+1:(n.v-1)]
		nc[which(nc[,"character"]=="A" & nc[,"number"]>1),"x1"]<- czw[2*(n.v-1)+1:(n.a-1)]
		nc[which(nc[,"character"]=="A" & nc[,"number"]>1),"x2"]<- czw[2*(n.v-1) + (n.a-1) + 1:(n.a-1)]
		
		return(nc)
	}
	
}





### when creating a distance matrix out of known coordinates, the ab/ab and virus/virus pairs should be omitted
# which positions in the distance matrix (object distance, created with dist() ) need to be set to NA (all virus/virus and ab/ab pairs)? 
# v.NA does this job
v.NA<-function(nv,na){
	# number of virions nv<-10
	# number of antibodies na<-5
	
	# first: add the virus/virus postions
	zw<-c()
	j<-0
	for(i in 1:(nv-1)){
		zw<-c(zw,j+1:(nv-i))
		j<-max(zw)+na
	}
	
	# second: add the ab/ab positions
	zw<-c(zw,j+na+1:( na*(na-1)/2 ))
	
	# return the positions
	zw
}





### function to create a symmetric matrix out of the distance matrix of structure "virus", "ab", "distance"
sym.mat.dis<-function(dis){
	
	# dis:			distance matrix of viruses and antibodies, dis<-read.csv("~/work/projects/mapping/testdatasets/virusantibodies/distances/sim11_nvi3_nab2_dist.csv")
	
	# number of virus strains
	nv<-length(unique(dis[,"virus"]))
	# number of antibodies
	na<-length(unique(dis[,"ab"]))
	
	out<-matrix(NA,nrow=nv+na,ncol=nv+na)
	
	# all diagonal elements are 0
	for(i in 1:(nv+na)){
		out[i,i]<-0
	}
	
	# distances for Ab and viruses added to the matrix
	for(i in 1:nv){
		
		for(j in 1:na){
			
			out[i,j+nv]<-dis[which(dis[,"virus"]==unique(dis[,"virus"])[i] & dis[,"ab"]==unique(dis[,"ab"])[j]), "distance"]
			out[j+nv,i]<-dis[which(dis[,"virus"]==unique(dis[,"virus"])[i] & dis[,"ab"]==unique(dis[,"ab"])[j]), "distance"]
			
		}
		
	}
	out
}




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




