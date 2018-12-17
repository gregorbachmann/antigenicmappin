########################################################################################
################## generation of test data sets for antigenic mapping ##################
########################################################################################






create.datasets.VA<-function(di,xmax,nvi,nab,fol.n,fil.n){
	
	
	# di:			dimension of the antigenic map
	# xmax:			dim-dimensional vector with maximum values for the axis
	# nvi:			number of viruses
	# nab: 			number of antibodies
	# fol.n:		folder name of where the coordinates/distances/plots-folders are stored, eg fol.n<-"~/Desktop/mapping/testdatasets/"
	# fil.n:		file name of the actual simulation, eg fil.n<-"sim1"
	
	# generate the coordinates of the viruses and antibodies, two didgits after comma
	virloc<-matrix(round(runif(nvi*di,0,xmax),2),byrow=T,nrow=nvi,ncol=di,dimnames=list(paste("V",1:nvi,sep=""),paste("x",1:di,sep="")))
	abloc<-matrix(round(runif(nab*di,0,xmax),2),byrow=T,nrow=nab,ncol=di,dimnames=list(paste("A",1:nab,sep=""),paste("x",1:di,sep="")))
	
	
	# generate output matrix for storing the coordinates
	out<-matrix(NA,nrow=nvi+nab,ncol=2+di,dimnames=list(rep("",nvi+nab),c("character","number",paste("x",1:di,sep=""))))
	out[,"character"]<-c(rep("V",nvi),rep("A",nab))
	out[,"number"]<-c(1:nvi,1:nab)
	out[1:nvi,3:(di+2)]<-virloc
	out[nvi+1:nab,3:(di+2)]<-abloc
	
	write.table(out , paste(fol.n,"coordinates/",fil.n,"_nvi",nvi,"_nab",nab,"_coord.csv",sep="") , quote=F,sep=",",row.names=F)
	
		
	# create a pdf only if 2D-dataset 
	if(di==2){
		hei<-4*xmax[2]/xmax[1]
	
	
		pdf(paste(fol.n,"plots/",fil.n,"_nvi",nvi,"_nab",nab,"_plot.pdf",sep=""),width=4,height=hei)
	
			plot(c(0,xmax[1]),c(0,xmax[2]),xlab="x1 coordinate",ylab="x2 coordinate",xlim=c(0,xmax[1]),ylim=c(0,xmax[2]),type="n")
		
			text(virloc[,1],virloc[,2],labels=paste("V",1:nvi,sep=""),col=1)
			text(abloc[,1],abloc[,2],labels=paste("A",1:nab,sep=""),col=3)
		
		dev.off()
	
	}
	
	# calculate difference between the points
	
	
	dist<-matrix(NA,nrow=0,ncol=3,dimnames=list(c(),c("virus","ab","distance")))
	for(i in 1:nvi){
		for(j in 1:nab){
			dist<-rbind(dist,c(i,j,(sum((virloc[i,]-abloc[j,])^di))^(1/di)))
			
		}
		
		
	}
	write.table(dist , paste(fol.n,"distances/",fil.n,"_nvi",nvi,"_nab",nab,"_dist.csv",sep="") , quote=F,sep=",",row.names=F)
	
}


create.datasets.VA<-function(di,xmax,nvi,nab,fol.n,fil.n){
	
	
	# di:			dimension of the antigenic map
	# xmax:			dim-dimensional vector with maximum values for the axis
	# nvi:			number of viruses
	# nab: 			number of antibodies
	# fol.n:		folder name of where the coordinates/distances/plots-folders are stored, eg fol.n<-"~/Desktop/mapping/testdatasets/virusantibodies"
	# fil.n:		file name of the actual simulation, eg fil.n<-"sim1"
	
	# generate the coordinates of the viruses and antibodies, two didgits after comma
	virloc<-matrix(round(runif(nvi*di,0,xmax),2),byrow=T,nrow=nvi,ncol=di,dimnames=list(paste("V",1:nvi,sep=""),paste("x",1:di,sep="")))
	abloc<-matrix(round(runif(nab*di,0,xmax),2),byrow=T,nrow=nab,ncol=di,dimnames=list(paste("A",1:nab,sep=""),paste("x",1:di,sep="")))
	
	
	# generate output matrix for storing the coordinates
	out<-matrix(NA,nrow=nvi+nab,ncol=2+di,dimnames=list(rep("",nvi+nab),c("character","number",paste("x",1:di,sep=""))))
	out[,"character"]<-c(rep("V",nvi),rep("A",nab))
	out[,"number"]<-c(1:nvi,1:nab)
	out[1:nvi,3:(di+2)]<-virloc
	out[nvi+1:nab,3:(di+2)]<-abloc
	
	write.table(out , paste(fol.n,"coordinates/",fil.n,"_nvi",nvi,"_nab",nab,"_coord.csv",sep="") , quote=F,sep=",",row.names=F)
	
		
	# create a pdf only if 2D-dataset 
	if(di==2){
		hei<-4*xmax[2]/xmax[1]
	
	
		pdf(paste(fol.n,"plots/",fil.n,"_nvi",nvi,"_nab",nab,"_plot.pdf",sep=""),width=4,height=hei)
	
			plot(c(0,xmax[1]),c(0,xmax[2]),xlab="x1 coordinate",ylab="x2 coordinate",xlim=c(0,xmax[1]),ylim=c(0,xmax[2]),type="n")
		
			text(virloc[,1],virloc[,2],labels=paste("V",1:nvi,sep=""),col=1)
			text(abloc[,1],abloc[,2],labels=paste("A",1:nab,sep=""),col=3)
		
		dev.off()
	
	}
	
	# calculate difference between the points
	
	
	dist<-matrix(NA,nrow=0,ncol=3,dimnames=list(c(),c("virus","ab","distance")))
	for(i in 1:nvi){
		for(j in 1:nab){
			dist<-rbind(dist,c(i,j,(sum((virloc[i,]-abloc[j,])^di))^(1/di)))
			
		}
		
		
	}
	write.table(dist , paste(fol.n,"distances/",fil.n,"_nvi",nvi,"_nab",nab,"_dist.csv",sep="") , quote=F,sep=",",row.names=F)
	
}




### function to simulate a data set with additional plasma samples


create.datasets.VAP<-function(di,xmax,nvi,nab,npl,fol.n,fil.n){
	
	
	# di:			dimension of the antigenic map
	# xmax:			dim-dimensional vector with maximum values for the axis
	# nvi:			number of viruses
	# nab: 			number of antibodies
	# npl:			number of plasma points
	# fol.n:		folder name of where the coordinates/distances/plots-folders are stored, eg fol.n<-"~/Desktop/mapping/testdatasets/virusantibodiesplasma/"
	# fil.n:		file name of the actual simulation, eg fil.n<-"sim1"
	
	# generate the coordinates of the viruses and antibodies, two didgits after comma
	virloc<-matrix(round(runif(nvi*di,0,xmax),2),byrow=T,nrow=nvi,ncol=di,dimnames=list(paste("V",1:nvi,sep=""),paste("x",1:di,sep="")))
	abloc<-matrix(round(runif(nab*di,0,xmax),2),byrow=T,nrow=nab,ncol=di,dimnames=list(paste("A",1:nab,sep=""),paste("x",1:di,sep="")))
	
	
	
	# generate the linear combinations of the plasma points:
	# generate random weights between 0 and 1
	# set 20% to 0
	# normalize the plasma weights
	
	alphas<-matrix(round(runif(npl*nab,0,1),5) * rbinom(npl*nab,1,0.8) ,nrow=npl,ncol=nab,dimnames=list(paste("P",1:npl,sep=""), paste("A",1:nab,sep="")))
	alphas<-alphas/rowSums(alphas)
	
	
	### store the alpha matrix
	write.table(alphas , paste(fol.n,"weights/",fil.n,"_nvi",nvi,"_nab",nab,"_weights.csv",sep="") , quote=F,sep=",")
	
	plloc<-matrix(NA,nrow=npl,ncol=di,dimnames=list(paste("P",1:npl,sep=""),paste("x",1:di,sep="")) )
	for(i in 1:npl){
		for(j in 1:di){			
			plloc[i,j]<-sum( alphas[i,]*abloc[,j] )	
		}		
	}
	
	
	# generate output matrix for storing the coordinates
	out<-matrix(NA,nrow=nvi+nab+npl,ncol=2+di,dimnames=list(rep("",nvi+nab+npl),c("character","number",paste("x",1:di,sep=""))))
	out[,"character"]<-c(rep("V",nvi),rep("A",nab),rep("P",npl))
	out[,"number"]<-c(1:nvi,1:nab,1:npl)
	out[1:nvi,3:(di+2)]<-virloc
	out[nvi+1:nab,3:(di+2)]<-abloc
	out[nvi+nab+1:npl,3:(di+2)]<-plloc
	
	write.table(out , paste(fol.n,"coordinates/",fil.n,"_nvi",nvi,"_nab",nab,"_coord.csv",sep="") , quote=F,sep=",",row.names=F)
	
		
	# create a pdf only if 2D-dataset 
	if(di==2){
		hei<-4*xmax[2]/xmax[1]
	
	
		pdf(paste(fol.n,"plots/",fil.n,"_nvi",nvi,"_nab",nab,"_plot.pdf",sep=""),width=4,height=hei)
			#### add the plasma points 
			plot(c( min(c(0,as.numeric(out[,"x1"]))), max(c(xmax[1],as.numeric(out[,"x1"]))) ), c( min(c(0,as.numeric(out[,"x2"]))), max(c(xmax[2],as.numeric(out[,"x2"]))) ), xlab="x1 coordinate",ylab="x2 coordinate",type="n")
		

			text(virloc[,1],virloc[,2],labels=paste("V",1:nvi,sep=""),col=1)
			text(abloc[,1],abloc[,2],labels=paste("A",1:nab,sep=""),col=3)
			text(plloc[,1],plloc[,2],labels=paste("P",1:npl,sep=""),col=4)
			
		dev.off()
	
	}
	
	# calculate difference between the points
	
	### distance also to plasma!
	dist<-matrix(NA,nrow=0,ncol=4,dimnames=list(c(),c("virus","inhibitor","number","distance")))
	for(i in 1:nvi){
		for(j in 1:nab){
			dist<-rbind(dist,c(i,"A",j,(sum((virloc[i,]-abloc[j,])^di))^(1/di)))
			
		}
		for(k in 1:npl){
			dist<-rbind(dist,c(i,"P",k,(sum((virloc[i,]-plloc[k,])^di))^(1/di)))
		}
		
		
	}
	write.table(dist , paste(fol.n,"distances/",fil.n,"_nvi",nvi,"_nab",nab,"_dist.csv",sep="") , quote=F,sep=",",row.names=F)
	
}




