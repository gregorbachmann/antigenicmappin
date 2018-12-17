############################################################
################ Protocol antigenic mapping ################
############################################################


### 03/12/2015
### 2D testdatasets

source("/Users/gregorbachmann/Desktop/CompBioProject")
di<-2
xmax<-c(10,10)
nvi<-10
nab<-5
fol.n<-"~/Desktop/mapping/testdatasets/virusantibodies/"


for(i in 1:10){

	fil.n<-paste("sim",i,sep="")
	create.datasets.VA(di,xmax,nvi,nab,fol.n,fil.n)
}


### 04/12/2015





source("~/Desktop/mapping/scripts/testdatasets.R")
di<-2
xmax<-c(10,10)
nvi<-10
nab<-5
npl<-4
fol.n<-"~/Desktop/mapping/testdatasets/virusantibodiesplasma/"


for(i in 1:10){

	fil.n<-paste("sim",i,sep="")
	create.datasets.VAP(di,xmax,nvi,nab,npl,fol.n,fil.n)
}





### 03/03/2016
# Testing mapping on real data 
source("~/work/projects/mapping/scripts/placingalgo.R")

dis<-read.csv("/Users/gregorbachmann/Desktop/CompBioProject/CAP257/CAP257_bnabdata.csv")
 k <- 2 
 niter<-1000
 magic <- 0.2
 tol <- 1e-04
fol.out<-"~/work/projects/mapping/results/CAP257/"
fn.out<-"CAP257testrun"
plot.estcoord.uncertainty(dis,k,niter,magic,tol,fol.out,fn.out)


dis2<-dis
dis2[,"distance"]<-10^dis[,"distance"]
fn.out<-"CAP257testrun2"
plot.estcoord.uncertainty(dis2,k,niter,magic,tol,fol.out,fn.out)


dis3<-dis
dis3[which(dis3[,"distance"]==26),"distance"]<-NA
fn.out<-"CAP257testrun3"
plot.estcoord.uncertainty(dis3,k,niter,magic,tol,fol.out,fn.out)

dis4<-dis
dis4[which(dis4[,"distance"]==26),"distance"]<-24
fn.out<-"CAP257testrun4"
plot.estcoord.uncertainty(dis4,k,niter,magic,tol,fol.out,fn.out)

dis5<-dis
dis5[which(dis5[,"distance"]>16),"distance"]<-17
fn.out<-"CAP257testrun5"
plot.estcoord.uncertainty(dis5,k,niter,magic,tol,fol.out,fn.out)


for(i in unique(dis[,"virus_week"])){

	fn.out<-paste("CAP257testrun6_week",i,sep="")
	plot.estcoord.uncertainty(dis[which(dis[,"virus_week"]==i),],k,niter,magic,tol,fol.out,fn.out)

}










### 08/03/2016

# graphs of reconstructed maps based on simulated data

source("~/work/projects/mapping/scripts/placingalgo.R")

k<-2
niter<-1000
magic<-0.2
tol<-1e-4
fol.out<-"~/work/projects/mapping/results/test_placing_unc/"

dis<-read.csv("~/work/projects/mapping/testdatasets/virusantibodies/distances/sim1_nvi10_nab5_dist.csv")
set.seed(11)
# dis[,"distance"]<-abs(dis[,"distance"]+runif(50,-.5,.5))
fn.out<-"test_error+-0"
plot.estcoord.uncertainty(dis,k,niter,magic,tol,fol.out,fn.out)


# graph of original data
coor<-read.csv("~/work/projects/mapping/testdatasets/virusantibodies/coordinates/sim1_nvi10_nab5_coord.csv")
	xlim<-range(coor[,"x1"])
		xlim[1]<-floor(xlim[1])
		xlim[2]<-ceiling(xlim[2])
		ylim<-range(coor[,"x2"])
		ylim[1]<-floor(ylim[1])
		ylim[2]<-ceiling(ylim[2])

pdf("~/work/projects/mapping/testdatasets/virusantibodies/plots/sim1_nvi10_nab5_niceplot.pdf",width=4,height=4)


plot(c(0,1),c(0,1),xlim=xlim,ylim=ylim,type="n", xlab="antigenic coordinate 1",ylab="antigenic coordinate 2", main="original data")
plot.coord(coor,xlim,ylim,add=TRUE)
		
		
	dev.off()



### using only antibodies that don't have any value ">25"
source("~/work/projects/mapping/scripts/placingalgo.R")
k<-2
niter<-1000
magic<-0.2
tol<-1e-4
fol.out<-"~/work/projects/mapping/results/CAP257/"
fn.out<-"CAP257_MPER+V3"

dis<-read.csv("~/work/projects/mapping/data/CAP257/CAP257_bnabdata.csv")

dis<-dis[which(dis[,"abname"] %in% c("x10E8","x4E10","PGT121","PGT128") ),]

plot.estcoord.uncertainty(dis,k,niter,magic,tol,fol.out,fn.out)


# 09/03/2016
# 10^IC50
fn.out<-"CAP257_MPER+V3_exp"
dis2<-dis
dis2[,"distance"]<-10^dis2[,"distance"]
plot.estcoord.uncertainty(dis2,k,niter,magic,tol,fol.out,fn.out)

### looking at the other antibodiessource("~/work/projects/mapping/scripts/placingalgo.R")
k<-2
niter<-1000
magic<-0.2
tol<-1e-4
fol.out<-"~/work/projects/mapping/results/CAP257/"
fn.out<-"CAP257_V2CD4"
dis<-read.csv("~/work/projects/mapping/data/CAP257/CAP257_bnabdata.csv")
dis<-dis[which(dis[,"abname"] %in% c("PG9","PGT145","VRC01","HJ16") ),]
plot.estcoord.uncertainty(dis,k,niter,magic,tol,fol.out,fn.out)

# 10^IC50
fn.out<-"CAP257_V2CD4_exp"
dis<-read.csv("~/work/projects/mapping/data/CAP257/CAP257_bnabdata.csv")
dis<-dis[which(dis[,"abname"] %in% c("PG9","PGT145","VRC01","HJ16") ),]
dis2<-dis
dis2[,"distance"]<-10^dis2[,"distance"]
plot.estcoord.uncertainty(dis2,k,niter,magic,tol,fol.out,fn.out)




# 10/03/2016
# analysing viral isolates and antibodies isolated from CAP256 (Doria-Rose, Nature 2014, Supplementary Figure 10)
### <0.02 is set to 0.001 (the following lower distance is 0.003)

source("~/work/projects/mapping/scripts/placingalgo.R")
k<-2
niter<-1000
magic<-0.2
tol<-1e-4
fol.out<-"~/work/projects/mapping/results/CAP256/"
fn.out<-"completedata"
dis<-read.csv("~/work/projects/mapping/data/CAP256/CAP256_bnabdata.csv")
plot.estcoord.uncertainty(dis,k,niter,magic,tol,fol.out,fn.out)



# all completely resistant viral strains are excluded

fn.out<-"noresistantstrains"
dis<-read.csv("~/work/projects/mapping/data/CAP256/CAP256_bnabdata.csv")
dis2<-dis[which(!(dis[,"virus"] %in% c(16,20,23:32) )),]
set.seed(101)
plot.estcoord.uncertainty(dis2,k,niter,magic,tol,fol.out,fn.out)

fn.out<-"noresistantstrains_onlypoints"
set.seed(101)
plot.estcoord.uncertainty(dis2,k,niter,magic,tol,fol.out,fn.out,FALSE)


### different seeds
seed<-matrix(NA,nrow=0,ncol=2,dimnames=list(c(),c("seed","stress")))
for(i in 1:1000){
	set.seed(i)
	zw<-plot.estcoord.uncertainty(dis2,k,niter,magic,tol,fol.out,fn.out,FALSE)	
	seed<-rbind(seed,c(i,zw$stress))
}

# minimum stress for seed 407
set.seed(407)
fn.out<-"noresistantstrains_seed407"
plot.estcoord.uncertainty(dis2,k,niter,magic,tol,fol.out,fn.out)
fn.out<-"noresistantstrains_onlypoints_seed407"
set.seed(407)
plot.estcoord.uncertainty(dis2,k,niter,magic,tol,fol.out,fn.out,F)

# maximum stress for seed 190
set.seed(190)
fn.out<-"noresistantstrains_seed190"
plot.estcoord.uncertainty(dis2,k,niter,magic,tol,fol.out,fn.out)
set.seed(190)
fn.out<-"noresistantstrains_onlypoints_seed190"
plot.estcoord.uncertainty(dis2,k,niter,magic,tol,fol.out,fn.out,F)




# all completely resistant viral strains and antibodies from time point 206 are excluded
fn.out<-"noresistantstrains_abweek59+119"
dis<-read.csv("~/work/projects/mapping/data/CAP256/CAP256_bnabdata.csv")
dis3<-dis[which(!(dis[,"virus"] %in% c(16,20,23:32) ) & dis[,"ab_week"]<206),]
set.seed(101)
plot.estcoord.uncertainty(dis3,k,niter,magic,tol,fol.out,fn.out)

fn.out<-"noresistantstrains_abweek59+119_onlypoints"
set.seed(101)
plot.estcoord.uncertainty(dis3,k,niter,magic,tol,fol.out,fn.out,FALSE)

### different seeds
seed<-matrix(NA,nrow=0,ncol=2,dimnames=list(c(),c("seed","stress")))
for(i in 1:1000){
	set.seed(i)
	zw<-plot.estcoord.uncertainty(dis3,k,niter,magic,tol,fol.out,fn.out,FALSE)	
	seed<-rbind(seed,c(i,zw$stress))
}
# minimum
seed[which(seed[,2]==min(seed[,2])),]
# maximum
seed[which(seed[,2]==max(seed[,2])),]

# minimum stress for seed 251
set.seed(seed[which(seed[,2]==min(seed[,2])),1])
fn.out<-paste("noresistantstrains_abweek59+119_seed", seed[which(seed[,2]==min(seed[,2])),1], sep="")
plot.estcoord.uncertainty(dis3,k,niter,magic,tol,fol.out,fn.out)
fn.out<-paste("noresistantstrains_abweek59+119_onlypoints_seed", seed[which(seed[,2]==min(seed[,2])),1], sep="")
set.seed(seed[which(seed[,2]==min(seed[,2])),1])
plot.estcoord.uncertainty(dis3,k,niter,magic,tol,fol.out,fn.out,F)

# maximum stress for seed 850
set.seed(seed[which(seed[,2]==max(seed[,2])),1])
fn.out<-paste("noresistantstrains_abweek59+119_seed", seed[which(seed[,2]==max(seed[,2])),1], sep="")
plot.estcoord.uncertainty(dis3,k,niter,magic,tol,fol.out,fn.out)
set.seed(seed[which(seed[,2]==max(seed[,2])),1])
fn.out<-paste("noresistantstrains_abweek59+119_onlypoints_seed", seed[which(seed[,2]==max(seed[,2])),1], sep="")
plot.estcoord.uncertainty(dis3,k,niter,magic,tol,fol.out,fn.out,F)





### all viruses but only antibodies isolated at week 59 and 119
dis<-read.csv("~/work/projects/mapping/data/CAP256/CAP256_bnabdata.csv")
dis4<-dis[which( dis[,"ab_week"]<206),]
fn.out<-"allviruses_abweek59+119"
### different seeds
seed<-matrix(NA,nrow=0,ncol=2,dimnames=list(c(),c("seed","stress")))
for(i in 1:1000){
	set.seed(i)
	zw<-plot.estcoord.uncertainty(dis4,k,niter,magic,tol,fol.out,fn.out,FALSE)	
	seed<-rbind(seed,c(i,zw$stress))
}
# minimum
seed[which(seed[,2]==min(seed[,2])),]
# maximum
seed[which(seed[,2]==max(seed[,2])),]

# minimum stress 
set.seed(seed[which(seed[,2]==min(seed[,2])),1])
fn.out<-paste("allviruses_abweek59+119_seed", seed[which(seed[,2]==min(seed[,2])),1], sep="")
plot.estcoord.uncertainty(dis4,k,niter,magic,tol,fol.out,fn.out)
fn.out<-paste("allviruses_abweek59+119_onlypoints_seed", seed[which(seed[,2]==min(seed[,2])),1], sep="")
set.seed(seed[which(seed[,2]==min(seed[,2])),1])
plot.estcoord.uncertainty(dis4,k,niter,magic,tol,fol.out,fn.out,F)

# maximum stress 
set.seed(seed[which(seed[,2]==max(seed[,2])),1])
fn.out<-paste("allviruses_abweek59+119_seed", seed[which(seed[,2]==max(seed[,2])),1], sep="")
plot.estcoord.uncertainty(dis4,k,niter,magic,tol,fol.out,fn.out)
set.seed(seed[which(seed[,2]==max(seed[,2])),1])
fn.out<-paste("allviruses_abweek59+119_onlypoints_seed", seed[which(seed[,2]==max(seed[,2])),1], sep="")
plot.estcoord.uncertainty(dis4,k,niter,magic,tol,fol.out,fn.out,F)




