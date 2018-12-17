################################################################################
# Script for running algos on real data
################################################################################

library(ggplot2)
library(gridExtra)
# --------------------------------- Readin the Data -------------------------------------- #
data1 <- read_in("/Users/gregorbachmann/Desktop/ETH/Semester9/CompBioProject/Data/CAP257/CAP257_bnabdata.csv",9,8,0)

# Read in the CAP256 data file, here some data points get removed due to cleaning
data2 <- read_in("/Users/gregorbachmann/Desktop/ETH/Semester9/CompBioProject/Data/CAP256/CAP256_bnabdata.csv",32,12,1)

# Set data to the set you are interested:
data <- data2

# ----------------------------------- Visualizing --------------------------------------- #
# Isomap
isodata <- data
isodata$symmetric <- log1p(isodata$symmetric)+mean(log1p(isodata$symmetric),na.rm=TRUE)/4
diag(isodata$symmetric) <- 0
isofit <- isomapvis(isodata, 2, 6, plot=FALSE,label=FALSE)


# Interval, Ordinal and Spline
intfit <- list()
ordfit <- list()
splinefit <- list()

intscore <- matrix(0,1000)
ordscore <- matrix(0,1000)
splinescore <- matrix(0,1000)
# Run 100 simulations with random initilization and keep the one with smallest stress
for (run in 1:1000){
  intfit[[run]] <- mdsvis(isodata,2,"interval",init="random",FALSE,FALSE,TRUE)
  intscore[run] <- intfit[[run]]$embedding$stress
  ordfit[[run]] <- mdsvis(isodata,2,"ordinal",init="random",FALSE,FALSE,TRUE)
  ordscore[run] <- ordfit[[run]]$embedding$stress
  splinefit[[run]] <- mdsvis(isodata,2,"mspline",init="random",FALSE,FALSE,TRUE)
  splinescore[run] <- splinefit[[run]]$embedding$stress
}

bestint <- which.min(intscore)
bestord <- which.min(ordscore)
bestspline <- which.min(splinescore)

grid.arrange(isofit$plot,intfit[[bestint]]$plot,ordfit[[bestord]]$plot,splinefit[[bestspline]]$plot,nrow=2)
