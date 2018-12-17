##############################################################
# Multidimensional Scaling for real antigenic data (CAP257) ##
##############################################################


# First read-in the data, automatically cleaning out all viruses with maximal distance to all
# antibodies:

# Read in the CAP257 data file
data1 <- read_in("/Users/gregorbachmann/Desktop/ETH/Semester9/CompBioProject/Data/CAP257/CAP257_bnabdata.csv",9,8,1)

# Read in the CAP256 data file, here some data points get removed due to cleaning
data2 <- read_in("/Users/gregorbachmann/Desktop/ETH/Semester9/CompBioProject/Data/CAP256/CAP256_bnabdata.csv",32,12,1)

# Set data to the set you are interested:
data <- data2

logdata <- data
logdata$standard <- log(1+data$standard)
logdata$symmetric <- log(1+data$symmetric)

#################
#####Isomap######
#################

# I use the package vegan, they implemented a version of isomap that takes distances
# as inputs and not the points from the higher dimensional space, as needed in our case

# 2d version
iso2D <- isomapvis(logdata, 2, 7, plot=TRUE,label=FALSE)

#3d version
iso3D <- isomapvis(logdata, 3, 4, plot=TRUE,label=FALSE)


############################
#Classical Metric Scaling###
############################

# I use the classical metric scaling approach as introduced by Leeuw and Mair.
# I take the "interval" version, meaning the distances will be linearly transformed

# 2d version
interval2d <- mdsvis(logdata, 2, "interval", init="torgerson", plot=TRUE,label=FALSE)

# 3d version
interval3d <- mdsvis(logdata, 3, "interval", init="torgerson", plot=TRUE,label=FALSE)

# Calculate stress level until dimension 10 (after then stress level remains same)
metric_stress <- matrix(0,10)
for (i in 1:10){
  metric_stress[i] <- mdsvis(data, i, "interval", init="torgerson", plot=FALSE,label=FALSE)$stress
}
plot(metric_stress,type="p")

# Dimension 2 seems to be optimal (very weird)


#############################
##Classical Ordinal Scaling##
#############################

# We can use the same weight matrix specified under metric scaling

# 2d version
ordinal2d <- mdsvis(logdata, 2, "ordinal", init="torgerson", plot=TRUE,label=FALSE)
cluster(ordinal2d,4)
# 3d version
ordinal_fit3d <- mdsvis(data, 3, "ordinal", init="torgerson", plot=TRUE,label=FALSE)

# Calculate stress level until dimension 10 (after then stress level remains same)
ordinal_stress <- matrix(0,10)
for (i in 1:15){
  ordinal_stress[i] <- mdsvis(data, i, "ordinal", init="torgerson", plot=FALSE,label=FALSE)$stress
}
plot(ordinal_stress,type="p")


######################
## My Sammon Mapping##
######################

sammon2d <- mysammon(data, dim=2, n.v=20, n.a=12, seed=37, label=FALSE, plot=TRUE)