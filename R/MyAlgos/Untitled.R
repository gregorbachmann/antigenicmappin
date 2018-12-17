################################################################################
# Script for evaluating algos
################################################################################

# Evaluate using these dimensions for true shape space
dimensions <- seq(2,20,by=1)

#------------------------------------Isomap------------------------------------#

iso_score <- matrix(0, length(dimensions))
index <- 1
for dim in 1:dimensions{
  iso_score <- test_algos(dim=dim, numtests=500, numpoints=60,method="Isomap")
  index <- index + 1
} 
#------------------------------------------------------------------------------#