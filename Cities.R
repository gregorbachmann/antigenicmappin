#######################################################################
# Illustration of Sammon/Metric with City Data
#######################################################################

library(MASS)
#-------------------------Load and prepare Data-------------------------------#
# Load data
data(world.cities)
all_coordinates <- cbind(world.cities$lat,world.cities$long)

# Select all capital cities
indices <- which(world.cities$capital==1)
coordinates <- all_coordinates[indices,]

# Calculate distance matrix
dist_matrix <- as.matrix(dist(coordinates))

#-----------------------------------------------------------------------------#


#------------------------------Sammon Mapping---------------------------------#

fit_sammon <- sammon(dist_matrix)

gg_sammon <- ggplot(as.data.frame(fit_sammon$points), aes(x = fit_sammon$points[,1],y= fit_sammon$points[,2],label = world.cities$name[indices]))
gg_sammon + geom_text(check_overlap = TRUE) + geom_point()  

#-----------------------------------------------------------------------------#


#----------------------------Metric Scaling-----------------------------------#

fit_metric <- cmdscale(dist_matrix)
gg_metric <- ggplot(as.data.frame(fit_metric), aes(x = fit_metric[,1],y= fit_metric[,2],label = world.cities$name[indices]))
gg_metric + geom_text(check_overlap = TRUE)  

#-----------------------------------------------------------------------------#