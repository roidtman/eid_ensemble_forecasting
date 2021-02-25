# install necessary packages
if(!require(fields)){install.packages('fields'); library(fields)}
if(!require(rgeos)){install.packages('rgeos'); library(rgeos)}
if(!require(raster)){install.packages('raster'); library(raster)}

# get command line arguments 
args <- commandArgs(trailingOnly = T)
index <- as.integer(args[1])

# load in shape file 
admin_2 <- shapefile('col_admbnda_adm2_unodc_ocha/col_admbnda_adm2_unodc_ocha.shp')
coordinates <- coordinates(admin_2)
coordinates = as.data.frame(coordinates)
names(coordinates) = c('lon', 'lat')

# calculate the distance to the centroids
distance_admin_2 = rdist.earth(coordinates[,1:2], coordinates[,1:2], miles = F)
diag(distance_admin_2) <- 0

# calculate the minimum distance from the centroid of the specified admin2 to all other admin2 polygons
names_admin_2 <- admin_2$admin2Name
min_distance_admin_2 <- rep(0, length(names_admin_2))

shp.subset <- subset(admin_2, admin2Name == names_admin_2[index])
for(ii in 1:length(min_distance_admin_2))
{
  if(ii != index)
  {
    min_distance_admin_2[ii] <- dist2Line(p = coordinates[ii,], line = shp.subset)[,'distance']/1000
  }
}

# write to file 
write.csv(min_distance_admin_2, file = paste('min_distances/min_distances_', index, '.csv', sep = ''), row.names = F)
