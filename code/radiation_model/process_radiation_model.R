# set working directory
# setwd('~/Documents/Research /Radiation Model/')
setwd('~/Desktop/Radiation Model/')

# clear existing workspace
rm(list = ls())

#=============================================================================#
# load in epi data
#=============================================================================#

load('processed_epi_data.RData')

#=============================================================================#
# install necessary packages 
#=============================================================================#

if(!require(fields)){install.packages('fields'); library(fields)}
if(!require(sf)){install.packages('sf'); library(sf)}
if(!require(maptools)){install.packages('maptools'); library(maptools)}
if(!require(rgeos)){install.packages('rgeos'); library(rgeos)}
if(!require(rgdal)){install.packages('rgdal'); library(rgdal)}

#=============================================================================#
# load in shape files 
#=============================================================================#

# load shapefile
# admin_2 <- shapefile('col_admbnda_adm2_unodc_ocha/col_admbnda_adm2_unodc_ocha.shp')
admin_2 <- readShapePoly('~/Desktop/Radiation Model/col_admbnda_adm2_unodc_ocha/col_admbnda_adm2_unodc_ocha.shp')
coordinates <- coordinates(admin_2)
coordinates = as.data.frame(coordinates)
names(coordinates) = c('lon', 'lat')

# calculate distance between all municipalities
distance_admin_2 = rdist.earth(coordinates[,1:2], coordinates[,1:2], miles = F)
diag(distance_admin_2) <- 0

# update naming conventions to match from Siraj et al. 
dept_shape_names = as.character(admin_2@data$admin1Name)
dept_shape_names = stringr::str_to_upper(dept_shape_names)
dept_shape_names = stringi::stri_trans_general(dept_shape_names, 'Latin-ASCII')

dept_shape_names[which(dept_shape_names == 'BOGOTA D.C.')] = 'SANTAFE DE BOGOTA D.C'
dept_shape_names[which(dept_shape_names == 'NARINO')] = 'NARIO'
admin_2@data$admin1Name <- dept_shape_names

muni_shape_names <- admin_2@data$admin2Name

#=============================================================================#
# populate radiation model
#=============================================================================#
# calculate the minimum distance from each admin 2 to all other admin2 polygons 
min_distance_admin_2 <- matrix(0, nrow = nrow(distance_admin_2), ncol = ncol(distance_admin_2))


for(ii in 1:ncol(min_distance_admin_2))
{
  dists <- read.csv(paste('min_distances/min_distances_', ii, '.csv', sep = ''), header = T)[,1]
  min_distance_admin_2[,ii] <- dists
}

# radiation model function 
radiation_model <- function(origin, destination, distance_admin_2, min_distance_admin_2, pop)
{
  # get the radius 
  radius <- distance_admin_2[origin,destination]
  
  # calculate total population with in the radius - excluding origin and destination
  ids <- which(min_distance_admin_2[origin,] <= radius)
  ids <- ids[ids != origin & ids != destination]
  
  pop.circle <- sum(pop[ids])
  
  # calculate the flux from origin to destination
  x <- (pop[origin] + pop.circle)
  y <- (pop[origin] + pop[destination] + pop.circle)
  flux <- (pop[origin] * pop[destination]) / (as.numeric(x) * as.numeric(y))
  
  # return flux
  return(flux)
  
}

# produce radiation model at municipatliy level
# loop over origin (oo) and destination (dd)
muni_connectivity_radiation = matrix(0, length(muni_shape_names), length(muni_shape_names))
pop <- admin_2@data$POPsum

for(oo in 1:length(muni_shape_names)){
  print(oo)
  for(dd in 1:length(muni_shape_names)){
    
    if(oo != dd){
      muni_connectivity_radiation[oo, dd] = radiation_model(origin = oo,
                                                       destination = dd,
                                                       distance_admin_2 = distance_admin_2,
                                                       min_distance_admin_2 = min_distance_admin_2,
                                                       pop = pop)
    }
    
  }
}

# aggregate to department level 
dept_shape_names_unique <- unique(dept_shape_names)
connectivity_radiation <- matrix(0, length(dept_shape_names_unique), length(dept_shape_names_unique))
for(oo in 1:length(dept_shape_names_unique))
{
  for(dd in 1:length(dept_shape_names_unique))
  {
    if(oo != dd)
    {
      indices = as.matrix(expand.grid(which(admin_2@data$admin1Name == dept_shape_names_unique[oo]), 
                            which(admin_2@data$admin1Name == dept_shape_names_unique[dd])))
      connectivity_radiation[oo,dd] <- sum(muni_connectivity_radiation[indices])
    }
  }
}



#=============================================================================#
# re-arrange indices in radiation model to match gravity model
#=============================================================================#

connectivity_radiation_indices = matrix(0, length(dept_names), length(dept_names))

for(oo in 1:length(dept_names)){
  for(dd in 1:length(dept_names)){
    
    oo_shape = which(dept_shape_names_unique == dept_names[oo])
    dd_shape = which(dept_shape_names_unique == dept_names[dd])
    
    if(oo != dd){
      connectivity_radiation_indices[oo, dd] = connectivity_radiation[oo_shape, dd_shape]
    }
    
  }
}



# assume diagonal is 0.95
diag_constant = 0.95
connectivity_radiation_normalized = matrix(0, nrow(connectivity_radiation_indices), ncol(connectivity_radiation_indices))
for(rr in 1:nrow(connectivity_radiation_indices)){
  
  connectivity_radiation_normalized[rr, ] = 
    connectivity_radiation_indices[rr,] / (sum(connectivity_radiation_indices[rr,]) / (1-diag_constant))
  
  connectivity_radiation_normalized[rr, rr] = diag_constant
}




save(connectivity_radiation_normalized, file = 'processed_radiation_model.RData')


