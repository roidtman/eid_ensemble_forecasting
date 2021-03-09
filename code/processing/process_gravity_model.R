#=============================================================================#
# Author: Rachel Oidtman

# Produce department-level gravity model from parameters from Kraemer et al. 2019
#=============================================================================#


#=============================================================================#
# load in epi data
#=============================================================================#

load('../data/processed/processed_epi_data.RData')


#=============================================================================#
# load in shape files 
#=============================================================================#

library(fields)
library(sf)
require(maptools)
require(stringr)

admin_1 <- readShapePoly('../radiation_model/col_admbnda_adm2_unodc_ocha/col_admbnda_adm1_unodc_ocha/col_admbnda_adm1_unodc_ocha.shp')
coordinates <- coordinates(admin_1)
coordinates = as.data.frame(coordinates)
names(coordinates) = c('lon', 'lat')


# update names from shape file to match dept_names
dept_shape_names = as.character(admin_1@data$admin1Name)
dept_shape_names = stringr::str_to_upper(dept_shape_names)
dept_shape_names = stringi::stri_trans_general(dept_shape_names, 'Latin-ASCII')

dept_shape_names[which(dept_shape_names == 'BOGOTA D.C.')] = 'SANTAFE DE BOGOTA D.C'
dept_shape_names[which(dept_shape_names == 'NARINO')] = 'NARIO'

distance_admin_1 = rdist.earth(coordinates[,1:2], coordinates[,1:2], miles = F)

#=============================================================================#
# populate gravity model
# parameter values pulled from Kraemer et al. (2019) 
# https://www.nature.com/articles/s41598-019-41192-3
#=============================================================================#

alpha = 0.78
beta = 0.75
gamma = 2.01


# gravity model function
gravity_model = function(alpha_in = alpha, beta_in = beta, gamma_in = gamma, k_in = 1,
                         pop_origin, pop_destination,
                         distance_in){
  flows_out = k_in * ((pop_origin ^ alpha_in) *  (pop_destination ^ beta_in)) / (distance_in ^ gamma_in)
}


# produce gravity model
# loop over origin (oo) and destination (dd)
connectivity_gravity = matrix(0, length(dept_names), length(dept_names))

for(oo in 1:length(dept_names)){
  for(dd in 1:length(dept_names)){
    
    oo_shape = which(dept_shape_names == dept_names[oo])
    dd_shape = which(dept_shape_names == dept_names[dd])
    
    if(oo != dd){
      connectivity_gravity[oo, dd] = gravity_model(pop_origin = df_stat$pop[oo],
                                                   pop_destination = df_stat$pop[dd],
                                                   distance_in = distance_admin_1[oo_shape, dd_shape])
    }
    
  }
}


# assume diagonal is 0.95
diag_constant = 0.95
connectivity_gravity_normalized = matrix(0, length(dept_names), length(dept_names))
for(rr in 1:nrow(connectivity_gravity)){
  
  connectivity_gravity_normalized[rr, ] = 
    connectivity_gravity[rr,] / (sum(connectivity_gravity[rr,]) / (1-diag_constant))
  
  connectivity_gravity_normalized[rr, rr] = diag_constant
}


save(connectivity_gravity_normalized, file = '../data/processed/processed_gravity_model.RData')

