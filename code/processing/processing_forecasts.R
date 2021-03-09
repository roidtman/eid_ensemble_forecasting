#=============================================================================#
# Author: Rachel Oidtman

# Processing spatial forecast output into RData files...
# HEAVY SCRIPT. TAKES A LONG TIME TO RUN!
#=============================================================================#


#=============================================================================#
# load in processed epi data
#=============================================================================#

load('../data/processed/processed_epi_data.RData')

  
# assign file path
file_path = 'no_tm_ini_2'


#=============================================================================#
# processing functions
#=============================================================================#

write_CI = function(param_posterior, lwr = 0.25, upr = 0.75){
  mat = matrix(0, ncol = ncol(param_posterior), nrow = 3)
  
  for(ii in 1:ncol(param_posterior)){
    mat[1, ii] = quantile(param_posterior[, ii], c(lwr, upr), na.rm = T)[1]
    mat[2, ii] = quantile(param_posterior[, ii], c(lwr, upr), na.rm = T)[2]  
    mat[3, ii] = median(param_posterior[, ii], na.rm = T)
  }
  return(mat)
}

#=============================================================================#
# load in particles
#=============================================================================#

time_btwn_assim = 4
when_to_assimilate = seq(week_first_case + time_btwn_assim, max(df_dyn$week), by = time_btwn_assim)
when_to_assimilate = c(week_first_case, when_to_assimilate)

for(tt in when_to_assimilate){
  
  if(tt == when_to_assimilate[1]){
    f_first_set = paste0('../output/', file_path, '/particles/particles_current_original.csv')
    particles_84 = read.csv(file = f_first_set)
  }else{
    f = paste0('../output/', file_path, '/particles/particles_current_', tt, '.csv')
    
    particles = read.csv(f)
    
    name_tmp = paste0('particles_', tt)
    assign(name_tmp, particles)
    rm(particles)
  }
}


#=============================================================================#
# load in I_F
#=============================================================================#

np = dim(particles_84)[1]

f_out = paste0('../output/', file_path, '/I_F_processed.RData')
if(!file.exists(f_out)){
  I_F_over_time = list()
  length(I_F_over_time) = length(when_to_assimilate)
  
  for(tt in when_to_assimilate){
    
    I_F_over_time[[which(tt == when_to_assimilate)]] = array(NA, dim = c(length(dept_names), np, length(tt:max(df_dyn$week))))
    
    l_time = 1
    for(tt_forecast in tt:max(df_dyn$week)){
      
      f = paste0('../output/', file_path, '/I_F/I_F_', tt, '_', tt_forecast, '.csv')
      I_F_tmp = read.csv(f)
      
      I_F_over_time[[which(tt == when_to_assimilate)]][,,l_time] = as.matrix(I_F_tmp)
      
      l_time = l_time + 1
    }
  }
  ## CrI 
  I_F_CI = list()
  length(I_F_CI) = length(I_F_over_time)
  
  for(tt in when_to_assimilate){
    
    I_F_CI[[which(tt == when_to_assimilate)]] = array(NA, dim = c(length(dept_names), 3, length(tt:max(df_dyn$week))))
    
    for(dd in 1:length(dept_names)){
      I_F_CI[[which(tt == when_to_assimilate)]][dd,,] = write_CI(I_F_over_time[[which(tt == when_to_assimilate)]][dd,,])
    }
  }
  save(I_F_over_time, I_F_CI, file = f_out)
}else{
  load(f_out)
}


#=============================================================================#
# load in I_A
#=============================================================================#

f_out = paste0('../output/', file_path, '/I_A_processed.RData')

# if(!file.exists(f_out)){
  I_A_over_time = array(NA, dim = c(length(dept_names), np, length(week_first_case:max(df_dyn$week))))
  
  for(tt in when_to_assimilate[-1]){
    f = paste0('../output/', file_path, '/I_A/I_A_', tt, '.csv')
    I_A_tmp = read.csv(f)
    
    # I_A_over_time[,,((tt-time_btwn_assim+1)-week_first_case) : (tt - week_first_case)] = as.matrix(I_A_tmp)
    I_A_over_time[,,(tt - week_first_case + 1) : (tt - week_first_case + time_btwn_assim)] = as.matrix(I_A_tmp)
  }
  
  ## CrI
  I_A_CI = array(NA, dim = c(length(dept_names), 3, length(week_first_case:max(df_dyn$week))))
  for(dd in 1:length(dept_names)){
    
    I_A_CI[dd,,] = write_CI(I_A_over_time[dd,,])
  }
  
  save(I_A_CI, I_A_over_time, file = f_out)
  
# }else{
#   load(f_out)
# }


#=============================================================================#
# apply reporting rate to I_F and I_A
#=============================================================================#


f_out = paste0('../output/', file_path, '/I_F_I_A_processed_rho.RData')

# if(!file.exists(f_out)){
  I_A_over_time_rho = I_A_over_time
  I_F_over_time_rho = I_F_over_time
  
  for(tt in when_to_assimilate){
    
    eval(parse(text = paste0(
      'particles = particles_', tt
    )))
    I_A_over_time_rho[,,(tt - week_first_case + 1) : (tt - week_first_case + time_btwn_assim)] = 
      I_A_over_time[,,(tt - week_first_case + 1) : (tt - week_first_case + time_btwn_assim)] * particles$rho
    
    I_F_over_time_rho[[which(tt == when_to_assimilate)]] = I_F_over_time_rho[[which(tt == when_to_assimilate)]] * particles$rho
  }
  
  ## CrI
  I_F_CI_rho = list()
  length(I_F_CI_rho) = length(I_F_over_time)
  I_F_CI_rho_large = list()
  length(I_F_CI_rho_large) = length(I_F_over_time)
  
  for(tt in when_to_assimilate){
    
    I_F_CI_rho[[which(tt == when_to_assimilate)]] = array(NA, dim = c(length(dept_names), 3, length(tt:max(df_dyn$week))))
    I_F_CI_rho_large[[which(tt == when_to_assimilate)]] = array(NA, dim = c(length(dept_names), 3, length(tt:max(df_dyn$week))))
    
    for(dd in 1:length(dept_names)){
      I_F_CI_rho[[which(tt == when_to_assimilate)]][dd,,] = write_CI(I_F_over_time_rho[[which(tt == when_to_assimilate)]][dd,,])
      I_F_CI_rho_large[[which(tt == when_to_assimilate)]][dd,,] = write_CI(I_F_over_time_rho[[which(tt == when_to_assimilate)]][dd,,], lwr = 0.125, upr = 0.875)
    }
  }
  
  I_A_CI_rho = array(NA, dim = c(length(dept_names), 3, length(week_first_case:max(df_dyn$week))))
  I_A_CI_rho_large = array(NA, dim = c(length(dept_names), 3, length(week_first_case:max(df_dyn$week))))
  for(dd in 1:length(dept_names)){
    
    I_A_CI_rho[dd,,] = write_CI(I_A_over_time_rho[dd,,])
    I_A_CI_rho_large[dd,,] = write_CI(I_A_over_time_rho[dd,,], lwr = 0.125, upr = 0.875)
  }
  
  
  # SAMPLE DOWN I_F_OVER_TIME_RHO
  I_F_over_time_rho_sample = list()
  length(I_F_over_time_rho_sample) = length(when_to_assimilate)
  for(tt in when_to_assimilate){
    
    s = sample(1:np, np / 10)
    I_F_over_time_rho_sample[[which(tt == when_to_assimilate)]] = I_F_over_time_rho[[which(tt == when_to_assimilate)]][,s,]
    
  }
  
  
  # SAMPLE DOWN I_A_OVER_TIME_RHO
  I_A_over_time_rho_sample = I_A_over_time_rho[,s,]
  
  
  save(
    # I_A_over_time_rho,
    # I_F_over_time_rho,
    I_F_over_time_rho_sample,
    I_A_over_time_rho_sample,
    I_A_CI_rho,
    I_A_CI_rho_large,
    I_F_CI_rho,
    I_F_CI_rho_large,
    file = f_out
  )
  
# }else{
#   load(f_out)
# }

# }
# 





