#=============================================================================#
# Author: Rachel Oidtman

# Processing non-spatial forecast output into RData files...
# HEAVY SCRIPT. TAKES AN EXTRA LONG TIME TO RUN!
#=============================================================================#



#=============================================================================#
# load in processed epi data
#=============================================================================#

load('../data/processed/processed_epi_data.RData')

# assign file path
file_path = 'no_tm_ini_8'

# dept models 
dept_models = c(1:21, 23:32)

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
# determine when_to_assimilate across different departments
#=============================================================================#

time_btwn_assim = 4
when_to_assimilate_spatial = seq(week_first_case + time_btwn_assim, max(df_dyn$week), by = time_btwn_assim)
when_to_assimilate_spatial = c(week_first_case, when_to_assimilate_spatial)

when_to_assimilate = list()
length(when_to_assimilate) = length(dept_names)

time_btwn_assim = 4

for(dd_curr in 1:length(dept_names)){
  if(dd_curr != 22){
    
    df_stat_tmp = df_stat[which(df_stat$dept_name == as.character(dept_names[dd_curr])),]
    df_dyn_tmp = df_dyn[which(df_dyn$dept_name == as.character(dept_names[dd_curr])),]
    week_first_case_tmp = df_stat_tmp$wk_first_case[which(df_stat_tmp$dept_name == as.character(dept_names[dd_curr]))]
    
    when_to_assimilate_tmp = seq(week_first_case_tmp + time_btwn_assim, max(df_dyn$week), by = time_btwn_assim)
    when_to_assimilate_tmp = c(week_first_case_tmp, when_to_assimilate_tmp)
    
    when_to_assimilate_dept_specific = when_to_assimilate_spatial[which((when_to_assimilate_spatial - week_first_case_tmp) > 0)]
    when_to_assimilate_dept_specific = c(week_first_case_tmp, when_to_assimilate_dept_specific)
    
    when_to_assimilate[[dd_curr]] = when_to_assimilate_dept_specific
  }
}

#=============================================================================#
# load in particles
#=============================================================================#

for(dd in 1:length(dept_names)){
  for(tt in when_to_assimilate[[dd]]){
    
    if(tt == when_to_assimilate[[dd]][1]){
      f = paste0('../output/', file_path, '/', dd, '/particles/particles_current_original.csv')
    }else{
      f = paste0('../output/', file_path, '/', dd, '/particles/particles_current_', tt, '.csv')
    }
    
    particles = read.csv(f)
    
    name_tmp = paste0('particles_', dd, '_', tt)
    assign(name_tmp, particles)
    rm(particles)
  }
}


#=============================================================================#
# load in I_F
#=============================================================================#


np = dim(particles_1_140)[1]

for(dd in dept_models){
  
  f_out = paste0('../output/', file_path,'/', dd, '/I_F_processed.RData')
  
  if(!file.exists(f_out)){
    I_F_over_time = list()
    length(I_F_over_time) = length(when_to_assimilate_spatial)
    
    for(tt in when_to_assimilate[[dd]]){
      
      if(tt == when_to_assimilate[[dd]][1]){
        I_F_over_time[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]] = 
          array(NA, dim = c(np, length(tt:max(df_dyn$week))))
      }else{
        I_F_over_time[[which(when_to_assimilate_spatial == tt)]] = array(NA, dim = c(np, length(tt:max(df_dyn$week))))
      }
      
      l_time = 1
      for(tt_forecast in tt:max(df_dyn$week)){
        
        f = paste0('../output/', file_path, '/', 
                   dd, '/I_F/I_F_', tt, '_', tt_forecast, '.csv')
        I_F_tmp = read.csv(f)
        
        if(tt == when_to_assimilate[[dd]][1]){
          I_F_over_time[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]][,l_time] = as.matrix(I_F_tmp)
        }else{
          I_F_over_time[[which(when_to_assimilate_spatial == tt)]][,l_time] = as.matrix(I_F_tmp)
        }
        
        l_time = l_time + 1
      }
    }
    
    # load in spin up and cbind it with the first forecast
    f = paste0('../output/', file_path, '/', dd, '/I_F/spin_up.rds')
    spin_up = readRDS(f)
    
    spin_up_mat = array(NA, dim = c(np, 36))
    for(pp in 1:np){
      spin_up_mat[pp, ] = c(rep(NA, 36 - length(spin_up[[pp]][1,])), spin_up[[pp]][1,])
    }
    
    mat_first = cbind(spin_up_mat, I_F_over_time[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]])
    I_F_over_time[[1]] = mat_first[,(ncol(mat_first) - length(week_first_case : max(df_dyn$week)) + 1) : ncol(mat_first)]
      
      
    ## CrI 
    I_F_CI = list()
    length(I_F_CI) = length(I_F_over_time)
    
    l = 1
    for(tt in when_to_assimilate[[dd]]){
      
      if(tt == when_to_assimilate[[dd]][1]){
        I_F_CI[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]] = 
          array(NA, dim = c(3, length(tt:max(df_dyn$week))))
        
        I_F_CI[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]] = 
          write_CI(I_F_over_time[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]])
        
      }else{
        I_F_CI[[which(when_to_assimilate_spatial == tt)]] = array(NA, dim = c(3, length(tt:max(df_dyn$week))))
        
        I_F_CI[[which(when_to_assimilate_spatial == tt)]] = write_CI(I_F_over_time[[which(when_to_assimilate_spatial == tt)]])
      }
    }
    save(I_F_over_time, I_F_CI, file = f_out)
  }else{
    load(f_out)
  }
}




#=============================================================================#
# load in I_A
#=============================================================================#

for(dd in dept_models){
  f_out = paste0('../output/', file_path,'/', dd, '/I_A_processed.RData')
  
  week_first_case_tmp = df_stat$wk_first_case[dd]
  
  # if(!file.exists(f_out)){
    I_A_over_time = array(NA, dim = c(np, length(week_first_case:max(df_dyn$week))))
    
    for(tt in when_to_assimilate[[dd]][-1]){
      
      f = paste0('../output/', file_path, '/', dd, 
                 '/I_A/I_A_', tt, '.csv')
      I_A_tmp = read.csv(f)
      I_A_over_time[,(tt - week_first_case + 1) : (tt - week_first_case + time_btwn_assim)] = as.matrix(I_A_tmp)
    }
    
    ## CrI
    I_A_CI = array(NA, dim = c(3, length(week_first_case:max(df_dyn$week))))
    I_A_CI = write_CI(I_A_over_time)
    
    save(I_A_CI, I_A_over_time, file = f_out)
    
  # }else{
  #   load(f_out)
  # }
  
}


#=============================================================================#
# apply reporting rate to I_F and I_A
#=============================================================================#

for(dd in dept_models){
  
  f_out = paste0('../output/', file_path, '/', dd, '/I_F_I_A_processed_rho.RData')
  
  f_load = paste0('../output/', file_path, '/', dd, '/I_A_processed.RData')
  load(f_load)
  
  f_load = paste0('../output/', file_path, '/', dd, '/I_F_processed.RData')
  load(f_load)
  
  # if(!file.exists(f_out)){
    I_A_over_time_rho = I_A_over_time
    I_F_over_time_rho = I_F_over_time
    
    for(tt in when_to_assimilate[[dd]]){
      
      eval(parse(text = paste0(
        'particles = particles_', dd, '_', tt
      )))
      
      # spin-up rho 
      I_F_over_time_rho[[1]] = I_F_over_time[[1]] * particles$rho
      
      if(tt != when_to_assimilate[[dd]][1]){
        I_A_over_time_rho[,(tt - week_first_case + 1) : (tt - week_first_case + time_btwn_assim)] = 
          I_A_over_time[,(tt - week_first_case + 1) : (tt - week_first_case + time_btwn_assim)] * particles$rho
      }
      
      if(tt == when_to_assimilate[[dd]][1]){
        I_F_over_time_rho[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]] = 
          I_F_over_time[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]] * particles$rho
      }else{
        I_F_over_time_rho[[which(when_to_assimilate_spatial == tt)]] = 
          I_F_over_time[[which(when_to_assimilate_spatial == tt)]] * particles$rho
      }
    }
    
    
    ## CrI
    I_F_CI_rho = list()
    length(I_F_CI_rho) = length(I_F_over_time)
    I_F_CI_rho_large = list()
    length(I_F_CI_rho_large) = length(I_F_over_time)
    
    for(tt in when_to_assimilate[[dd]]){
      
      if(tt == when_to_assimilate[[dd]][1]){
        I_F_CI_rho[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]] = array(NA, dim = c(3, length(tt : max(df_dyn$week))))
        I_F_CI_rho_large[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]] = array(NA, dim = c(3, length(tt : max(df_dyn$week))))
        
        I_F_CI_rho[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]] = 
          write_CI(I_F_over_time_rho[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]])
        I_F_CI_rho_large[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]] = 
          write_CI(I_F_over_time_rho[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]],
                   lwr = 0.125, upr = 0.875)
        
      }else{
        I_F_CI_rho[[which(when_to_assimilate_spatial == tt)]] = array(NA, dim = c(3, length(tt : max(df_dyn$week))))
        I_F_CI_rho_large[[which(when_to_assimilate_spatial == tt)]] = array(NA, dim = c(3, length(tt : max(df_dyn$week))))
        
        I_F_CI_rho[[which(when_to_assimilate_spatial == tt)]] = write_CI(I_F_over_time_rho[[which(when_to_assimilate_spatial == tt)]])
        I_F_CI_rho_large[[which(when_to_assimilate_spatial == tt)]] = write_CI(I_F_over_time_rho[[which(when_to_assimilate_spatial == tt)]],
                                                                               lwr = 0.125, upr = 0.875)
      }
    }
    
    # spin-up CI
    I_F_CI_rho[[1]] = write_CI(I_F_over_time_rho[[1]])
    I_F_CI_rho_large[[1]] = write_CI(I_F_over_time_rho[[1]], lwr  = 0.125, upr = 0.875)
    
    
    I_A_CI_rho = array(NA, dim = c(3, length(week_first_case:max(df_dyn$week))))
    I_A_CI_rho_large = array(NA, dim = c(3, length(week_first_case:max(df_dyn$week))))
    
    I_A_CI_rho = write_CI(I_A_over_time_rho)
    I_A_CI_rho_large = write_CI(I_A_over_time_rho, lwr = 0.125, upr = 0.875)

    
    
    # SAMPLE DOWN I_F_OVER_TIME_RHO
    I_F_over_time_rho_sample = list()
    length(I_F_over_time_rho_sample) = length(when_to_assimilate_spatial)
    
    for(tt in when_to_assimilate[[dd]]){
      
      s = sample(1:np, np / 10)
      
      if(tt == when_to_assimilate[[dd]][1]){
        I_F_over_time_rho_sample[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]] = 
          I_F_over_time_rho[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]][s,]
      }else{
        I_F_over_time_rho_sample[[which(when_to_assimilate_spatial == tt)]] = 
          I_F_over_time_rho[[which(when_to_assimilate_spatial == tt)]][s,]
      }
    }
    
    s = sample(1:np, np / 10)
    I_F_over_time_rho_sample[[1]] = I_F_over_time_rho[[1]][s,]
    
    
    # SAMPLE DOWN I_A_OVER_TIME_RHO
    I_A_over_time_rho_sample = I_A_over_time_rho[s,]
    
    
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
}


#=============================================================================#
# combine output for I_F and I_A across departments
#=============================================================================#


I_A_CI_rho_full = array(NA, dim = c(length(dept_names), 3, length(week_first_case:max(df_dyn$week))))
I_A_CI_rho_large_full = array(NA, dim = c(length(dept_names), 3, length(week_first_case:max(df_dyn$week))))

I_A_over_time_rho_sample_full = array(NA, dim = c(length(dept_names), np / 10, length(week_first_case:max(df_dyn$week))))

I_F_over_time_rho_sample_full = list()
length(I_F_over_time_rho_sample_full) = length(when_to_assimilate_spatial)
for(tt in 1:length(when_to_assimilate_spatial)){
  I_F_over_time_rho_sample_full[[tt]] = array(NA, dim = c(length(dept_names), np / 10, 
                                                          length(when_to_assimilate_spatial[tt] : max(df_dyn$week))))
}

## CrI
I_F_CI_rho_full = list()
length(I_F_CI_rho_full) = length(when_to_assimilate_spatial)

I_F_CI_rho_large_full = list()
length(I_F_CI_rho_large_full) = length(when_to_assimilate_spatial)

for(tt in 1:length(when_to_assimilate_spatial)){
  I_F_CI_rho_full[[tt]] = array(NA, dim = c(length(dept_names), 3, length(when_to_assimilate_spatial[tt] : max(df_dyn$week))))
  I_F_CI_rho_large_full[[tt]] = array(NA, dim = c(length(dept_names), 3, length(when_to_assimilate_spatial[tt] : max(df_dyn$week))))
}




for(dd in 1:length(dept_names)){
  if(dd != 22){
    
    f_load = paste0('../output/', file_path, '/', dd, '/I_F_I_A_processed_rho.RData')
    load(f_load)
    
    I_A_CI_rho_full[dd,,] = I_A_CI_rho
    I_A_CI_rho_large_full[dd,,] = I_A_CI_rho_large
    
    I_A_over_time_rho_sample_full[dd,,] = I_A_over_time_rho_sample
    
    for(tt in when_to_assimilate[[dd]]){
      print(tt)
      
      I_F_over_time_rho_sample_full[[1]][dd,,] = I_F_over_time_rho_sample[[1]]
      
      if(tt == when_to_assimilate[[dd]][1]){
        
        I_F_tmp = cbind(matrix(NA,ncol = (ncol(I_F_over_time_rho_sample_full[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]][dd,,]) - 
                          ncol(I_F_over_time_rho_sample[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]])), nrow =  np / 10),
                        I_F_over_time_rho_sample[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]])
        
        I_F_over_time_rho_sample_full[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]][dd,,] = I_F_tmp
        
        I_F_CI_tmp = cbind(matrix(NA,ncol = (ncol(I_F_over_time_rho_sample_full[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]][dd,,]) - 
                                               ncol(I_F_over_time_rho_sample[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]])), nrow =  3),
                           I_F_CI_rho[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]])
        I_F_CI_rho_full[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]][dd,,] = as.matrix(I_F_CI_tmp)
        
        
        I_F_CI_tmp = cbind(matrix(NA,ncol = (ncol(I_F_over_time_rho_sample_full[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]][dd,,]) - 
                                               ncol(I_F_over_time_rho_sample[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]])), nrow =  3),
                           I_F_CI_rho_large[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]])
        I_F_CI_rho_large_full[[which(when_to_assimilate_spatial == when_to_assimilate[[dd]][2]) - 1]][dd,,] = as.matrix(I_F_CI_tmp)
        
        
      }else{
        I_F_over_time_rho_sample_full[[which(when_to_assimilate_spatial == tt)]][dd,,] = 
          I_F_over_time_rho_sample[[which(when_to_assimilate_spatial == tt)]]
        
        I_F_CI_rho_full[[which(when_to_assimilate_spatial == tt)]][dd,,] = I_F_CI_rho[[which(when_to_assimilate_spatial == tt)]]
        I_F_CI_rho_large_full[[which(when_to_assimilate_spatial == tt)]][dd,,] = I_F_CI_rho_large[[which(when_to_assimilate_spatial == tt)]]
      }
      
      I_F_CI_rho_full[[1]][dd,,] = I_F_CI_rho[[1]]
      I_F_CI_rho_large_full[[1]][dd,,] = I_F_CI_rho_large[[1]]
    }
  }
}


f_out = paste0('../output/', file_path, '/I_F_I_A_processed_rho.RData')

save(I_F_over_time_rho_sample_full, 
     I_A_over_time_rho_sample_full,
     I_A_CI_rho_full,
     I_A_CI_rho_large_full,
     I_F_CI_rho_full,
     I_F_CI_rho_large_full,
     file = f_out)

