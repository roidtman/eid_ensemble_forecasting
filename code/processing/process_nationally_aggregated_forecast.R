#=============================================================================#
# Author: Rachel Oidtman

# produce nationally aggregated forecasts
#=============================================================================#



#=============================================================================#
# load in processed epi data
#=============================================================================#

require(stats)

load('../data/processed/processed_epi_data.RData')

# assign file path
runs = c(1:16)
file_paths = paste0('no_tm_ini_',runs)


#=============================================================================#
# load in processed forecast data
#=============================================================================#


for(fs in 1:length(file_paths)){
  f = paste0('../output/', file_paths[fs], '/I_F_I_A_processed_rho.RData')
  load(f); rm(f)
  
  if(file_paths[fs] %in% c('no_tm_ini_5', 'no_tm_ini_6', 'no_tm_ini_7', 'no_tm_ini_8')){
    I_A_over_time_rho_sample = I_A_over_time_rho_sample_full
    I_A_CI_rho_large = I_A_CI_rho_large_full
    I_A_CI_rho = I_A_CI_rho_full
    
    I_F_CI_rho  = I_F_CI_rho_full
    I_F_CI_rho_large = I_F_CI_rho_large_full
    I_F_over_time_rho_sample = I_F_over_time_rho_sample_full
    
    rm(I_A_over_time_rho_sample_full, I_A_CI_rho_full, I_A_CI_rho_large_full,
       I_F_CI_rho_full, I_F_CI_rho_large_full, I_F_over_time_rho_sample_full)
    
  }
  assign(paste0('I_A_CI_rho_', runs[fs]), I_A_CI_rho); rm(I_A_CI_rho)
  assign(paste0('I_F_CI_rho_', runs[fs]), I_F_CI_rho); rm(I_F_CI_rho)
  assign(paste0('I_A_over_time_rho_sample_', runs[fs]), I_A_over_time_rho_sample); rm(I_A_over_time_rho_sample)
  assign(paste0('I_F_over_time_rho_sample_', runs[fs]), I_F_over_time_rho_sample); rm(I_F_over_time_rho_sample)
  assign(paste0('I_A_CI_rho_large_', runs[fs]), I_A_CI_rho_large); rm(I_A_CI_rho_large)
  assign(paste0('I_F_CI_rho_large_', runs[fs]), I_F_CI_rho_large); rm(I_F_CI_rho_large)
}


time_btwn_assim = 4
when_to_assimilate = seq(week_first_case + time_btwn_assim, max(df_dyn$week), by = time_btwn_assim)
when_to_assimilate = c(week_first_case, when_to_assimilate)


wks_assimilated_at_when_to_assimilate = matrix(NA, nrow = length(when_to_assimilate),
                                               ncol = time_btwn_assim)
for(ii in 1:length(when_to_assimilate)){
  wks_assimilated_at_when_to_assimilate[ii, ] = (when_to_assimilate[ii] - week_first_case + 1) : 
    (when_to_assimilate[ii] - week_first_case + time_btwn_assim)
}


# non-spatial model when to assimilate
when_to_assimilate_NS = list()
length(when_to_assimilate_NS) = length(dept_names)

time_btwn_assim = 4

for(dd_curr in 1:length(dept_names)){
  if(dd_curr != 22){
    
    df_stat_tmp = df_stat[which(df_stat$dept_name == as.character(dept_names[dd_curr])),]
    df_dyn_tmp = df_dyn[which(df_dyn$dept_name == as.character(dept_names[dd_curr])),]
    week_first_case_tmp = df_stat_tmp$wk_first_case[which(df_stat_tmp$dept_name == as.character(dept_names[dd_curr]))]
    
    when_to_assimilate_tmp = seq(week_first_case_tmp + time_btwn_assim, max(df_dyn$week), by = time_btwn_assim)
    when_to_assimilate_tmp = c(week_first_case_tmp, when_to_assimilate_tmp)
    
    when_to_assimilate_dept_specific = when_to_assimilate[which((when_to_assimilate - week_first_case_tmp) > 0)]
    when_to_assimilate_dept_specific = c(week_first_case_tmp, when_to_assimilate_dept_specific)
    
    when_to_assimilate_NS[[dd_curr]] = when_to_assimilate_dept_specific
  }
}
rm(df_stat_tmp, df_dyn_tmp, week_first_case_tmp)


np = dim(I_A_over_time_rho_sample_1)[2]

dept_models = c(1:21, 23:32)

#=============================================================================#
# necessary functions
#=============================================================================#

# create ensemble weights function

produce_ensemble_wts = function(num_models, LS_in, browse = F){
  
  if(browse) browser()
  
  # one data point per model for national incidence
  # number of departments * number of targets assessed 
  # num_data_points = num_models + length(dept_names)
  num_data_points = ncol(LS_in)
  
  wts = rep(1, num_models) 
  wts = wts / sum(wts)
  wts.new = wts.old = wts
  
  t = 0
  Delta = 1
  epsilon = 1e-12
  
  while(Delta > epsilon){
    t = t + 1
    wts.old = wts.new
    fzpi.old = rep(NA, num_data_points)
    
    for(ii_fxn in 1:num_data_points){
      fzpi.old[ii_fxn] = sum(wts.old * exp(as.numeric(LS_in[,ii_fxn])))
    }
    wts.new = matrix(NA, num_data_points, num_models)
    
    for(ii_fxn in 1:num_data_points){
      wts.new[ii_fxn,] = exp(as.numeric(LS_in[,ii_fxn])) / fzpi.old[ii_fxn]
    }
    wts.new = wts.old * colMeans(wts.new)
    fzpi.new = rep(NA, num_data_points)
    
    for(ii_fxn in 1:num_data_points){
      fzpi.new[ii_fxn] = sum(wts.new * exp(as.numeric(LS_in[,ii_fxn])))
    }
    Delta = (mean(log(fzpi.new)) - mean(log(fzpi.old))) / abs(mean(log(fzpi.new)))
  }
  
  return(wts.new)
}


# incidence logscore
incidence_LS = function(distribution_in, data_in){
  
  d = density(distribution_in)
  fxn = stats::approxfun(d$x, d$y)
  LS = fxn(data_in)
  
  return(log10(LS))
}


# color bar function
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  # dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  # axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}


# credible interval
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
# assess logscores for different models at different points in time
# for example: when first 4 weeks have been assimilated into the model, 
# assess how well the model performed at weeks 1:4
#=============================================================================#

I_A_nationally_aggregated = array(NA, dim = c(length(runs),
                                              np,
                                              length(week_first_case:max(df_dyn$week))))

I_F_nationally_aggregated = list()
length(I_F_nationally_aggregated) = length(when_to_assimilate)
for(ii in 1:length(when_to_assimilate)){
  
  I_F_nationally_aggregated[[ii]] = array(NA, dim = c(length(runs),
                                                      np,
                                                      dim(I_F_over_time_rho_sample_1[[ii]])[3]))
}




for(rr in 1:length(runs)){
  
  eval(parse(text = paste0(
    'I_F_over_time_rho_sample = I_F_over_time_rho_sample_', runs[rr]
  )))
  
  eval(parse(text = paste0(
    'I_A_over_time_rho_sample = I_A_over_time_rho_sample_', runs[rr]
  )))
  
  for(tt in 1:length(when_to_assimilate)){
    # for(tt in 2:length(when_to_assimilate)){
    
    forecast_current = matrix(NA, nrow = np, ncol = dim(I_F_over_time_rho_sample[[tt]])[3])
    
    if(rr %in% 5:8){
      forecast_NS = array(NA, dim = c(length(dept_names), np, dim(I_F_over_time_rho_sample[[tt]])[3]))
      for(dd in 1:length(dept_names)){
        if(dd != 22){
          if(which(when_to_assimilate == when_to_assimilate_NS[[dd]][2]) - 1 <= tt){
            forecast_NS[dd,,] = I_F_over_time_rho_sample[[tt]][dd,,]
          }else{
            forecast_NS[dd,,] = I_F_over_time_rho_sample[[1]][dd,,((when_to_assimilate[tt] - week_first_case) + 1):
                                                                (max(df_dept$week) - week_first_case + 1)]
          }
          
          for(tt_tmp in 1:dim(forecast_NS)[3]){
            forecast_current[,tt_tmp] = sapply(1:np, function(ff) sum(forecast_NS[,ff,tt_tmp], na.rm = T))
          }
        }
      }
    }else{
      for(tt_tmp in 1:dim(I_F_nationally_aggregated[[tt]])[3]){
        forecast_current[,tt_tmp] = sapply(1:np, function(ff) sum(I_F_over_time_rho_sample[[tt]][,ff,tt_tmp], na.rm = T))
      }
    }
    I_F_nationally_aggregated[[tt]][rr,,] = forecast_current
  }
  if(rr %in% 5:8){
    for(tt_tmp in 1:60){
      I_A_nationally_aggregated[rr,,tt_tmp] = sapply(1:np, function(ff) sum(I_A_over_time_rho_sample[,ff,tt_tmp], na.rm = T))
    }
  }else{
    for(tt_tmp in 1:60){
      I_A_nationally_aggregated[rr,,tt_tmp] = sapply(1:np, function(ff) sum(I_A_over_time_rho_sample[,ff,tt_tmp], na.rm = T))
    }
  }
}  


## CrI
I_F_CI_rho_national = list()
length(I_F_CI_rho_national) = length(I_F_nationally_aggregated)
I_F_CI_rho_large_national = list()
length(I_F_CI_rho_large_national) = length(I_F_nationally_aggregated)

for(tt in when_to_assimilate){
  I_F_CI_rho_national[[which(tt == when_to_assimilate)]] = array(NA, dim = c(length(runs), 3, length(tt:max(df_dyn$week))))
  I_F_CI_rho_large_national[[which(tt == when_to_assimilate)]] = array(NA, dim = c(length(runs), 3, length(tt:max(df_dyn$week))))
  for(rr in 1:length(runs)){
    I_F_CI_rho_national[[which(tt == when_to_assimilate)]][rr,,] = write_CI(I_F_nationally_aggregated[[which(tt == when_to_assimilate)]][rr,,])
    I_F_CI_rho_large_national[[which(tt == when_to_assimilate)]][rr,,] = write_CI(I_F_nationally_aggregated[[which(tt == when_to_assimilate)]][rr,,], lwr = 0.125, upr = 0.875)
  }
}


I_A_CI_rho_national = array(NA, dim = c(length(runs), 3, length(week_first_case:max(df_dyn$week))))
I_A_CI_rho_large_national = array(NA, dim = c(length(runs), 3, length(week_first_case:max(df_dyn$week))))
for(rr in 1:length(runs)){
  I_A_CI_rho_national[rr,,] = write_CI(I_A_nationally_aggregated[rr,,])
  I_A_CI_rho_large_national[rr,,] = write_CI(I_A_nationally_aggregated[rr,,], lwr = 0.125, upr = 0.875)  
}




save(I_F_nationally_aggregated, I_A_nationally_aggregated,
     I_A_CI_rho_national, I_A_CI_rho_large_national,
     I_F_CI_rho_national, I_F_CI_rho_large_national,
     file = '../output/nationally_agg_forecasts.RData')

