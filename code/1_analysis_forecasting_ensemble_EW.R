#=============================================================================#
# Author: Rachel Oidtman

# Producing equally-weighted forecast output
#=============================================================================#



#=============================================================================#
# load in processed epi data
#=============================================================================#

load('../data/processed/processed_epi_data.RData')

# assign file path
runs = 1:16
file_paths = paste0('no_tm_ini_',runs)


#=============================================================================#
# load in processed forecast data
#=============================================================================#

# load ensemble forecast
load('../output/main_text_figures/eq_wts_ensemble_forecast.RData')
I_F_CI_rho_ensemble = I_F_CI_rho ; rm(I_F_CI_rho)
I_F_CI_rho_large_ensemble = I_F_CI_rho_large ; rm(I_F_CI_rho_large)
I_F_over_time_rho_sample_ensemble = I_F_over_time_rho_sample; rm(I_F_over_time_rho_sample)
I_A_CI_rho_ensemble = I_A_CI_rho ; rm(I_A_CI_rho)
I_A_CI_rho_large_ensemble = I_A_CI_rho_large ; rm(I_A_CI_rho_large)
I_A_over_time_rho_sample_ensemble = I_A_over_time_rho_sample ; rm(I_A_over_time_rho_sample)

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


wks_assimilated_at_when_to_assimilate = matrix(NA, nrow = length(when_to_assimilate),
                                               ncol = time_btwn_assim)
for(ii in 1:length(when_to_assimilate)){
  wks_assimilated_at_when_to_assimilate[ii, ] = (when_to_assimilate[ii] - week_first_case + 1) : (when_to_assimilate[ii] - week_first_case + time_btwn_assim)
}

make_figures = F

runs = c(runs, 'ensemble')

#=============================================================================#
# necessary functions
#=============================================================================#

# CrI
write_CI = function(param_posterior, lwr = 0.25, upr = 0.75){
  mat = matrix(0, ncol = ncol(param_posterior), nrow = 3)
  
  for(ii in 1:ncol(param_posterior)){
    mat[1, ii] = quantile(param_posterior[, ii], c(lwr, upr), na.rm = T)[1]
    mat[2, ii] = quantile(param_posterior[, ii], c(lwr, upr), na.rm = T)[2]  
    mat[3, ii] = median(param_posterior[, ii], na.rm = T)
  }
  return(mat)
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


# center data around some center point for plotting
center_data_fxn = function(vec_in, zero_pt_in = 0, 
                           col_less_in = 'yellow', 
                           col_middle_in = 'white',
                           col_greater_in = 'dodgerblue3',
                           col_vec_length = 100, browse = F){
  
  if(browse){browser()}
  
  less_than_zero_prop = length(seq(min(vec_in, na.rm = T), zero_pt_in - 0.001, by = 0.001))
  greater_than_zero_prop = length(seq(zero_pt_in + 0.001, max(vec_in, na.rm = T), by = 0.001))
  
  colfunc_less = colorRampPalette(c(col_less_in, col_middle_in))
  colfunc_greater = colorRampPalette(c(col_middle_in, col_greater_in))
  
  cols = c(colfunc_less(less_than_zero_prop + 1), colfunc_greater(greater_than_zero_prop + 1))
  
  return(cols)
}



#=============================================================================#
# peak week target
#=============================================================================#


# empirical peak week for each department
peak_week_by_dept = sapply(1:length(dept_names), function(dd) 
  which.max(df_dept$cases_new[which(df_dept$dept_name == dept_names[dd])])) - week_first_case
peak_week_by_dept[22] = NA



# peak week forecast CrI
peak_week_forecast = list()
length(peak_week_forecast) = length(runs)
for(rr in 1:length(runs)){
  peak_week_forecast[[rr]] = list()
  length(peak_week_forecast[[rr]]) = length(when_to_assimilate)
}


for(rr in 1:length(runs[-17])){
  
  eval(parse(text = paste0(
    'I_F_over_time_rho_sample = I_F_over_time_rho_sample_', runs[rr]
  )))
  
  eval(parse(text = paste0(
    'I_A_over_time_rho_sample = I_A_over_time_rho_sample_', runs[rr]
  )))
  
  for(tt in 1:length(when_to_assimilate)){
    
    peak_week_forecast[[rr]][[tt]] = matrix(NA, ncol = np, nrow = length(dept_names))
    
    for(dd in 1:length(dept_names)){
      if(dd != 22){
        
        if(file_paths[rr] %in% c('no_tm_ini_5', 'no_tm_ini_6', 'no_tm_ini_7', 'no_tm_ini_8')){
          
          if(tt <= which(when_to_assimilate == when_to_assimilate_NS[[dd]][2]) - 1){
            forecast = round(I_F_over_time_rho_sample[[1]][dd,,])
          }else{
            forecast = round(cbind(I_A_over_time_rho_sample[dd, , 5:(when_to_assimilate[tt] - week_first_case + time_btwn_assim)], 
                             I_F_over_time_rho_sample[[tt]][dd,,]))
          }
        }else{
          if(tt > 1){
            forecast = round(cbind(I_A_over_time_rho_sample[dd, , 5:(when_to_assimilate[tt] - week_first_case + time_btwn_assim)], 
                                   I_F_over_time_rho_sample[[tt]][dd, ,]))
          }else{
            forecast = round(I_F_over_time_rho_sample[[tt]][dd,,])
          }
        }
        peak_week_forecast[[rr]][[tt]][dd, ] = sapply(1:np, function(ff) which.max(forecast[ff, ]))
      }
    }
  }
}


# ensemble forecast
rr = 17
eval(parse(text = paste0(
  'I_F_over_time_rho_sample = I_F_over_time_rho_sample_', runs[rr]
)))

eval(parse(text = paste0(
  'I_A_over_time_rho_sample = I_A_over_time_rho_sample_', runs[rr]
)))

for(tt in 1:length(when_to_assimilate)){
  
  peak_week_forecast[[rr]][[tt]] = matrix(NA, ncol = np, nrow = length(dept_names))
  for(dd in 1:length(dept_names)){
    if(dd != 22){
      if(tt == 1){
        forecast = round(I_F_over_time_rho_sample[[tt]][dd,,])
      }else{
        forecast = matrix(NA, nrow = np, ncol = 60)
        
        for(pp in 1:np){
          
          if(sum(is.na(I_F_over_time_rho_sample[[tt]][dd,pp,])) > 0){ 
            forecast[pp,] = round(I_F_over_time_rho_sample[[1]][dd,pp,])
          }else{
            forecast[pp,] = round(c(I_A_over_time_rho_sample[dd, pp, 5:(when_to_assimilate[tt] - week_first_case + time_btwn_assim)], 
                                        I_F_over_time_rho_sample[[tt]][dd,pp,]))
          }
        }
      } 
      forecast[which(is.na(forecast))] = 0
      peak_week_forecast[[rr]][[tt]][dd, ] = sapply(1:np, function(ff) which.max(forecast[ff, ]))
    }
  }
}


#=============================================================================#
# peak week logscoring
#=============================================================================#

peak_week_by_dept_range = matrix(NA, nrow = length(peak_week_by_dept), ncol = 3)
peak_week_by_dept_range[,1] = peak_week_by_dept - 2
peak_week_by_dept_range[,2] = peak_week_by_dept
peak_week_by_dept_range[,3] = peak_week_by_dept + 2

peak_week_logscore_across_models = array(NA, dim= c(length(runs), 
                                                    length(peak_week_forecast[[1]]), 
                                                    length(dept_names)))
for(rr in 1:length(runs)){
  print(rr)
  for(tt in 1:(length(peak_week_forecast[[1]]))){
    for(dd in 1:length(dept_names)){
      if(dd != 22){
        peak_week_logscore_across_models[rr, tt, dd] = length(which(peak_week_forecast[[rr]][[tt]][dd, ] 
                                                                    %in% peak_week_by_dept_range[dd, 1] : 
                                                                      peak_week_by_dept_range[dd, 3])) / np 
      }
    }
  }
}


if(make_figures){
  
  # pdf('../output/ensemble/peak_week_logscores_ensemble.pdf', height = 6, width = 5)
  cols = colorRampPalette(c('aliceblue', 'navy'))(25)
  par(mar = c(1,1,1,1),
      oma = c(2, 4, 0,1))
  
  layout(matrix(1:32, ncol = 8))
  
  labs = paste0('Model ', runs)
  
  for(dd in 1:length(dept_names)){
    if(dd != 22){
      image(t(peak_week_logscore_across_models[,,dd]), col = cols, 
            yaxt = 'n', xaxt = 'n', 
            zlim = c(0, max(peak_week_logscore_across_models, na.rm = T)))
      axis(side = 1, at = seq(0, 1, length.out = length(when_to_assimilate)),
           labels = NA)
      mtext(side = 3, text = dept_names[dd], cex = 0.5)
      
      if(dd %in% 1:4){
        axis(side = 2, at = seq(0, 1, length.out = length(runs)), las = 1, 
             labels = labs,  cex.axis= 0.7)
      }else{
        axis(side = 2, at = seq(0, 1, length.out = length(runs)), las = 1, 
             labels = NA)
      }
      
    }
  }
  mtext(side = 1, text = 'Weeks of data assimilated into model',line = 3,
        outer = T)
  
  color.bar(cols, min = -1)
  
  dev.off()
}


#=============================================================================#
# peak incidence target
#=============================================================================#

# empirical peak incidence for each department
pk_wk_inds = sapply(1:length(dept_names), function(dd) 
  which.max(df_dept$cases_new[which(df_dept$dept_name == dept_names[dd])]))

peak_incidence_by_dept = sapply(1:length(dept_names), 
                                function(dd) 
                                  df_dept$cases_new[which(df_dept$dept_name == dept_names[dd])][pk_wk_inds[dd]])



# peak week forecast CrI
peak_incidence_forecast = list()
length(peak_incidence_forecast) = length(runs)
for(rr in 1:length(runs)){
  peak_incidence_forecast[[rr]] = list()
  length(peak_incidence_forecast[[rr]]) = length(when_to_assimilate)
}


for(rr in 1:length(runs[-17])){
  
  eval(parse(text = paste0(
    'I_F_over_time_rho_sample = I_F_over_time_rho_sample_', runs[rr]
  )))
  
  eval(parse(text = paste0(
    'I_A_over_time_rho_sample = I_A_over_time_rho_sample_', runs[rr]
  )))
  
  
  for(tt in 1:length(when_to_assimilate)){
    
    peak_incidence_forecast[[rr]][[tt]] = matrix(NA, ncol = np, nrow = length(dept_names))
    
    
    for(dd in 1:length(dept_names)){
      
      if(dd != 22){
        if(file_paths[rr] %in% c('no_tm_ini_5', 'no_tm_ini_6', 'no_tm_ini_7', 'no_tm_ini_8')){
          
          if(tt <= which(when_to_assimilate == when_to_assimilate_NS[[dd]][2]) - 1){
            forecast = round(I_F_over_time_rho_sample[[1]][dd,,])
          }else{
            forecast = round(cbind(I_A_over_time_rho_sample[dd, , 1:(when_to_assimilate[tt] - week_first_case)], 
                             I_F_over_time_rho_sample[[tt]][dd,,]))
          }
        }else{
          if(tt > 1){
            forecast = round(cbind(I_A_over_time_rho_sample[dd, , 5:(when_to_assimilate[tt] - week_first_case + time_btwn_assim)], 
                                   I_F_over_time_rho_sample[[tt]][dd, ,]))
          }else{
            forecast = round(I_F_over_time_rho_sample[[tt]][dd,,])
          }
        }
        
        forecasted_pk_wk_inds = sapply(1:np, function(ff) which.max(forecast[ff, ]))
        
        for(pp in 1:np){
          peak_incidence_forecast[[rr]][[tt]][dd, pp] = forecast[pp, forecasted_pk_wk_inds[pp]]
        }
      }
    }
  }
}

# ensemble forecast
rr = 17
eval(parse(text = paste0(
  'I_F_over_time_rho_sample = I_F_over_time_rho_sample_', runs[rr]
)))

eval(parse(text = paste0(
  'I_A_over_time_rho_sample = I_A_over_time_rho_sample_', runs[rr]
)))

for(tt in 1:length(when_to_assimilate)){
  
  peak_incidence_forecast[[rr]][[tt]] = matrix(NA, ncol = np, nrow = length(dept_names))
  for(dd in 1:length(dept_names)){
    if(dd != 22){
      if(tt == 1){
        forecast = round(I_F_over_time_rho_sample[[tt]][dd,,])
      }else{
        forecast = matrix(NA, nrow = np, ncol = 60)
        
        for(pp in 1:np){
          
          if(sum(is.na(I_F_over_time_rho_sample[[tt]][dd,pp,])) > 0){ 
            forecast[pp,] = round(I_F_over_time_rho_sample[[1]][dd,pp,])
          }else{
            forecast[pp,] = round(c(I_A_over_time_rho_sample[dd, pp, 5:(when_to_assimilate[tt] - week_first_case + time_btwn_assim)], 
                                    I_F_over_time_rho_sample[[tt]][dd,pp,]))
          }
        }
      } 
      forecast[which(is.na(forecast))] = 0
      forecasted_pk_wk_inds = sapply(1:np, function(ff) which.max(forecast[ff, ]))
      
      for(pp in 1:np){
        peak_incidence_forecast[[rr]][[tt]][dd, pp] = forecast[pp, forecasted_pk_wk_inds[pp]]
      }
      
    }
  }
}


#=============================================================================#
# peak incidence logscoring
# to get an approximation of logscores across departments, need to take sum
#=============================================================================#

# empirical peak incidence
empirical_peak_incidence = sapply(1:length(dept_names), function(ff) df_dept$cases[which(df_dept$dept_name == dept_names[ff])][pk_wk_inds[ff]])
peak_incidence_logscore_across_models = array(NA, dim= c(length(runs), 
                                                         length(peak_week_forecast[[1]]), 
                                                         length(dept_names)))


for(rr in 1:length(runs)){
  print(rr)
  for(tt in 1:(length(peak_week_forecast[[1]]))){
    for(dd in 1:length(dept_names)){
      if(dd != 22){
        
        x = peak_incidence_forecast[[rr]][[tt]][dd, ]
        d = density(x)
        fxn = stats::approxfun(d$x, d$y)
        peak_incidence_logscore_across_models[rr, tt, dd] = fxn(empirical_peak_incidence[dd])
        
        # print(peak_incidence_logscore_across_models[rr, tt, dd])

      }
    }
  }
}

# setting all infinite logscores equal to the minimum logscore
peak_incidence_logscore_across_models[which(peak_incidence_logscore_across_models == 0)] = 
  min(peak_incidence_logscore_across_models[-which(peak_incidence_logscore_across_models == 0)], na.rm = T) / 1e4
peak_incidence_logscore_across_models[which(is.na(peak_incidence_logscore_across_models))] = 
  min(peak_incidence_logscore_across_models, na.rm = T) / 1e7

peak_incidence_logscore_across_models = log(peak_incidence_logscore_across_models)


if(make_figures){
  
  # pdf('../output/ensemble/peak_week_logscores_ensemble.pdf', height = 6, width = 5)
  cols = colorRampPalette(c('aliceblue', 'navy'))(25)
  par(mar = c(1,1,1,1),
      oma = c(2, 4, 0,1))
  
  # layout(matrix(c(1,1,1,1,2,
  #                 1,1,1,1,2,
  #                 1,1,1,1,2,
  #                 1,1,1,1,2), 
  #               byrow = T, ncol = 5))
  
  layout(matrix(1:32, ncol = 8))
  
  labs = paste0('Model ', runs)
  
  for(dd in 1:length(dept_names)){
    if(dd != 22){
      image(t(peak_incidence_logscore_across_models[,,dd]), col = cols, 
            yaxt = 'n', xaxt = 'n', 
            zlim = c(min(peak_incidence_logscore_across_models, na.rm = T), 
                     max(peak_incidence_logscore_across_models, na.rm = T)))
      axis(side = 1, at = seq(0, 1, length.out = length(when_to_assimilate)),
           labels = NA)
      mtext(side = 3, text = dept_names[dd], cex = 0.5)
      
      if(dd %in% 1:4){
        axis(side = 2, at = seq(0, 1, length.out = length(runs)), las = 1, 
             labels = labs,  cex.axis= 0.7)
      }else{
        axis(side = 2, at = seq(0, 1, length.out = length(runs)), las = 1, 
             labels = NA)
      }
      
    }
  }
  mtext(side = 1, text = 'Weeks of data assimilated into model',line = 3,
        outer = T)
  
  color.bar(cols, min = -1)
  # axis(side = 4, at = seq(-1, 1, lenght.out = )
  #        seq(0, max(peak_week_logscore_across_models, na.rm = T),
  #                         length.out = length(cols)))
  
  dev.off()
}



#=============================================================================#
# onset week target - week when 10 cumulative cases have occurred
#=============================================================================#


# empirical week by which 10 cases occurred
week_of_first_10_cases = rep(0, length(dept_names))
for(dd in 1:length(dept_names)){
  df_tmp = df_dyn[which(df_dyn$dept_name == dept_names[dd]),]
  week_of_first_10_cases[dd] = df_tmp$week[which(cumsum(df_tmp$cases_new)>=10)][1] - week_first_case
}
week_of_first_10_cases[which(is.na(week_of_first_10_cases))] = 60


# setting up storage
week_of_first_10_cases_forecast = list()
length(week_of_first_10_cases_forecast) = length(runs)
for(rr in 1:length(runs)){
  week_of_first_10_cases_forecast[[rr]] = list()
  length(week_of_first_10_cases_forecast[[rr]]) = length(when_to_assimilate)
}


# week of first case forecast
for(rr in 1:length(runs[-17])){
  
  eval(parse(text = paste0(
    'I_F_over_time_rho_sample = I_F_over_time_rho_sample_', runs[rr]
  )))
  
  eval(parse(text = paste0(
    'I_A_over_time_rho_sample = I_A_over_time_rho_sample_', runs[rr]
  )))
  
  for(tt in 1:length(when_to_assimilate)){
    
    week_of_first_10_cases_forecast[[rr]][[tt]] = matrix(NA, ncol = np, nrow = length(dept_names))
    
    for(dd in 1:length(dept_names)){
      
      if(dd != 22){
        if(file_paths[rr] %in% c('no_tm_ini_5', 'no_tm_ini_6', 'no_tm_ini_7', 'no_tm_ini_8')){
          
          if(tt <= which(when_to_assimilate == when_to_assimilate_NS[[dd]][2]) - 1){
            forecast = round(I_F_over_time_rho_sample[[1]][dd,,])
          }else{
            forecast = round(cbind(I_A_over_time_rho_sample[dd, , 1:(when_to_assimilate[tt] - week_first_case)], 
                             I_F_over_time_rho_sample[[tt]][dd,,]))
          }
        }else{
          if(tt > 1){
            forecast = round(cbind(I_A_over_time_rho_sample[dd, , 5:(when_to_assimilate[tt] - week_first_case + time_btwn_assim)], 
                             I_F_over_time_rho_sample[[tt]][dd, ,]))
          }else{
            forecast = round(I_F_over_time_rho_sample[[tt]][dd,,])
          }
        }
        
        week_of_first_10_cases_forecast[[rr]][[tt]][dd, ] = sapply(1:np, function(ff) 
          which(cumsum(forecast[ff,]) >= 10)[1])
      }
    }
  }
}


# ensemble forecast
rr = 17
eval(parse(text = paste0(
  'I_F_over_time_rho_sample = I_F_over_time_rho_sample_', runs[rr]
)))

eval(parse(text = paste0(
  'I_A_over_time_rho_sample = I_A_over_time_rho_sample_', runs[rr]
)))

for(tt in 1:length(when_to_assimilate)){
  
  week_of_first_10_cases_forecast[[rr]][[tt]] = matrix(NA, ncol = np, nrow = length(dept_names))
  for(dd in 1:length(dept_names)){
    if(dd != 22){
      if(tt == 1){
        forecast = round(I_F_over_time_rho_sample[[tt]][dd,,])
      }else{
        forecast = matrix(NA, nrow = np, ncol = 60)
        
        for(pp in 1:np){
          
          if(sum(is.na(I_F_over_time_rho_sample[[tt]][dd,pp,])) > 0){ 
            forecast[pp,] = round(I_F_over_time_rho_sample[[1]][dd,pp,])
          }else{
            forecast[pp,] = round(c(I_A_over_time_rho_sample[dd, pp, 5:(when_to_assimilate[tt] - week_first_case + time_btwn_assim)], 
                                    I_F_over_time_rho_sample[[tt]][dd,pp,]))
          }
        }
      } 
      forecast[which(is.na(forecast))] = 0
      
      week_of_first_10_cases_forecast[[rr]][[tt]][dd, ] = sapply(1:np, function(ff) 
        which(cumsum(forecast[ff,]) >= 10)[1])
    }
  }
}


#=============================================================================#
# onset week log scoring 
#=============================================================================#

onset_week_by_dept_range = matrix(NA, nrow = length(week_of_first_10_cases), ncol = 3)
onset_week_by_dept_range[,1] = week_of_first_10_cases - 2
onset_week_by_dept_range[,2] = week_of_first_10_cases
onset_week_by_dept_range[,3] = week_of_first_10_cases + 2

onset_week_logscore_across_models = array(NA, dim= c(length(runs), 
                                                    length(week_of_first_10_cases_forecast[[1]]), 
                                                    length(dept_names)))
for(rr in 1:length(runs)){
  print(rr)
  for(tt in 1:length(week_of_first_10_cases_forecast[[1]])){
    for(dd in 1:length(dept_names)){
      if(dd != 22){
        onset_week_logscore_across_models[rr, tt, dd] = length(which(week_of_first_10_cases_forecast[[rr]][[tt]][dd, ] 
                                                                    %in% onset_week_by_dept_range[dd, 1] : 
                                                                      onset_week_by_dept_range[dd, 3])) / np 
      }
    }
  }
}



if(make_figures){
  
  # pdf('../output/ensemble/peak_week_logscores_ensemble.pdf', height = 6, width = 5)
  cols = colorRampPalette(c('aliceblue', 'navy'))(25)
  par(mar = c(1,1,1,1),
      oma = c(2, 4, 0,1))
  
  layout(matrix(1:32, ncol = 8))
  
  labs = paste0('Model ', runs)
  
  for(dd in 1:length(dept_names)){
    if(dd != 22){
      image(t(onset_week_logscore_across_models[,,dd]), col = cols, 
            yaxt = 'n', xaxt = 'n', 
            zlim = c(min(onset_week_logscore_across_models, na.rm = T), 
                     max(onset_week_logscore_across_models, na.rm = T)))
      axis(side = 1, at = seq(0, 1, length.out = length(when_to_assimilate)),
           labels = NA)
      mtext(side = 3, text = dept_names[dd], cex = 0.5)
      
      if(dd %in% 1:4){
        axis(side = 2, at = seq(0, 1, length.out = length(runs)), las = 1, 
             labels = labs,  cex.axis= 0.7)
      }else{
        axis(side = 2, at = seq(0, 1, length.out = length(runs)), las = 1, 
             labels = NA)
      }
      
    }
  }
  mtext(side = 1, text = 'Weeks of data assimilated into model',line = 3,
        outer = T)
  
  color.bar(cols, min = -1)
  
  # dev.off()
}


# save output for easier figure production
save(peak_week_logscore_across_models, 
     onset_week_logscore_across_models,
     peak_incidence_logscore_across_models,
     file = '../output/main_text_figures/output_for_forecast_target_figure.RData')

#=============================================================================#
# Forecast scores averaged over time
#=============================================================================#

dept_names_to_plot = c('La Guajira', 'Magdalena', 'Atlántico', 'Cesar', 'Bolívar',
                       'Sucre', 'Córdoba', 'Norte de Santander', 'Antioquia', 
                       'Chocó', 'Santander', 'Arauca', 'Boyacá', 'Vichada', 
                       'Casanare', 'Cudinamarca', 'Caldas', 'Risaralda', 'Tolima', 
                       'Valle del Cauca', 'Meta', 'Bogota', 'Quindío', 'Guainía', 'Huila', 
                       'Cauca', 'Caquetá', 'Guaviare', 'Nariño', 'Vaupés', 
                       'Putumayo', 'Amazonas')

num_targets = 3

peak_week_logscore_across_models[which(peak_week_logscore_across_models == 0)] = 
  min(peak_week_logscore_across_models[-which(peak_week_logscore_across_models == 0)], na.rm = T) / 1e4
peak_week_logscore_across_models[which(is.na(peak_week_logscore_across_models))] = 
  min(peak_week_logscore_across_models, na.rm = T) / 1e7
logged_peak_week_logscore = log(peak_week_logscore_across_models)

onset_week_logscore_across_models[which(onset_week_logscore_across_models == 0)] = 
  min(onset_week_logscore_across_models[-which(onset_week_logscore_across_models == 0)], na.rm = T) / 1e4
onset_week_logscore_across_models[which(is.na(onset_week_logscore_across_models))] = 
  min(onset_week_logscore_across_models, na.rm = T) / 1e7
logged_onset_week_logscore = log(onset_week_logscore_across_models)


# forecast score by target

forecast_scores_by_target = array(NA, dim = c(num_targets, length(runs), length(dept_names)))

for(tt in 1:3){
  if(tt == 1){
    target = logged_peak_week_logscore
  }else if(tt == 2){
    target = peak_incidence_logscore_across_models
  }else if(tt == 3){
    target = logged_onset_week_logscore
  }
  
  for(dd in 1:length(dept_names)){
    if(dd != 22){
      for(mm in 1:length(runs)){
        forecast_scores_by_target[tt, mm, dd] = sum(target[mm,,dd]) / length(when_to_assimilate)
      }
    }
  }
}

forecast_scores_by_target = exp(forecast_scores_by_target)


# forecast score AVERAGED ACROSS TARGETS
average_LS_model_dept = array(NA, dim = c(length(runs), length(dept_names)))

for(dd in 1:length(dept_names)){
  if(dd != 22){

    for(mm in 1:length(runs)){
      average_LS_model_dept[mm, dd] = sum(
        logged_peak_week_logscore[mm, , dd]
        ,
        peak_incidence_logscore_across_models[mm,,dd]
        ,
        logged_onset_week_logscore[mm,,dd]
      ) /
        (length(when_to_assimilate) * num_targets)
    }
  }
}

forecast_score_model_dept = exp(average_LS_model_dept)


# look at forecast scores relative to baseline ensemble model : target specific
average_FS_relative_to_baseline_targets = array(NA, dim = c(num_targets, 
                                                            nrow(forecast_score_model_dept), 
                                                            ncol(forecast_score_model_dept)))
for(tt in 1:num_targets){
  for(mm in 1:length(runs)){
    average_FS_relative_to_baseline_targets[tt,mm,] = forecast_scores_by_target[tt,mm,] - forecast_scores_by_target[tt,17,]
  }
}
average_FS_relative_to_baseline_targets = average_FS_relative_to_baseline_targets[,,-22]

model_names = c('R-1-CDRs', 'R-2-CDRs', 'Rt-1-CDRs', 'Rt-2-CDRs',
                'R-1-nonspatial', 'R-2-nonspatial', 'Rt-1-nonspatial', 'Rt-2-nonspatial',
                'R-1-gravity', 'R-2-gravity', 'Rt-1-gravity', 'Rt-2-gravity',
                'R-1-radiation', 'R-2-radiation', 'Rt-1-radiation', 'Rt-2-radiation',
                'Ensemble')


pdf('../output/main_text_figures/relative_forecast_score.pdf', height = 8.5, width = 6.5)
par(mar = c(1,6,2,1), oma = c(10,4,1,0))
layout(matrix(1:3, nrow = 3))
layout(matrix(c(
  1,1,1,1,1,4,
  2,2,2,2,2,5,
  3,3,3,3,3,6
),nrow = 3, byrow = T))

# plot x axis ordered by cumulative incidence in each department
cum_cases_by_dept = sapply(1:length(dept_names), function(ff) sum(df_dept$cases_new[which(df_dept$dept_name == dept_names[ff])][week_first_case : max(df_dyn$week)]))
cum_cases_by_dept = cum_cases_by_dept[-22]
xs_to_plot = order(cum_cases_by_dept)

# target names
target_names = c('Timing of peak week', 'Incidence at peak week', 'Week by which 10 cases occurred')

for(tt in 1:num_targets){
  # plot y axis ordered by forecast score
  ys_to_plot = order(rowSums(average_FS_relative_to_baseline_targets[tt,,]))
  
  cols = center_data_fxn(vec_in = as.vector(average_FS_relative_to_baseline_targets[tt,,]), browse = F,
                         col_less_in = 'deeppink4', col_greater_in = 'darkcyan')
  image(t(average_FS_relative_to_baseline_targets[tt, ys_to_plot, rev(xs_to_plot)]), col = cols, xaxt = 'n', yaxt = 'n')
  axis(side = 2, at = seq(0, 1, length.out = length(runs)), labels = model_names[ys_to_plot], las = 1)
  axis(side = 2, at = seq(0, 1, length.out = length(runs))[which(ys_to_plot == 17)], 
       labels = 'Ensemble', las = 1, col = 'red', col.axis = 'red')
  axis(side = 1, at = seq(0, 1, length.out = 31), labels = NA, las = 2) 
  mtext(side = 3, text = target_names[tt], cex = 0.8)
  
  abline(h = seq(0, 1, length.out = length(runs))[which(ys_to_plot == 17)] + (0.0625 / 2), col = 'navy', lwd = 2)
  abline(h = seq(0, 1, length.out = length(runs))[which(ys_to_plot == 17)] - (0.0625 / 2), col = 'navy', lwd = 2)
}
axis(side = 1, at = seq(0, 1, length.out = 31), labels = dept_names_to_plot[-22][rev(xs_to_plot)], las = 2) 
mtext(side = 1, line = 9.5, text = "Departments ordered by overall incidence (high to low)", cex = 0.8)
mtext(side = 2, line = 2, text = 'Models ordered by performance in each metric (low to high)', outer = T, cex = 0.8)


# add colorbars
par(mar = c(1,0,2,5))
for(tt in 1:num_targets){
  cols = center_data_fxn(vec_in = as.vector(average_FS_relative_to_baseline_targets[tt,,]), browse = F,col_less_in = 'deeppink4', col_greater_in = 'darkcyan')
  color.bar(cols, min = -1)
  if(tt == 1){
    mtext(side = 3, text = '% change from', cex = 0.8, line = 1)
    mtext(side = 3, text = 'ensemble', line = 0, cex = 0.8)
  }
  
  y_0 = seq(-1, 1, length.out = length(as.vector(average_FS_relative_to_baseline_targets[tt,,])))[min(which(as.vector(average_FS_relative_to_baseline_targets[tt,,])[order(as.vector(average_FS_relative_to_baseline_targets[tt,,]))] == 0))]
  
  axis(side =4, labels = NA, lwd.ticks = 0)
  axis(side = 4, at = 1, label = round(max(average_FS_relative_to_baseline_targets[tt,,]) * 100), las = 1)
  axis(side = 4, at = -1, label = round(min(average_FS_relative_to_baseline_targets[tt,,]) * 100), las = 1)
  # axis(side = 4, at = y_0, label = 0, las = 1)
  
}

dev.off()



#=============================================================================#
# Forecast scores averaged over departments
#=============================================================================#

num_targets = 3


# forecast score by target
forecast_scores_by_target = array(NA, dim = c(num_targets, length(runs), length(when_to_assimilate)))

for(tt in 1:3){
  if(tt == 1){
    target = logged_peak_week_logscore
  }else if(tt == 2){
    target = peak_incidence_logscore_across_models
  }else if(tt == 3){
    target = logged_onset_week_logscore
  }
  
  for(dd in 1:length(when_to_assimilate)){
    for(mm in 1:length(runs)){
      forecast_scores_by_target[tt, mm, dd] = sum(target[mm,,dd]) / length(dept_names)
    }
  }
}

forecast_scores_by_target = exp(forecast_scores_by_target)


# forecast score AVERAGED ACROSS TARGETS

# look at forecast scores relative to baseline ensemble model : target specific
average_FS_relative_to_baseline_targets = array(NA, dim = c(num_targets, 
                                                            length(model_names), 
                                                            length(when_to_assimilate)))
for(tt in 1:num_targets){
  for(mm in 1:length(runs)){
    average_FS_relative_to_baseline_targets[tt,mm,] = forecast_scores_by_target[tt,mm,] - forecast_scores_by_target[tt,17,]
  }
}



pdf('../output/main_text_figures/relative_forecast_score_time.pdf', height = 8.5, width = 6.5)
par(mar = c(2,6,2,1), oma = c(3,4,1,0))
layout(matrix(1:3, nrow = 3))
layout(matrix(c(
  1,1,1,1,1,4,
  2,2,2,2,2,5,
  3,3,3,3,3,6
),nrow = 3, byrow = T))

# plot x axis ordered by cumulative incidence in each department

# target names
target_names = c('Timing of peak week', 'Incidence at peak week', 'Week by which 10 cases occurred')

for(tt in 1:num_targets){
  # plot y axis ordered by forecast score
  ys_to_plot = order(rowSums(average_FS_relative_to_baseline_targets[tt,,]))
  
  cols = center_data_fxn(vec_in = as.vector(average_FS_relative_to_baseline_targets[tt,,]), browse = F,
                         col_less_in = 'deeppink4', col_greater_in = 'darkcyan')
  image(t(average_FS_relative_to_baseline_targets[tt, ys_to_plot, ]), col = cols, xaxt = 'n', yaxt = 'n')
  axis(side = 2, at = seq(0, 1, length.out = length(runs)), labels = model_names[ys_to_plot], las = 1)
  axis(side = 2, at = seq(0, 1, length.out = length(runs))[which(ys_to_plot == 17)], 
       labels = 'Ensemble', las = 1, col = 'red', col.axis = 'red')
  axis(side = 1, at = seq(0, 1, length.out = length(when_to_assimilate)), 
       labels = when_to_assimilate - week_first_case, las = 2) 
  
  mtext(side = 3, text = target_names[tt], cex = 0.8)
  
  abline(h = seq(0, 1, length.out = length(runs))[which(ys_to_plot == 17)] + (0.0625 / 2), col = 'navy', lwd = 2)
  abline(h = seq(0, 1, length.out = length(runs))[which(ys_to_plot == 17)] - (0.0625 / 2), col = 'navy', lwd = 2)
}
axis(side = 1, at = seq(0, 1, length.out = length(when_to_assimilate)), 
     labels = when_to_assimilate - week_first_case, las = 2) 
mtext(side = 1, line = 11, text = "Departments ordered by overall incidence (high to low)", cex = 0.8)
mtext(side = 2, line = 2, text = 'Models ordered by performance in each metric (low to high)', outer = T, cex = 0.8)


# add colorbars
par(mar = c(1,0,2,5))
for(tt in 1:num_targets){
  cols = center_data_fxn(vec_in = as.vector(average_FS_relative_to_baseline_targets[tt,,]), browse = F,col_less_in = 'deeppink4', col_greater_in = 'darkcyan')
  color.bar(cols, min = -1)
  if(tt == 1){
    mtext(side = 3, text = '% change from', cex = 0.8, line = 1)
    mtext(side = 3, text = 'ensemble', line = 0, cex = 0.8)
  }
  
  y_0 = seq(-1, 1, length.out = length(as.vector(average_FS_relative_to_baseline_targets[tt,,])))[min(which(as.vector(average_FS_relative_to_baseline_targets[tt,,])[order(as.vector(average_FS_relative_to_baseline_targets[tt,,]))] == 0))]
  
  axis(side =4, labels = NA, lwd.ticks = 0)
  axis(side = 4, at = 1, label = round(max(average_FS_relative_to_baseline_targets[tt,,]) * 100), las = 1)
  axis(side = 4, at = -1, label = round(min(average_FS_relative_to_baseline_targets[tt,,]) * 100), las = 1)
  # axis(side = 4, at = y_0, label = 0, las = 1)
  
}

dev.off()

