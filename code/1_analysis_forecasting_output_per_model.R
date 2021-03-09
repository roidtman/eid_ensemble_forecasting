#=============================================================================#
# Author: Rachel Oidtman

# Sanity checks and model-specific analyses
# No supplementary figures
#=============================================================================#



#=============================================================================#
# load in processed epi data
#=============================================================================#

library(viridis)
library(lubridate)
library(wesanderson)

load('../data/processed/processed_epi_data.RData')
date_first_case = as_date('2015-08-09')

# assign file path
file_path = 'no_tm_ini_14'
file_out_extension = '_14'

#=============================================================================#
# load in processed forecast data
#=============================================================================#

f = paste0('../output/', file_path, '/I_F_I_A_processed_rho.RData')
load(f); rm(f)
if(file_path %in% c('no_tm_ini_5', 'no_tm_ini_6', 'no_tm_ini_7', 'no_tm_ini_8')){
  I_A_over_time_rho_sample = I_A_over_time_rho_sample_full
  I_A_CI_rho_large = I_A_CI_rho_large_full
  I_A_CI_rho = I_A_CI_rho_full
  
  I_F_CI_rho  = I_F_CI_rho_full
  I_F_CI_rho_large = I_F_CI_rho_large_full
  I_F_over_time_rho_sample = I_F_over_time_rho_sample_full
  
  rm(I_A_over_time_rho_sample_full, I_A_CI_rho_full, I_A_CI_rho_large_full,
     I_F_CI_rho_full, I_F_CI_rho_large_full, I_F_over_time_rho_sample_full)
}


time_btwn_assim = 4
when_to_assimilate = seq(week_first_case + time_btwn_assim, max(df_dyn$week), by = time_btwn_assim)
when_to_assimilate = c(week_first_case, when_to_assimilate)


if(file_path %in% c('no_tm_ini_5', 'no_tm_ini_6', 'no_tm_ini_7', 'no_tm_ini_8')){
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
}



np = dim(I_A_over_time_rho_sample)[2]


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

# color bar horiz function
color.bar.horiz <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  # dev.new(width=1.75, height=5)
  plot(c(min,max), c(0,10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  # axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(y,0,y+1/scale,10, col=lut[i], border=NA)
  }
}


# center data around some center point for plotting
center_data_fxn = function(vec_in, zero_pt_in = median(range(vec_in, na.rm = T)), 
                           col_less_in = 'yellow', 
                           col_middle_in = 'white',
                           col_greater_in = 'dodgerblue3',
                           col_vec_length = 100, browse = F){
  
  if(browse){browser()}
  
  less_than_zero_prop = length(seq(min(vec_in, na.rm = T), zero_pt_in - 0.01, by = 0.01))
  greater_than_zero_prop = length(seq(zero_pt_in + 0.01, max(vec_in, na.rm = T), by = 0.01))
  
  colfunc_less = colorRampPalette(c(col_less_in, col_middle_in))
  colfunc_greater = colorRampPalette(c(col_middle_in, col_greater_in))
  
  cols = c(colfunc_less(less_than_zero_prop + 1), colfunc_greater(greater_than_zero_prop + 1))
  
  return(cols)
}


# get sunday of each week
# for labeling
yearSunday <- function(year) {
  dates <- as.Date(paste(year, "-01-01", sep = "")) + 0:6
  days <- weekdays(dates) == "Sunday"
  seq(dates[days], as.Date(paste(year, "-12-31", sep = "")),  by = "+7 day")
}


#=============================================================================#
# department-specific forecasts
#=============================================================================#


f = paste0('../output/', file_path, '/dept-specific_forecasts', file_out_extension, '.pdf')
pdf(file = f, height = 8, width = 6)

if(file_path %in% c('no_tm_ini_5', 'no_tm_ini_6', 'no_tm_ini_7', 'no_tm_ini_8')){
  
  
  for(dd in 1:length(dept_names)){
    if(dd != 22){
      # dd = 8
      
      I_A_to_plot = log10(round(I_A_CI_rho[dd,,]))
      I_A_to_plot[which(is.infinite(I_A_to_plot))] = 0
      I_A_to_plot[which(is.na(I_A_to_plot))] = 0
      
      I_A_to_plot_large = log10(round(I_A_CI_rho_large[dd,,]))
      I_A_to_plot_large[which(is.infinite(I_A_to_plot_large))] = 0
      I_A_to_plot_large[which(is.na(I_A_to_plot_large))] = 0
      
      cases_now = log10(df_dept$cases_new[which(df_dept$dept_name == dept_names[dd] & df_dept$week %in% week_first_case: max(df_dyn$week))])
      
      layout(matrix(1:15, nrow = 15))
      par(mar = c(0.5,2,0,1), oma = c(3,3,2,1))
      
      for(tt in 1:length(when_to_assimilate)){
        
        
        plot(-100, -100, xlim = c(0, dim(I_A_to_plot)[2]), 
             ylim = c(0, max(I_F_to_plot_large, I_A_to_plot_large, cases_now, na.rm = T)), 
             xlab = '', ylab = '', xaxt = 'n', las = 1)
        
        if(tt <= which(when_to_assimilate == when_to_assimilate_NS[[dd]][2]) - 1){
          
          I_F_to_plot = log10(round(I_F_CI_rho[[1]][dd,,]))
          I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
          
          I_F_to_plot_large = log10(round(I_F_CI_rho_large[[1]][dd,,]))
          I_F_to_plot_large[which(is.infinite(I_F_to_plot_large))] = 0
          
          # 75% CI I_F
          polygon(x = c((1:ncol(I_F_to_plot_large)),
                        rev(1:ncol(I_F_to_plot_large))),
                  y = c(I_F_to_plot_large[1,],
                        rev(I_F_to_plot_large[2,])), col = adjustcolor('slateblue', alpha.f = 0.1), 
                  border = adjustcolor('slateblue', alpha.f = 0.1))
          
          # 50% CI I_F
          polygon(x = c(1:ncol(I_F_to_plot),
                        rev(1:ncol(I_F_to_plot))),
                  y = c(I_F_to_plot[1,],
                        rev(I_F_to_plot[2,])), col = adjustcolor('slateblue', alpha.f = 0.3), 
                  border = adjustcolor('slateblue', alpha.f = 0.3))
          
          
          lines(x = (1:ncol(I_F_to_plot_large)), 
                y = I_F_to_plot[3,], col = 'slateblue', lwd = 4)
          
        }else{
          
          I_F_to_plot = log10(round(I_F_CI_rho[[tt]][dd,,]))
          I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
          
          I_F_to_plot_large = log10(round(I_F_CI_rho_large[[tt]][dd,,]))
          I_F_to_plot_large[which(is.infinite(I_F_to_plot_large))] = 0
          
          
          # 75% CI I_A
          polygon(x = c(1:(when_to_assimilate[tt] - week_first_case),
                        rev(1:(when_to_assimilate[tt] - week_first_case))),
                  y = c(I_A_to_plot_large[1, 1:(when_to_assimilate[tt] - week_first_case)],
                        rev(I_A_to_plot_large[2, 1:(when_to_assimilate[tt] - week_first_case)])),
                  col = adjustcolor('darkolivegreen1', alpha.f = 0.3), 
                  border = adjustcolor('darkolivegreen1', alpha.f = 0.3))
          
          # 50% CI I_A
          polygon(x = c(1:(when_to_assimilate[tt] - week_first_case),
                        rev(1:(when_to_assimilate[tt] - week_first_case))),
                  y = c(I_A_to_plot[1, 1:(when_to_assimilate[tt] - week_first_case)],
                        rev(I_A_to_plot[2, 1:(when_to_assimilate[tt] - week_first_case)])),
                  col = adjustcolor('darkolivegreen1', alpha.f = 0.6), 
                  border = adjustcolor('darkolivegreen1', alpha.f = 0.6))
          
          # 75% CI I_F
          polygon(x = c((when_to_assimilate[tt] - week_first_case) : 
                          (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)),
                        rev((when_to_assimilate[tt] - week_first_case) : 
                              (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)))),
                  y = c(I_F_to_plot_large[1,],
                        rev(I_F_to_plot_large[2,])), 
                  col = adjustcolor('slateblue', alpha.f = 0.1), 
                  border = adjustcolor('slateblue', alpha.f = 0.1))
          
          # 50% CI I_F
          polygon(x = c((when_to_assimilate[tt] - week_first_case) : 
                          (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)),
                        rev((when_to_assimilate[tt] - week_first_case) : 
                              (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)))),
                  y = c(I_F_to_plot[1,],
                        rev(I_F_to_plot[2,])), 
                  col = adjustcolor('slateblue', alpha.f = 0.3), 
                  border = adjustcolor('slateblue', alpha.f = 0.3))
          
          lines(I_A_to_plot[3, 1:(when_to_assimilate[tt] - week_first_case)], col = 'darkolivegreen3', lwd = 4)
          lines(x = (when_to_assimilate[tt] - week_first_case) : 
                  (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)), 
                y = I_F_to_plot[3,], col = 'slateblue', lwd = 4)
          
        }
        
        points(cases_now, cex = 1, pch = 4, col = rgb(0,0,0,0.8))
        # points(cases_now[1:(when_to_assimilate[tt] - week_first_case)], cex = 0.8, pch = 4)
        # points(c(rep(NA, length(1:(when_to_assimilate[tt] - week_first_case))),
        #          cases_now[(when_to_assimilate[tt] - week_first_case):length(cases_now)]),
        #        cex = 0.8, pch = 4, col = rgb(0,0,0,0.3))
        
      }
      
      axis(side = 1)
      
      mtext(side = 1, line = 2.5, text = 'Time')
      mtext(side = 2, line = 1, text = 'log10(Incidence)', outer = T)
      mtext(side = 3, outer = T, text = dept_names[dd])
      
    }
  }
  
  dev.off()
}else{
  
  
  for(dd in 1:length(dept_names)){
    
    # dd = 8
    
    I_A_to_plot = log10(round(I_A_CI_rho[dd,,]))
    I_A_to_plot[which(is.infinite(I_A_to_plot))] = 0
    I_A_to_plot[which(is.na(I_A_to_plot))] = 0
    
    I_A_to_plot_large = log10(round(I_A_CI_rho_large[dd,,]))
    I_A_to_plot_large[which(is.infinite(I_A_to_plot_large))] = 0
    I_A_to_plot_large[which(is.na(I_A_to_plot_large))] = 0
    
    cases_now = log10(df_dept$cases_new[which(df_dept$dept_name == dept_names[dd] & df_dept$week %in% week_first_case: max(df_dyn$week))])
    
    layout(matrix(1:15, nrow = 15))
    par(mar = c(0.5,2,0,1), oma = c(3,3,2,1))
    
    for(tt in 1 : length(when_to_assimilate)){
      
      I_F_to_plot = log10(round(I_F_CI_rho[[tt]][dd,,]))
      I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
      
      I_F_to_plot_large = log10(round(I_F_CI_rho_large[[tt]][dd,,]))
      I_F_to_plot_large[which(is.infinite(I_F_to_plot_large))] = 0
      
      plot(-100, -100, xlim = c(0, dim(I_A_to_plot)[2]), 
           ylim = c(0, max(I_F_to_plot_large, I_A_to_plot_large, cases_now, na.rm = T)), 
           xlab = '', ylab = '', xaxt = 'n', las = 1)
      
      if(tt > 1){
        # 75% CI I_A
        polygon(x = c(1:(when_to_assimilate[tt] - week_first_case),
                      rev(1:(when_to_assimilate[tt] - week_first_case))),
                y = c(I_A_to_plot_large[1, 1:(when_to_assimilate[tt] - week_first_case)],
                      rev(I_A_to_plot_large[2, 1:(when_to_assimilate[tt] - week_first_case)])),
                col = adjustcolor('darkolivegreen1', alpha.f = 0.3), 
                border = adjustcolor('darkolivegreen1', alpha.f = 0.3))
        
        # 50% CI I_A
        polygon(x = c(1:(when_to_assimilate[tt] - week_first_case),
                      rev(1:(when_to_assimilate[tt] - week_first_case))),
                y = c(I_A_to_plot[1, 1:(when_to_assimilate[tt] - week_first_case)],
                      rev(I_A_to_plot[2, 1:(when_to_assimilate[tt] - week_first_case)])),
                col = adjustcolor('darkolivegreen1', alpha.f = 0.6), 
                border = adjustcolor('darkolivegreen1', alpha.f = 0.6))
      }
      
      # 75% CI I_F
      polygon(x = c((when_to_assimilate[tt] - week_first_case) : 
                      (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)),
                    rev((when_to_assimilate[tt] - week_first_case) : 
                          (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)))),
              y = c(I_F_to_plot_large[1,],
                    rev(I_F_to_plot_large[2,])), 
              col = adjustcolor('slateblue', alpha.f = 0.1), 
              border = adjustcolor('slateblue', alpha.f = 0.1))
      
      # 50% CI I_F
      polygon(x = c((when_to_assimilate[tt] - week_first_case) : 
                      (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)),
                    rev((when_to_assimilate[tt] - week_first_case) : 
                          (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)))),
              y = c(I_F_to_plot[1,],
                    rev(I_F_to_plot[2,])), 
              col = adjustcolor('slateblue', alpha.f = 0.3), 
              border = adjustcolor('slateblue', alpha.f = 0.3))
      
      lines(I_A_to_plot[3, 1:(when_to_assimilate[tt] - week_first_case)], col = 'darkolivegreen3', lwd = 4)
      lines(x = (when_to_assimilate[tt] - week_first_case) : 
              (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)), 
            y = I_F_to_plot[3,], col = 'slateblue', lwd = 4)
      
      points(cases_now, cex = 1, pch = 4, col = rgb(0,0,0,0.8))
    }
    
    axis(side = 1)
    
    mtext(side = 1, line = 2.5, text = 'Time')
    mtext(side = 2, line = 1, text = 'log10(Incidence)', outer = T)
    mtext(side = 3, outer = T, text = dept_names[dd])
    
  }
  dev.off()
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
length(peak_week_forecast) = length(when_to_assimilate)

for(tt in 1:length(when_to_assimilate)){
  
  peak_week_forecast[[tt]] = matrix(NA, ncol = np, nrow = length(dept_names))
  
  for(dd in 1:length(dept_names)){
    if(dd != 22){
      
      if(file_path %in% c('no_tm_ini_5', 'no_tm_ini_6', 'no_tm_ini_7', 'no_tm_ini_8')){
        
        if(tt <= which(when_to_assimilate == when_to_assimilate_NS[[dd]][2]) - 1){
          forecast = I_F_over_time_rho_sample[[1]][dd,,]
        }else{
          forecast = cbind(I_A_over_time_rho_sample[dd, , 1:(when_to_assimilate[tt] - week_first_case)], 
                           I_F_over_time_rho_sample[[tt]][dd,,])
        }
      }else{
        if(tt > 1){
          forecast = cbind(I_A_over_time_rho_sample[dd, , 1:(when_to_assimilate[tt] - week_first_case)], 
                           I_F_over_time_rho_sample[[tt]][dd, ,])
        }else{
          forecast = I_F_over_time_rho_sample[[tt]][dd,,]
        }
      }
      peak_week_forecast[[tt]][dd, ] = sapply(1:np, function(ff) which.max(forecast[ff, ]))
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

peak_week_logscore = matrix(NA, nrow = length(peak_week_forecast), ncol = length(dept_names))
for(tt in 1:(length(peak_week_forecast))){
  
  for(dd in 1:length(dept_names)){
    if(dd != 22){
      peak_week_logscore[tt, dd] = length(which(peak_week_forecast[[tt]][dd, ] %in% 
                                                  peak_week_by_dept_range[dd, 1] : 
                                                  peak_week_by_dept_range[dd, 3])) / np 
    }
  }
}


f = paste0('../output/', file_path, '/peak_week_logscore_by_dept', file_out_extension, '.pdf')
pdf(file = f, height = 6, width = 5.5)
cols = colorRampPalette(c('aliceblue', 'navy'))(25)
par(mar = c(2,2,2,2),
    oma = c(4, 10, 0,1))

layout(matrix(c(1,1,1,1,2,
                1,1,1,1,2,
                1,1,1,1,2,
                1,1,1,1,2), 
              byrow = T, ncol = 5))

image((peak_week_logscore), col= cols,
      xaxt = 'n', 
      yaxt = 'n' )
axis(side = 1, at = seq(-0.1, 1.1, length.out = length(when_to_assimilate) + 1),
     labels = c(when_to_assimilate, 143) - week_first_case, las = 2)
mtext(side = 1, text = 'Weeks of data assimilated into model',line = 3)
mtext(side = 3, text = 'Peak week (+/- 1 week)')

labs = dept_names
axis(side = 2, at = seq(0, 1, length.out = length(dept_names)), las = 1, 
     labels = labs)

# add a box for each peak week 
for(dd in 1:length(dept_names)){
  if(dd != 22){
    ind_tmp = which.min(abs(when_to_assimilate - week_first_case - peak_week_by_dept_range[dd,2]))
    points(y = seq(0, 1, length.out = length(dept_names))[dd],
           x = seq(0, 1, length.out = length(when_to_assimilate))[ind_tmp], 
           col = 'deeppink2', pch = 19)
  }
}

color.bar(colorRampPalette(c('aliceblue', 'navy'))(25), min = -1)
mtext(side = 1, text = 'Less accurate')
mtext(side = 3, text = 'More accurate')
mtext(side = 2, text = 'Logscore')
dev.off()



#=============================================================================#
# probability that peak is coming 
# try to do level plots and color region
#=============================================================================#

# empirical peak week 
peak_week_by_dept = sapply(1:length(dept_names), function(dd) 
  which.max(df_dept$cases_new[which(df_dept$dept_name == dept_names[dd])])) - week_first_case
peak_week_by_dept[22] = NA


peak_week_by_dept_range = matrix(NA, nrow = length(peak_week_by_dept), ncol = 3)
peak_week_by_dept_range[,1] = peak_week_by_dept - 2
peak_week_by_dept_range[,2] = peak_week_by_dept
peak_week_by_dept_range[,3] = peak_week_by_dept + 2


# peak week forecast CrI
prob_peak_week_in_future = matrix(NA, nrow = length(dept_names), ncol = length(when_to_assimilate))
modal_peak_week_in_future = matrix(NA, nrow = length(dept_names), ncol = length(when_to_assimilate))


for(tt in 1:length(when_to_assimilate)){
  
  for(dd in 1:length(dept_names)){
    if(dd != 22){
      if(file_path %in% c('no_tm_ini_5', 'no_tm_ini_6', 'no_tm_ini_7', 'no_tm_ini_8')){
        
        if(tt <= which(when_to_assimilate == when_to_assimilate_NS[[dd]][2]) - 1){
          forecast = I_F_over_time_rho_sample[[1]][dd,,]
        }else{
          forecast = cbind(I_A_over_time_rho_sample[dd, , 1:(when_to_assimilate[tt] - week_first_case)], 
                           I_F_over_time_rho_sample[[tt]][dd,,])
        }
      }else{
        if(tt > 1){
          forecast = cbind(I_A_over_time_rho_sample[dd, , 1:(when_to_assimilate[tt] - week_first_case)], 
                           I_F_over_time_rho_sample[[tt]][dd, ,])
        }else{
          forecast = I_F_over_time_rho_sample[[tt]][dd,,]
        }
      }
      
      # prob_peak_week_in_future[dd, l] = length(which(sapply(1:np, function(ff) which.max(forecast[ff, ])) > 
      #                                                      when_to_assimilate[l] - week_first_case)) / np 
      
      modal_peak_week_in_future[dd, tt] = which.max(prop.table(table(sapply(1:np, function(ff) which.max(forecast[ff, ])))))
      prob_peak_week_in_future[dd, tt] = prop.table(table(sapply(1:np, function(ff) which.max(forecast[ff, ]))))[
        which.max(prop.table(table(sapply(1:np, function(ff) which.max(forecast[ff, ])))))
        ]
    }
  }
}

modal_peak_week_targets = matrix(NA, nrow = length(when_to_assimilate),
                                 ncol = time_btwn_assim)
for(ii in 1:length(when_to_assimilate)){
  modal_peak_week_targets[ii, ] = (when_to_assimilate[ii] - week_first_case + 1) : (when_to_assimilate[ii]  + time_btwn_assim - week_first_case)
}

modal_peak_week_range_in_future = modal_peak_week_in_future
for(rr in 1:nrow(modal_peak_week_range_in_future)){
  if(rr != 22){
    for(cc in 1:ncol(modal_peak_week_range_in_future)){
      logs_tmp = sapply(1:nrow(modal_peak_week_targets), 
                        function(ff) modal_peak_week_in_future[rr,cc] %in% modal_peak_week_targets[ff,])
      modal_peak_week_range_in_future[rr,cc] = which(logs_tmp)
    }
  }
}


# f = paste0('../output/', file_path, '/prob_peak_week_in_future_with_forecasts_', file_out_extension, '.pdf')
f = paste0('../output/', file_path, '/modal_week_in_future_with_forecasts_', file_out_extension, '.pdf')
pdf(file = f, height = 6, width = 8)

# color_func = functin(alpha_in){
#   colorRampPalette(c('aliceblue', 'khaki', 'hotpink4'))(32)
# }

# cols = c(rep('gray97', 10), colorRampPalette(c('darkolivegreen3', 
#                                                'lightslateblue', 
#                                                'khaki', 
#                                                'deeppink3'))(35-10))

cols = colorRampPalette(c('gray97',
                          'darkolivegreen3', 
                          'lightslateblue',
                          'khaki', 
                          'deeppink3'))(15)
cols_to_plot_in_color_bar = c(colorRampPalette(c('gray97',
                                                 'darkolivegreen3', 
                                                 'lightslateblue', 
                                                 'khaki', 
                                                 'deeppink3'))(15), 'black')

# image_to_plot = t(modal_peak_week_in_future[order(peak_week_by_dept),])
image_to_plot = t(modal_peak_week_range_in_future[order(peak_week_by_dept),])
# image_to_plot = t(modal_peak_week_range_in_future)

par(mar = c(3.5,0.5,0,2.5),
    oma = c(4,12,2,1))

# layout(matrix(c(1,1,1,1,2,3,3,3,
#                 1,1,1,1,2,4,4,4,
#                 1,1,1,1,2,5,5,5,
#                 1,1,1,1,2,6,6,6),
#               byrow = T, nrow = 4))

layout(matrix(c(
  1,1,1,1,3,3,3,
  2,2,2,2,3,3,3,
  2,2,2,2,4,4,4,
  2,2,2,2,4,4,4,
  2,2,2,2,5,5,5,
  2,2,2,2,5,5,5,
  2,2,2,2,6,6,6,
  2,2,2,2,6,6,6),
  byrow = T, nrow = 8))

# add color bar
color.bar.horiz(cols_to_plot_in_color_bar, min = -1)
# axis(side = 1, las = 1, at = seq(-1,1, length.out = max(modal_peak_week_in_future, na.rm = T)),
#      labels = 1:max(modal_peak_week_in_future, na.rm = T))
axis(side = 1, las = 1, at = seq(-1,1, length.out = max(modal_peak_week_range_in_future, na.rm = T) + 1),
     labels = NA)
axis(side = 1, las = 1, at = seq(-0.95,0.95, length.out = max(modal_peak_week_range_in_future, na.rm = T)),
     labels = c('1-4', '5-8', '9-12', '13-16',
                         '17-20', '21-24', '25-28', '29-32',
                         '33-36', '36-39', '40-43', '44-47', '48-51', '52-55', '56-59'),
     tick = F, las = 2
     # , labels = F
     )

mtext(side = 3, text = 'Modal peak week')

par(mar = c(1,0.5,0,2.5))


# add peak wk plot
image(image_to_plot,
      col = cols,
      xaxt = 'n', 
      yaxt = 'n' )
axis(side = 1, at = seq(0, 1, length.out = length(when_to_assimilate)),
     labels = when_to_assimilate - week_first_case, las = 2)
mtext(side = 1, text = 'Weeks of data assimilated into model',line = 3)
# mtext(side = 3, text = 'Probability that peak week is in the future')

labs = dept_names
axis(side = 2, at = seq(0, 1, length.out = length(dept_names)), las = 1, 
     labels = labs[order(peak_week_by_dept)])

# add a box for each peak week 
l = 1
for(dd in order(peak_week_by_dept)){
  if(dd != 22){
    ind_tmp = which.min(abs(when_to_assimilate - week_first_case - peak_week_by_dept[dd]))
    
    text(y = seq(0, 1, length.out = length(dept_names))[l],
           x = seq(0, 1, length.out = length(when_to_assimilate))[ind_tmp], 
          labels = modal_peak_week_in_future[dd, ind_tmp], col = 'navy')
  }
  l = l + 1
}


# add plots at peak week
par(mar = c(2,5,0.5,1))
cols_to_plot = c('darkolivegreen1', 'darkolivegreen3', 'darkolivegreen4')

depts_to_plot = c(8,9,10,32)

for(dd in depts_to_plot){
  
  if(dd != 22){
    I_A_to_plot = log10(round(I_A_CI_rho[dd,,]))
    I_A_to_plot[which(is.infinite(I_A_to_plot))] = 0
    I_A_to_plot[which(is.na(I_A_to_plot))] = 0
    
    I_A_to_plot_large = log10(round(I_A_CI_rho_large[dd,,]))
    I_A_to_plot_large[which(is.infinite(I_A_to_plot_large))] = 0
    I_A_to_plot_large[which(is.na(I_A_to_plot_large))] = 0
    
    cases_now = log10(df_dept$cases_new[which(df_dept$dept_name == dept_names[dd] & df_dept$week %in% week_first_case: max(df_dyn$week))])
    
    tt = which.min(abs(when_to_assimilate - week_first_case - peak_week_by_dept_range[dd,2]))
    l = tt
    
    I_F_to_plot = log10(round(I_F_CI_rho[[tt]][dd,,]))
    I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
    
    I_F_to_plot_large = log10(round(I_F_CI_rho_large[[tt]][dd,,]))
    I_F_to_plot_large[which(is.infinite(I_F_to_plot_large))] = 0
    
    plot(-100, -100, xlim = c(0, dim(I_A_to_plot)[2]), 
         ylim = c(0, max(I_F_to_plot_large, I_A_to_plot_large, cases_now, na.rm = T)), 
         xlab = '', ylab = '', xaxt = 'n', las = 1, yaxt = 'n', bty = 'n')
    abline(v = (time_btwn_assim * l), col = 'slateblue', lwd = 2)
    mtext(side = 3, text = dept_names[dd], cex = 0.5)
    axis(side = 2, las = 1, at = seq(0, max(I_F_to_plot_large, na.rm = T), length.out = 4),
         labels = signif(seq(0, max(I_F_to_plot_large, na.rm = T), length.out = 4), 2))
    
    # 75% CI I_A
    polygon(x = c(1:(time_btwn_assim * l),
                  rev(1:(time_btwn_assim * l))),
            y = c(I_A_to_plot_large[1, 1:(time_btwn_assim * l)],
                  rev(I_A_to_plot_large[2, 1:(time_btwn_assim * l)])),
            col = adjustcolor(cols_to_plot[1], alpha.f = 0.3), 
            border = adjustcolor(cols_to_plot[1], alpha.f = 0.3))
    
    # 50% CI I_A
    polygon(x = c(1:(time_btwn_assim * l),
                  rev(1:(time_btwn_assim * l))),
            y = c(I_A_to_plot[1, 1:(time_btwn_assim * l)],
                  rev(I_A_to_plot[2, 1:(time_btwn_assim * l)])),
            col = adjustcolor(cols_to_plot[2], alpha.f = 0.6), 
            border = adjustcolor(cols_to_plot[2], alpha.f = 0.6))
    
    # 75% CI I_F
    polygon(x = c((1:ncol(I_F_to_plot_large)) + (time_btwn_assim * l - 1),
                  rev((1:ncol(I_F_to_plot_large)) + (time_btwn_assim * l - 1))),
            y = c(I_F_to_plot_large[1,],
                  rev(I_F_to_plot_large[2,])), col = adjustcolor(cols_to_plot[1], alpha.f = 0.2), 
            border = adjustcolor(cols_to_plot[1], alpha.f = 0.8))
    
    # 50% CI I_F
    polygon(x = c((1:ncol(I_F_to_plot)) + (time_btwn_assim * l - 1),
                  rev((1:ncol(I_F_to_plot)) + (time_btwn_assim * l - 1))),
            y = c(I_F_to_plot[1,],
                  rev(I_F_to_plot[2,])), col = adjustcolor(cols_to_plot[2], alpha.f = 0.3), 
            border = adjustcolor(cols_to_plot[2], alpha.f = 0.3))
    
    lines(I_A_to_plot[3, 1:(time_btwn_assim * l)], col = cols_to_plot[3], lwd = 2.5)
    lines(c(rep(NA,(time_btwn_assim * l - 1)), I_F_to_plot[3,]) , col = cols_to_plot[3], lwd = 2.5)
    
    points(cases_now[1:(when_to_assimilate[l] - week_first_case)], cex = 0.8, pch = 4, col = rgb(0,0,0,0.8))
    points(c(rep(NA, length(1:(when_to_assimilate[l] - week_first_case))), 
             cases_now[(when_to_assimilate[l] - week_first_case):length(cases_now)]), 
           cex = 0.8, pch = 4
           , col = rgb(0,0,0,0.8)
    )
    
    
  }
  
  if(dd == 9){
    mtext(side = 2, line = 3, at = -1, text = paste0('Weekly incidence (log10)'))
  }
  
  
  if(dd %in% c(32)){
    axis(side = 1, las = 2)
    mtext(side = 1, line = 2.5, text = 'Number of weeks since first')
    mtext(side = 1, line = 4, text = 'case reported in Colombia')
  }
}


dev.off()


#=============================================================================#
# forecast for depts of focus
#=============================================================================#


f = paste0('../output/', file_path, '/depts_forecasted_to_plot', file_out_extension, '.pdf')
pdf(file = f, height = 4, width = 7)

times_to_plot = seq(1, length(when_to_assimilate), by = 4)

layout(matrix(1:(length(times_to_plot)* length(depts_to_plot)) , 
              nrow = length(times_to_plot), ncol = 4, byrow = F))
par(mar = c(0.5,3.5,0,1), oma = c(3,3,2,1))

cols_to_plot = c('darkolivegreen1', 'darkolivegreen3', 'darkolivegreen4')
for(dd in depts_to_plot){
  
  # dd = 8
  
  I_A_to_plot = log10(round(I_A_CI_rho[dd,,]))
  I_A_to_plot[which(is.infinite(I_A_to_plot))] = 0
  I_A_to_plot[which(is.na(I_A_to_plot))] = 0
  
  I_A_to_plot_large = log10(round(I_A_CI_rho_large[dd,,]))
  I_A_to_plot_large[which(is.infinite(I_A_to_plot_large))] = 0
  I_A_to_plot_large[which(is.na(I_A_to_plot_large))] = 0
  
  cases_now = log10(df_dept$cases_new[which(df_dept$dept_name == dept_names[dd] & df_dept$week %in% week_first_case: max(df_dyn$week))])
  
  
  for(tt in 1:length(when_to_assimilate)){
    
    if(tt %in% times_to_plot){
      I_F_to_plot = log10(round(I_F_CI_rho[[tt]][dd,,]))
      I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
      
      I_F_to_plot_large = log10(round(I_F_CI_rho_large[[tt]][dd,,]))
      I_F_to_plot_large[which(is.infinite(I_F_to_plot_large))] = 0
      
      plot(-100, -100, xlim = c(0, dim(I_A_to_plot)[2]), 
           # ylim = c(0, max(I_F_to_plot_large, I_A_to_plot_large, cases_now, na.rm = T)), 
           ylim = c(0, 4.5),
           xlab = '', ylab = '', xaxt = 'n', las = 1, yaxt= 'n')
      # axis(side = 2, las = 1, at = seq(0, max(I_F_to_plot_large, I_A_to_plot_large, na.rm = T), length.out = 4),
      #      labels = signif(seq(0, max(I_F_to_plot_large, I_A_to_plot_large, na.rm = T), length.out = 4), 2))
      axis(side = 2, las = 1, at = 1:4, labels= c(10, 100, 1000, 10000))
      
      abline(v = when_to_assimilate[tt] - week_first_case, col = 'slateblue', lwd = 2)
      
      if(tt == 1){
        mtext(side = 3, text = dept_names[dd], cex = 0.75)
      }
      
      
      if(tt > 1){
        # 75% CI I_A
        polygon(x = c(1:(when_to_assimilate[tt] - week_first_case),
                      rev(1:(when_to_assimilate[tt] - week_first_case))),
                y = c(I_A_to_plot_large[1, 1:(when_to_assimilate[tt] - week_first_case)],
                      rev(I_A_to_plot_large[2, 1:(when_to_assimilate[tt] - week_first_case)])),
                col = adjustcolor('darkolivegreen1', alpha.f = 0.3), 
                border = adjustcolor('darkolivegreen1', alpha.f = 0.3))
        
        # 50% CI I_A
        polygon(x = c(1:(when_to_assimilate[tt] - week_first_case),
                      rev(1:(when_to_assimilate[tt] - week_first_case))),
                y = c(I_A_to_plot[1, 1:(when_to_assimilate[tt] - week_first_case)],
                      rev(I_A_to_plot[2, 1:(when_to_assimilate[tt] - week_first_case)])),
                col = adjustcolor('darkolivegreen1', alpha.f = 0.6), 
                border = adjustcolor('darkolivegreen1', alpha.f = 0.6))
      }
      
      # 75% CI I_F
      polygon(x = c((when_to_assimilate[tt] - week_first_case) : 
                      (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)),
                    rev((when_to_assimilate[tt] - week_first_case) : 
                          (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)))),
              y = c(I_F_to_plot_large[1,],
                    rev(I_F_to_plot_large[2,])), 
              col = adjustcolor('slateblue', alpha.f = 0.1), 
              border = adjustcolor('slateblue', alpha.f = 0.1))
      
      # 50% CI I_F
      polygon(x = c((when_to_assimilate[tt] - week_first_case) : 
                      (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)),
                    rev((when_to_assimilate[tt] - week_first_case) : 
                          (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)))),
              y = c(I_F_to_plot[1,],
                    rev(I_F_to_plot[2,])), 
              col = adjustcolor('slateblue', alpha.f = 0.3), 
              border = adjustcolor('slateblue', alpha.f = 0.3))
      
      lines(I_A_to_plot[3, 1:(when_to_assimilate[tt] - week_first_case)], col = 'darkolivegreen3', lwd = 4)
      lines(x = (when_to_assimilate[tt] - week_first_case) : 
              (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)), 
            y = I_F_to_plot[3,], col = 'slateblue', lwd = 4)
      
      points(cases_now, cex = 1, pch = 4, col = rgb(0,0,0,0.8))
    }
  }
  axis(side = 1)
    
}

mtext(side = 1, line = 2, text = 'Number of weeks since first case reported in Colombia', outer = T)
mtext(side = 2, line = 1, text = 'New weekly cases', outer = T)

dev.off()

#=============================================================================#
# peak incidence forecasting
#=============================================================================#


# empirical peak incidence for each department
pk_wk_inds = sapply(1:length(dept_names), function(dd) 
  which.max(df_dept$cases_new[which(df_dept$dept_name == dept_names[dd])]))

peak_incidence_by_dept = sapply(1:length(dept_names), 
                                function(dd) 
                                  df_dept$cases_new[which(df_dept$dept_name == dept_names[dd])][pk_wk_inds[dd]])



# peak week forecast CrI
peak_incidence_forecast = list()
length(peak_incidence_forecast) = length(when_to_assimilate)
  
  
for(tt in 1:length(when_to_assimilate)){
  
  peak_incidence_forecast[[tt]] = matrix(NA, ncol = np, nrow = length(dept_names))
  
  
  for(dd in 1:length(dept_names)){
    
    if(dd != 22){
      if(file_path %in% c('no_tm_ini_5', 'no_tm_ini_6', 'no_tm_ini_7', 'no_tm_ini_8')){
        
        if(tt <= which(when_to_assimilate == when_to_assimilate_NS[[dd]][2]) - 1){
          forecast = I_F_over_time_rho_sample[[1]][dd,,]
        }else{
          forecast = cbind(I_A_over_time_rho_sample[dd, , 1:(when_to_assimilate[tt] - week_first_case)], 
                           I_F_over_time_rho_sample[[tt]][dd,,])
        }
      }else{
        if(tt > 1){
          forecast = cbind(I_A_over_time_rho_sample[dd, , 1:(when_to_assimilate[tt] - week_first_case)], 
                           I_F_over_time_rho_sample[[tt]][dd, ,])
        }else{
          forecast = I_F_over_time_rho_sample[[tt]][dd,,]
        }
      }
      
      forecasted_pk_wk_inds = sapply(1:np, function(ff) which.max(forecast[ff, ]))
    }
    
    for(pp in 1:np){
      peak_incidence_forecast[[tt]][dd, pp] = forecast[pp, forecasted_pk_wk_inds[pp]]
    }
  }
}


#=============================================================================#
# peak incidence logscoring
# to get an approximation of logscores across departments, need to take sum
#=============================================================================#

# empirical peak incidence
empirical_peak_incidence = sapply(1:length(dept_names), function(ff) df_dept$cases[which(df_dept$dept_name == dept_names[ff])][pk_wk_inds[ff]])
peak_incidence_logscore = array(NA, dim= c(length(when_to_assimilate), length(dept_names)))


for(tt in 1:(length(when_to_assimilate))){
  for(dd in 1:length(dept_names)){
    if(dd != 22){
      
      x = peak_incidence_forecast[[tt]][dd, ]
      d = density(x)
      fxn = stats::approxfun(d$x, d$y)
      peak_incidence_logscore[tt, dd] = fxn(empirical_peak_incidence[dd])
      
    }
  }
}

# setting all infinite logscores equal to the minimum logscore
peak_incidence_logscore = log10(peak_incidence_logscore)
peak_incidence_logscore[which(is.infinite(peak_incidence_logscore))] = 
  peak_incidence_logscore[-which(is.infinite(peak_incidence_logscore))][which.min(peak_incidence_logscore[-which(is.infinite(peak_incidence_logscore))])]


f = paste0('../output/', file_path, '/peak_incidence_logscore_by_dept', file_out_extension, '.pdf')
pdf(file = f, height = 6, width = 5.5)
cols = colorRampPalette(c('aliceblue', 'navy'))(25)
par(mar = c(2,2,2,2),
    oma = c(4, 10, 0,1))

layout(matrix(c(1,1,1,1,2,
                1,1,1,1,2,
                1,1,1,1,2,
                1,1,1,1,2), 
              byrow = T, ncol = 5))

image(peak_incidence_logscore, col= cols,
      xaxt = 'n', 
      yaxt = 'n' )

axis(side = 1, at = seq(0, 1, length.out = length(when_to_assimilate)),
     labels = when_to_assimilate - week_first_case)
mtext(side = 1, text = 'Weeks of data assimilated into model',line = 3)
mtext(side = 3, text = 'Peak incidence')

labs = dept_names
axis(side = 2, at = seq(0, 1, length.out = length(dept_names)), las = 1, 
     labels = labs)

# add a box for each 
for(dd in 1:length(dept_names)){
  if(dd != 22){
    ind_tmp = which.min(abs(when_to_assimilate - week_first_case - peak_week_by_dept_range[dd,2]))
    points(y = seq(0, 1, length.out = length(dept_names))[dd],
           x = seq(0, 1, length.out = length(when_to_assimilate))[ind_tmp], 
           col = 'deeppink2', pch = 19)
  }
}

color.bar(cols, min = -1)
mtext(side = 1, text = 'Less accurate')
mtext(side = 3, text = 'More accurate')
mtext(side = 2, text = 'Logscore')
dev.off()




#=============================================================================#
# week when 10 cumulative cases have occurred 
#=============================================================================#

if(file_path %in% c('no_tm_ini_5', 'no_tm_ini_6', 'no_tm_ini_7', 'no_tm_ini_8')){
  print('nope')
}else{
  # week of first case forecast
  week_of_first_10_cases_forecast = list()
  length(week_of_first_10_cases_forecast) = length(when_to_assimilate)
  
  for(tt in 1:length(when_to_assimilate)){
    
    week_of_first_10_cases_forecast[[tt]] = matrix(NA, ncol = np, nrow = length(dept_names))
    
    for(dd in 1:length(dept_names)){
      
      if(tt > 1){
        forecast = cbind(I_A_over_time_rho_sample[dd, , 1:(when_to_assimilate[tt] - week_first_case)], 
                         I_F_over_time_rho_sample[[tt]][dd, ,])
      }else{
        forecast = I_F_over_time_rho_sample[[tt]][dd,,]
      }
      
      week_of_first_10_cases_forecast[[tt]][dd, ] = sapply(1:np, function(ff) which(cumsum(forecast[ff,]) >= 10)[1])
      # week_of_first_10_cases_forecast[[l]][dd,][which(is.na(week_of_first_10_cases_forecast[[l]][dd, ]))] = 0 
    }
  }
  
  
  # empirical week by which 10 cases occurred
  week_of_first_10_cases = rep(0, length(dept_names))
  for(dd in 1:length(dept_names)){
    df_tmp = df_dyn[which(df_dyn$dept_name == dept_names[dd]),]
    week_of_first_10_cases[dd] = df_tmp$week[which(cumsum(df_tmp$cases_new)>=10)][1] - week_first_case
  }
  week_of_first_10_cases[which(is.na(week_of_first_10_cases))] = 60
  
  
  # only look at first 3 forecasts
  f_out = paste0('../output/', file_path, '/week_by_which_10_cases_occurred', file_out_extension, '.pdf')
  pdf(file = f_out, height = 7, width = 10)
  xs_to_plot = sort(week_of_first_10_cases, index.return = T)$ix
  layout(matrix(1:6, ncol = 3, byrow = T))
  par(mar = c(1,1,2,1), oma = c(8, 7, 2, 1))
  
  for(tt in 1:6){
    
    plot(-100, -100, 
         ylim = c(min(week_of_first_10_cases, na.rm = T, week_of_first_10_cases_forecast[[tt]]),
                  max(week_of_first_10_cases, na.rm = T, week_of_first_10_cases_forecast[[tt]])),
         xlim = c(1, length(dept_names)),
         xaxt = 'n', yaxt = 'n', ylab = '', xlab = '')
    
    txt_tmp = paste0(when_to_assimilate[tt] - week_first_case, ' wks after first case in Colombia')
    mtext(side = 3, text = txt_tmp, line = 0.5)
    
    if(tt == 1){
      mtext(side = 2, line = 5, outer = T, text = 'Week by which 10 cases occurred')
    }
    
    if(tt == 6){
      legend('topleft', legend = c('10 cases in dept. have occurred',
                                   '10 cases in dept. have not occurred',
                                   '10 empirical cases occurred in dept.'),
             fill = c('darkolivegreen3', 'darkslategray3', NA), 
             border = c('darkolivegreen3', 'darkslategray3', NA),
             pch = c(NA, NA, 16),
             col = c(NA, NA, 'blueviolet'))
    }
    
    if(tt > 3){
      axis(side = 1, at = 1:32, labels = dept_names[xs_to_plot], las = 2, cex.axis = 0.75)
    }
    
    if(tt %in% c(1, 4)){
      start_of_week = c(yearSunday(2015), yearSunday(2016))[32 : (32 + 60)]
      labs = paste(month.abb[month(start_of_week)], day(start_of_week), year(start_of_week))
      axis(side = 2, las = 1, at = seq(0, 60, by = 10), labels = labs[c(1, 10, 20, 30, 40, 50, 60)])
    }
    
    for(dd in 1:length(dept_names)){
      if(dd != 22){
        
        if(week_of_first_10_cases[dd] > (when_to_assimilate[tt] - week_first_case)){
          boxplot(x = week_of_first_10_cases_forecast[[tt]][dd, ], at = which(xs_to_plot == dd), add = T,
                  col = adjustcolor('darkslategray3', alpha.f = 0.8)
                  ,
                  outcol = adjustcolor('darkslategray3', alpha.f = 0.1)
                  , outpch = 20
                  , border = adjustcolor('darkslategray4', alpha.f = 0.9)
                  , medcol = 'navy'
                  , pars = list(xaxt = 'n', yaxt = 'n', axes = F)
          )
        }else{
          boxplot(x = week_of_first_10_cases_forecast[[tt]][dd, ], at = which(xs_to_plot == dd), add = T,
                  col = adjustcolor('darkolivegreen3', alpha.f = 0.8)
                  ,
                  outcol = adjustcolor('darkolivegreen3', alpha.f = 0.1)
                  , outpch = 20
                  , border = adjustcolor('darkolivegreen4', alpha.f = 0.9)
                  , medcol = 'navy'
                  , pars = list(xaxt = 'n', yaxt = 'n', axes = F)
          )
        }
        
        segments(x0 =  which(xs_to_plot == dd) - 0.2,
                 x1 = which(xs_to_plot == dd) + 0.2,
                 y0 =  week_of_first_10_cases[dd], 
                 y1 = week_of_first_10_cases[dd],
                 lwd = 3, col = 'blueviolet')
        
      }
    }
  }
  
  dev.off()
}
