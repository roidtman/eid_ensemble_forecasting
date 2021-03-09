#=============================================================================#
# Author: Rachel Oidtman

# ANALYSIS: produce forecasting EM weights through time

# Supplementary Figures
# Log scores through time and departments (multi-panel)
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


#=============================================================================#
# assess logscores for different models at different points in time
# for example: when first 4 weeks have been assimilated into the model, 
# assess how well the model performed at weeks 1:4
#=============================================================================#


# AGGREGATED, CUMULATIVE NATIONAL INCIDENCE 
LS_national_by_model_over_time = array(NA, dim = c(length(runs),
                                                   length(when_to_assimilate)))

# CUMULATIVE DEPARTMENT INCIDENCE
LS_dept_by_model_over_time = array(NA, dim = c(length(runs),
                                               length(when_to_assimilate),
                                               length(dept_names)))


for(rr in 1:length(runs)){
  
  eval(parse(text = paste0(
    'I_F_over_time_rho_sample = I_F_over_time_rho_sample_', runs[rr]
  )))
  
  eval(parse(text = paste0(
    'I_A_over_time_rho_sample = I_A_over_time_rho_sample_', runs[rr]
  )))
  
  for(tt in 1:length(when_to_assimilate)){
    # for(tt in 2:length(when_to_assimilate)){
    for(dd in dept_models){
      
      if(rr %in% 5:8){
        
        if(which(when_to_assimilate == when_to_assimilate_NS[[dd]][2]) - 1 <= tt){
          # forecast_current = round(rowSums(
          #   I_A_over_time_rho_sample[dd,,wks_assimilated_at_when_to_assimilate[tt,]], na.rm = T))
          
          forecast_current = round(rowSums(
            I_F_over_time_rho_sample[[tt]][dd, , 1:time_btwn_assim], na.rm = T))
          
          dat_current = sum(df_dyn$cases_new[which(
            df_dyn$week %in% (wks_assimilated_at_when_to_assimilate[tt,] + week_first_case - 1) & 
              df_dyn$dept_name == dept_names[dd])])
          
          LS_dept_by_model_over_time[rr, tt, dd] = incidence_LS(distribution_in = forecast_current, 
                                                                data_in = dat_current)
        }
      }else{
        # forecast_current = round(rowSums(
        #   I_A_over_time_rho_sample[dd,,wks_assimilated_at_when_to_assimilate[tt,]], na.rm = T))
        
        forecast_current = round(rowSums(
          I_F_over_time_rho_sample[[tt]][dd, , 1:time_btwn_assim], na.rm = T))
        
        dat_current = sum(df_dyn$cases_new[which(
          df_dyn$week %in% (wks_assimilated_at_when_to_assimilate[tt,] + week_first_case - 1) & 
            df_dyn$dept_name == dept_names[dd])])
        
        LS_dept_by_model_over_time[rr, tt, dd] = incidence_LS(distribution_in = forecast_current, 
                                                              data_in = dat_current)
      }
    }
  }
  for(tt in 1:length(when_to_assimilate)){
    # for(tt in 2:length(when_to_assimilate)){
    
    if(rr %in% 5:8){
      forecast_NS = array(NA, dim = c(length(dept_names), np, time_btwn_assim))
      for(dd in 1:length(dept_names)){
        if(dd != 22){
          if(which(when_to_assimilate == when_to_assimilate_NS[[dd]][2]) - 1 <= tt){
            forecast_NS[dd,,] = I_F_over_time_rho_sample[[tt]][dd,, 1:time_btwn_assim]
          }else{
            forecast_NS[dd,,] = I_F_over_time_rho_sample[[1]][dd,, 1:time_btwn_assim]
          }
          forecast_current = sapply(1:np, function(ff) sum(forecast_NS[,ff,], na.rm = T))
          dat_current = sum(df_dyn$cases_new[which(
            df_dyn$week %in% (wks_assimilated_at_when_to_assimilate[tt,] + week_first_case - 1))])
          
          LS_national_by_model_over_time[rr, tt] = incidence_LS(distribution_in = forecast_current,
                                                                data_in = dat_current)
        }
      }
      
    }else{
      # forecast_current = sapply(1:np, function(ff)
      #   sum(I_A_over_time_rho_sample[,ff,wks_assimilated_at_when_to_assimilate[tt,]], na.rm = T))
      
      forecast_current = sapply(1:np, function(ff)
        sum(I_F_over_time_rho_sample[[tt]][, ff, 1:time_btwn_assim], na.rm = T))
      
      dat_current = sum(df_dyn$cases_new[which(
        df_dyn$week %in% (wks_assimilated_at_when_to_assimilate[tt,] + week_first_case - 1))])
      
      LS_national_by_model_over_time[rr, tt] = incidence_LS(distribution_in = forecast_current,
                                                            data_in = dat_current)
      
    }
  }
}  

# REMOVE BOGOTA
dim(LS_dept_by_model_over_time)
LS_dept_by_model_over_time = LS_dept_by_model_over_time[,,-22]
dim(LS_dept_by_model_over_time)

#=============================================================================#
# assess ensemble weights at different points in time
#=============================================================================#

# assigning 1000-fold reduction for -infinity and NA values
LS_dept_by_model_over_time[which(is.infinite(LS_dept_by_model_over_time))] = 
  min(LS_dept_by_model_over_time[-which(is.infinite(LS_dept_by_model_over_time))], na.rm = T) * 1e4
LS_dept_by_model_over_time[which(is.na(LS_dept_by_model_over_time))] = min(LS_dept_by_model_over_time, na.rm = T)


# LS_national_by_model_over_time[which(is.infinite(LS_national_by_model_over_time))] = min(LS_national_by_model_over_time[-which(is.infinite(LS_national_by_model_over_time))], na.rm = T) * 1e4
LS_national_by_model_over_time[which(is.na(LS_national_by_model_over_time))] = 
  min(LS_national_by_model_over_time, na.rm = T) * 1e4


# one object to hold all log scores
LS_over_time = array(NA, dim = c(length(runs), length(when_to_assimilate), length(dept_models) + 1))

for(ii in 1:length(when_to_assimilate)){
  LS_over_time[,ii, 1:length(dept_models)] = LS_dept_by_model_over_time[,ii,]
  LS_over_time[,ii,(length(dept_models) + 1)] = LS_national_by_model_over_time[,ii]
}
LS_over_time_all = LS_over_time


ensemble_wts_through_time = matrix(NA, nrow = length(runs), ncol = length(when_to_assimilate))

for(ii in 1:length(when_to_assimilate)){
  ensemble_wts_through_time[, ii] = produce_ensemble_wts(num_models = length(runs),
                                                         LS_in = LS_over_time[,ii,], browse = F)
}

save(ensemble_wts_through_time, 
     file = '../output/ensemble/ensemble_weights_through_time_I-F.RData')


#=============================================================================#
# figures
#=============================================================================#

dept_names_to_plot = c('La Guajira', 'Magdalena', 'Atlántico', 'Cesar', 'Bolívar',
                       'Sucre', 'Córdoba', 'Norte de Santander', 'Antioquia', 
                       'Chocó', 'Santander', 'Arauca', 'Boyacá', 'Vichada', 
                       'Casanare', 'Cudinamarca', 'Caldas', 'Risaralda', 'Tolima', 
                       'Valle del Cauca', 'Meta', 'Bogota', 'Quindío', 'Guainía', 'Huila', 
                       'Cauca', 'Caquetá', 'Guaviare', 'Nariño', 'Vaupés', 
                       'Putumayo', 'Amazonas')



# pdf('../output/ensemble/ensemble_weights_through_time.pdf', height = 5, width = 5)
pdf('../output/main_text_figures/supplement/ensemble_weights_through_time.pdf', height = 5, width = 5)

model_names = c('R-1-CDRs', 'R-2-CDRs', 'Rt-1-CDRs', 'Rt-2-CDRs',
                'R-1-nonspatial', 'R-2-nonspatial', 'Rt-1-nonspatial', 'Rt-2-nonspatial',
                'R-1-gravity', 'R-2-gravity', 'Rt-1-gravity', 'Rt-2-gravity', 
                'R-1-radiation', 'R-2-radiation', 'Rt-1-radiation', 'Rt-2-radiation')

cols = colorRampPalette(c('aliceblue', 'navy'))(25)
par(mar = c(2,2,2,2),
    oma = c(4, 6, 0,1))

layout(matrix(c(1,1,1,1,2,
                1,1,1,1,2,
                1,1,1,1,2,
                1,1,1,1,2), 
              byrow = T, ncol = 5))

image(t(ensemble_wts_through_time), col= cols,
      xaxt = 'n', 
      yaxt = 'n' )
axis(side = 1, at = seq(0, 1, length.out = length(when_to_assimilate)),
     labels = when_to_assimilate - week_first_case, las = 2)
mtext(side = 1, text = 'Weeks of data assimilated into model',line = 3)
mtext(side = 3, text = 'Ensemble weights through time')


labs = paste0('Model ', runs)
axis(side = 2, at = seq(0, 1, length.out = length(runs)), las = 1, 
     labels = model_names)

par(mar = c(0.25,2,0.25,2))
color.bar(cols, min = -1)
axis(side = 4, at = seq(-1, 1, length.out = length(cols)), las = 1, 
     labels = sprintf(seq(0, max(ensemble_wts_through_time, na.rm = T), length.out = 25), fmt = '%#.2f'))
mtext(side = 2, text = 'Ensemble weight')

dev.off()


#=============================================================================#
# figures
#=============================================================================#

# log scores by
library(viridis)


pdf('../output/main_text_figures/supplement/median_log_score_real-time_analysis.pdf',
    height = 5, width = 6)
par(mar = c(3,3,1, 0), oma = c(1,3,1,7))
layout(matrix(c(1,1,1,1,3,
                1,1,1,1,3,
                2,2,2,2,3,
                2,2,2,2,3), 
              byrow = T, ncol = 5))
cols = magma(16)

# department-specific
plot(-100, 100, xlim = c(1, 15), ylim = c(-5, 0), las = 1, xaxt = 'n')
axis(side = 1, at = 1:length(when_to_assimilate), labels = when_to_assimilate - week_first_case)
for(ii in 1:15){
  for(dd in 1:16){
    if(dd %in% 7:8){
      points(x = ii, y = median(LS_dept_by_model_over_time[dd,ii,], na.rm = T), cex = 1.25,
             pch = 17, col = cols[dd]) 
    }else if(dd %in% 5:6){
      points(x = ii, y = median(LS_dept_by_model_over_time[dd,ii,], na.rm = T), cex = 1.25,
             pch = 15, col = cols[dd]) 
    }else{
      points(x = ii, y = median(LS_dept_by_model_over_time[dd,ii,], na.rm = T), cex = 1.25,
             pch = 16, col = cols[dd])
    }
  }
}
mtext(side = 2, text = 'Median (over departments) log score', cex = 0.8, line = 2.25)

# national-level
plot(-100, 100, xlim = c(1, 15), ylim = c(-7, 0), las = 1, xaxt = 'n')
axis(side = 1, at = 1:length(when_to_assimilate), labels = when_to_assimilate - week_first_case)
for(ii in 1:15){
  for(dd in 1:16){
    if(dd %in% 7:8){
      points(x = ii, y =LS_national_by_model_over_time[dd,ii], cex = 1.25,
             pch = 17, col = cols[dd]) 
    }else if(dd %in% 5:6){
      points(x = ii, y = LS_national_by_model_over_time[dd,ii], cex = 1.25,
             pch = 15, col = cols[dd]) 
    }else{
      points(x = ii, y = LS_national_by_model_over_time[dd,ii], cex = 1.25,
             pch = 16, col = cols[dd])
    }
  }
}
mtext(side = 2, text = 'National log score', cex = 0.8, line = 2.25)

mtext(side = 1, text = 'Week of forecast (relative to wk. of first case)', cex = 0.8, line = 2.25)
color.bar(c(cols, 'white'), min = -1)
axis(side = 4, at = seq(-0.95, 0.95, length.out = length(cols)), las = 1, 
     labels = model_names)
dev.off()


#=============================================================================#
# figures
#=============================================================================#

# plot x axis ordered by cumulative incidence in each department
cum_cases_by_dept = sapply(1:length(dept_names), 
                           function(ff) 
                             sum(df_dept$cases_new[which(df_dept$dept_name == dept_names[ff])][
                               week_first_case : max(df_dyn$week)]))
cum_cases_by_dept = cum_cases_by_dept[-22]
xs_to_plot = c(order(cum_cases_by_dept), 32)

# attack rates
pop_by_dept = sapply(1:length(dept_names), function(ff) 
  df_dept$pop[which(df_dept$dept_name == dept_names[ff])][1])
pop_by_dept = pop_by_dept[-22]

cum_cases_by_dept / pop_by_dept


pdf('../output/main_text_figures/supplement/logscore_by_assimilation_period_by_model_by_dept.pdf',
    height = 5, width = 11)
# cols = colorRampPalette('deeppink4', 'white', 'darkcyan')(40)
cols = colorRampPalette(c('aliceblue', 'navy'))(40)
layout(matrix(c(1,2,3,4,5,16,
                6,7,8,9,10,16,
                11,12,13,14,15,16), byrow = T, nrow = 3))
# layout(matrix(1:15, 3, 5, byrow = T))
par(mar = c(2,0.5,1,0), oma = c(3,5,1,0))
for(tt in 1:15){
  image(t(exp(LS_over_time_all[,tt,rev(xs_to_plot)])), col= cols,
        xaxt = 'n', 
        yaxt = 'n' )
  mtext(side = 3, text = when_to_assimilate[tt] - week_first_case, cex = 0.75)
  mtext(side = 3, text = letters[tt], 
        at = 0, cex = 0.75)
  
  if(tt %in% 11:15){
    axis(side = 1, seq(0.0, 1, length.out = 32), 
         labels = c(dept_names_to_plot[-22], 'National')[rev(xs_to_plot)],
         las = 2, cex.axis=0.4)
  }else{
    axis(side = 1, seq(0.0, 1, length.out = 32), labels = NA)
  }
  
  if(tt %in% c(1,6,11)){
    axis(side = 2, at = seq(0, 1, length.out = length(runs)),
         labels = model_names, las = 1, cex.axis = 0.4)
  }
  
}
par(mar = c(0.25,2,0.25,8))
color.bar(c(cols, 'white'), min = -1)
mtext(side = 3, text = 'log scores', line = -1, cex = 0.75)
mtext(side = 3, text = '4-wk ahead forecast', line = 0, at = 13, cex = 0.75)
axis(side = 4, at = c(-0.975, 0.975), labels = c('Least accurate', 'Most accurate'), las = 1)
dev.off()

