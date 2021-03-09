#=============================================================================#
# Author: Rachel Oidtman

# FIGURE 4 and 5
#=============================================================================#


#=============================================================================#
# load in processed epi data
#=============================================================================#

load('../data/processed/processed_epi_data.RData')

# assign file path
runs = 1:16
file_paths = paste0('no_tm_ini_',runs)

load('../output/main_text_figures/output_for_forecast_target_figure.RData')

#=============================================================================#
# load in processed forecast data
#=============================================================================#

time_btwn_assim = 4
when_to_assimilate = seq(week_first_case + time_btwn_assim, max(df_dyn$week), by = time_btwn_assim)
when_to_assimilate = c(week_first_case, when_to_assimilate)

np = 2000

wks_assimilated_at_when_to_assimilate = matrix(NA, nrow = length(when_to_assimilate),
                                               ncol = time_btwn_assim)
for(ii in 1:length(when_to_assimilate)){
  wks_assimilated_at_when_to_assimilate[ii, ] = (when_to_assimilate[ii] - week_first_case + 1) : (when_to_assimilate[ii] - week_first_case + time_btwn_assim)
}

runs = c(runs, 'ensemble')

dept_names_to_plot = c('La Guajira', 'Magdalena', 'Atlántico', 'Cesar', 'Bolívar',
                       'Sucre', 'Córdoba', 'Norte de Santander', 'Antioquia', 
                       'Chocó', 'Santander', 'Arauca', 'Boyacá', 'Vichada', 
                       'Casanare', 'Cudinamarca', 'Caldas', 'Risaralda', 'Tolima', 
                       'Valle del Cauca', 'Meta', 'Quindío', 'Guainía', 'Huila', 
                       'Cauca', 'Caquetá', 'Guaviare', 'Nariño', 'Vaupés', 
                       'Putumayo', 'Amazonas')

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
                           col_vec_length = 100, 
                           add_in = 0.001, by_in = 0.001, 
                           browse = F){
  
  if(browse){browser()}
  
  less_than_zero_prop = length(seq(min(vec_in, na.rm = T), zero_pt_in - add_in, by = by_in))
  greater_than_zero_prop = length(seq(zero_pt_in + add_in, max(vec_in, na.rm = T), by = by_in))
  
  colfunc_less = colorRampPalette(c(col_less_in, col_middle_in))
  colfunc_greater = colorRampPalette(c(col_middle_in, col_greater_in))
  
  cols = c(colfunc_less(less_than_zero_prop + 1), colfunc_greater(greater_than_zero_prop + 1))
  
  return(cols)
}


#=============================================================================#
# Forecast scores averaged over time
# FIGURE 4
#=============================================================================#

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
par(mar = c(1,6,2,1), oma = c(12,4,1,0))
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
axis(side = 1, at = seq(0, 1, length.out = 31), labels = dept_names_to_plot[rev(xs_to_plot)], las = 2) 
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


#=============================================================================#
# Forecast scores averaged over departments
# FIGURE 5
#=============================================================================#

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
forecast_scores_by_target = array(NA, dim = c(num_targets, length(runs), length(when_to_assimilate)))

for(tt in 1:3){
  if(tt == 1){
    target = logged_peak_week_logscore
  }else if(tt == 2){
    target = peak_incidence_logscore_across_models
  }else if(tt == 3){
    target = logged_onset_week_logscore
  }
  
  for(ww in 1:length(when_to_assimilate)){
    for(mm in 1:length(runs)){
      forecast_scores_by_target[tt, mm, ww] = sum(target[mm,ww,]) / length(dept_names)
    }
  }
}

forecast_scores_by_target = exp(forecast_scores_by_target)


# forecast score AVERAGED ACROSS TARGETS

model_names = c('R-1-CDRs', 'R-2-CDRs', 'Rt-1-CDRs', 'Rt-2-CDRs',
                'R-1-nonspatial', 'R-2-nonspatial', 'Rt-1-nonspatial', 'Rt-2-nonspatial',
                'R-1-gravity', 'R-2-gravity', 'Rt-1-gravity', 'Rt-2-gravity',
                'R-1-radiation', 'R-2-radiation', 'Rt-1-radiation', 'Rt-2-radiation',
                'Ensemble')

# look at forecast scores relative to baseline ensemble model : target specific
average_FS_relative_to_baseline_targets = array(NA, dim = c(num_targets, 
                                                            length(model_names), 
                                                            length(when_to_assimilate)))
for(tt in 1:num_targets){
  for(mm in 1:length(runs)){
    average_FS_relative_to_baseline_targets[tt,mm,] = forecast_scores_by_target[tt,mm,] - forecast_scores_by_target[tt,17,]
  }
}



pdf('../output/main_text_figures/relative_forecast_score_time.pdf', 
    height = 7.5, width = 6)
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
  
  if(tt == 2){
    cols = center_data_fxn(vec_in = as.vector(average_FS_relative_to_baseline_targets[tt,,]), browse = F,
                           col_less_in = 'deeppink4', col_greater_in = 'darkcyan', add_in = 0.000001,
                           by_in = 0.000001)
  }else{
    cols = center_data_fxn(vec_in = as.vector(average_FS_relative_to_baseline_targets[tt,,]), browse = F,
                           col_less_in = 'deeppink4', col_greater_in = 'darkcyan',
                           add_in = 0.00001, by_in = 0.00001)
  }
  
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
mtext(side = 1, line = 3, text = "Weeks since 1st reported case in Colombia", cex = 0.8)
mtext(side = 2, line = 2, text = 'Models ordered by performance in each metric (low to high)', outer = T, cex = 0.8)


# add colorbars
par(mar = c(1,0,2,5))
for(tt in 1:num_targets){
  if(tt == 2){
    cols = center_data_fxn(vec_in = as.vector(average_FS_relative_to_baseline_targets[tt,,]), browse = F,
                           col_less_in = 'deeppink4', col_greater_in = 'darkcyan', add_in = 0.000001,
                           by_in = 0.000001)
  }else{
    cols = center_data_fxn(vec_in = as.vector(average_FS_relative_to_baseline_targets[tt,,]), browse = F,
                           col_less_in = 'deeppink4', col_greater_in = 'darkcyan',
                           add_in = 0.00001, by_in = 0.00001)
  }
  color.bar(cols, min = -1)
  if(tt == 1){
    mtext(side = 3, text = '% change from', cex = 0.8, line = 1)
    mtext(side = 3, text = 'ensemble', line = 0, cex = 0.8)
  }
  
  y_0 = seq(-1, 1, length.out = length(as.vector(average_FS_relative_to_baseline_targets[tt,,])))[min(which(as.vector(average_FS_relative_to_baseline_targets[tt,,])[order(as.vector(average_FS_relative_to_baseline_targets[tt,,]))] == 0))]
  
  axis(side =4, labels = NA, lwd.ticks = 0)
  
  if(tt == 2){
    axis(side = 4, at = 1, label = signif(max(average_FS_relative_to_baseline_targets[tt,,]) * 100, 1), 
         las = 1)
    axis(side = 4, at = -1, label = signif(min(average_FS_relative_to_baseline_targets[tt,,]) * 100, 1), 
         las = 1)
  }else{
    axis(side = 4, at = 1, label = round(max(average_FS_relative_to_baseline_targets[tt,,]) * 100), 
         las = 1)
    axis(side = 4, at = -1, label = round(min(average_FS_relative_to_baseline_targets[tt,,]) * 100), 
         las = 1)
  }
}

dev.off()


