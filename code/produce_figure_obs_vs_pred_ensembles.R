#=============================================================================#
# Author: Rachel Oidtman

# Supplementary Figure
# Observed verus predicted values for EW-ensemble and EM-ensemble forecast
#=============================================================================#


#=============================================================================#
# load in processed epi data
#=============================================================================#

require(stats)
require(RColorBrewer)

load('../data/processed/processed_epi_data.RData')
load('../data/processed/temp_mosq_data_frame.RData')

# assign file path
runs = 1:16
file_paths = paste0('no_tm_ini_',runs)

# load weights through time
load('../output/ensemble/ensemble_weights_through_time_I-F.RData')

# load post hoc weights
load('../output/ensemble/ensemble_weights_post_hoc.RData')

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

library(viridis)

np = dim(I_A_over_time_rho_sample_1)[2]


wks_assimilated_at_when_to_assimilate = matrix(NA, nrow = length(when_to_assimilate),
                                               ncol = time_btwn_assim)
for(ii in 1:length(when_to_assimilate)){
  wks_assimilated_at_when_to_assimilate[ii, ] = (when_to_assimilate[ii] - week_first_case + 1) : (when_to_assimilate[ii] - week_first_case + time_btwn_assim)
}

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


#=============================================================================#
# EM-weighted forecast
#=============================================================================#

num_weights = 2e4

I_A = array(NA, dim = c(length(dept_names),
                        num_weights,
                        dim(I_F_over_time_rho_sample_1[[1]])[3]))

I_F = list()
length(I_F) = length(when_to_assimilate)

for(tt in 1:length(when_to_assimilate)){
  
  I_F[[tt]] = array(NA, dim = c(length(dept_names),
                                num_weights,
                                length((wks_assimilated_at_when_to_assimilate[tt,1]) : 
                                         dim(I_A)[3])))
  
  num_samps_per_model = round(num_weights * ensemble_wts_through_time[,tt])
  models_with_weight = which(num_samps_per_model != 0)
  
  for(rr in 1:length(models_with_weight)){
    print(rr)
    samps_tmp = sample(1:np, num_samps_per_model[models_with_weight[rr]], replace = T)
    
    eval(parse(text = paste0(
      'I_A_over_time_rho_sample = I_A_over_time_rho_sample_', runs[models_with_weight[rr]]
    )))
    
    eval(parse(text = paste0(
      'I_F_over_time_rho_sample = I_F_over_time_rho_sample_', runs[rr]
    )))
    
    if(rr == 1){
      # try tt instead of tt + 1 here if things look weird
      I_F[[tt]][,1:num_samps_per_model[models_with_weight[rr]],] = 
        I_F_over_time_rho_sample[[tt]][, samps_tmp, ]
      
      I_A[, 1:num_samps_per_model[models_with_weight[rr]], 
          wks_assimilated_at_when_to_assimilate[tt,]] = 
        I_A_over_time_rho_sample[, samps_tmp, wks_assimilated_at_when_to_assimilate[tt,]]
    }else{
      
      I_F[[tt]][, (num_samps_per_model[models_with_weight[rr - 1]] + 1) : 
                  (num_samps_per_model[models_with_weight[rr - 1]] + 
                     num_samps_per_model[models_with_weight[rr]]), ] = 
        I_F_over_time_rho_sample[[tt]][, samps_tmp, ]
      
      I_A[, (num_samps_per_model[models_with_weight[rr - 1]] + 1) : 
            (num_samps_per_model[models_with_weight[rr - 1]] + 
               num_samps_per_model[models_with_weight[rr]]), 
          wks_assimilated_at_when_to_assimilate[tt,]] = 
        I_A_over_time_rho_sample[, samps_tmp, wks_assimilated_at_when_to_assimilate[tt,]]
    }
    
  }
}

## CrI
I_F_CI_rho = list()
length(I_F_CI_rho) = length(when_to_assimilate)
I_F_CI_rho_large = list()
length(I_F_CI_rho_large) = length(when_to_assimilate)

for(tt in 1:length(when_to_assimilate)){
  
  I_F_CI_rho[[tt]] = array(NA, dim = c(length(dept_names), 3, length(when_to_assimilate[tt]:
                                                                       max(df_dyn$week))))
  I_F_CI_rho_large[[tt]] = array(NA, dim = c(length(dept_names), 3, length(when_to_assimilate[tt]:
                                                                             max(df_dyn$week))))
  
  for(dd in 1:length(dept_names)){
    if(dd != 22){
      I_F_CI_rho[[tt]][dd,,] = write_CI(I_F[[tt]][dd,,])
      I_F_CI_rho_large[[tt]][dd,,] = write_CI(I_F[[tt]][dd,,], lwr = 0.125, upr = 0.875)
    }
  }
}

# assigning new names
I_F_CI_rho_EM = I_F_CI_rho
I_F_CI_rho_large_EM = I_F_CI_rho_large

I_A_CI_rho = array(NA, dim = c(length(dept_names), 3,  length(week_first_case:max(df_dyn$week))))
I_A_CI_rho_large = array(NA, dim = c(length(dept_names), 3, length(week_first_case:max(df_dyn$week))))
for(dd in 1:length(dept_names)){
  
  if(dd != 22){
    I_A_CI_rho[dd,,] = write_CI(I_A[dd,,])
    I_A_CI_rho_large[dd,,] = write_CI(I_A[dd,,], lwr = 0.125, upr = 0.875)
  }
}


# assigning new names
I_A_CI_rho_EM = I_A_CI_rho
I_A_CI_rho_large_EM = I_A_CI_rho_large

#=============================================================================#
# equally-weighted forecast
#=============================================================================#

ensemble_wts_through_time = matrix(1, nrow = length(runs), ncol = length(when_to_assimilate))
ensemble_wts_through_time = prop.table(ensemble_wts_through_time, 1)

num_weights = 2e4

I_A = array(NA, dim = c(length(dept_names),
                        num_weights,
                        dim(I_F_over_time_rho_sample_1[[1]])[3]))

I_F = list()
length(I_F) = length(when_to_assimilate)

for(tt in 1:length(when_to_assimilate)){
  
  I_F[[tt]] = array(NA, dim = c(length(dept_names),
                                num_weights,
                                length((wks_assimilated_at_when_to_assimilate[tt,1]) : 
                                         dim(I_A)[3])))
  
  num_samps_per_model = round(num_weights * ensemble_wts_through_time[,tt])
  models_with_weight = which(num_samps_per_model != 0)
  
  for(rr in 1:length(models_with_weight)){
    print(rr)
    samps_tmp = sample(1:np, num_samps_per_model[models_with_weight[rr]], replace = T)
    
    eval(parse(text = paste0(
      'I_A_over_time_rho_sample = I_A_over_time_rho_sample_', runs[models_with_weight[rr]]
    )))
    
    eval(parse(text = paste0(
      'I_F_over_time_rho_sample = I_F_over_time_rho_sample_', runs[rr]
    )))
    
    if(rr == 1){
      # try tt instead of tt + 1 here if things look weird
      I_F[[tt]][,1:num_samps_per_model[models_with_weight[rr]],] = 
        I_F_over_time_rho_sample[[tt]][, samps_tmp, ]
      
      I_A[, 1:num_samps_per_model[models_with_weight[rr]], 
          wks_assimilated_at_when_to_assimilate[tt,]] = 
        I_A_over_time_rho_sample[, samps_tmp, wks_assimilated_at_when_to_assimilate[tt,]]
    }else{
      
      I_F[[tt]][, (num_samps_per_model[models_with_weight[rr - 1]] + 1) : 
                  (num_samps_per_model[models_with_weight[rr - 1]] + 
                     num_samps_per_model[models_with_weight[rr]]), ] = 
        I_F_over_time_rho_sample[[tt]][, samps_tmp, ]
      
      I_A[, (num_samps_per_model[models_with_weight[rr - 1]] + 1) : 
            (num_samps_per_model[models_with_weight[rr - 1]] + 
               num_samps_per_model[models_with_weight[rr]]), 
          wks_assimilated_at_when_to_assimilate[tt,]] = 
        I_A_over_time_rho_sample[, samps_tmp, wks_assimilated_at_when_to_assimilate[tt,]]
    }
    
  }
}

## CrI
I_F_CI_rho = list()
length(I_F_CI_rho) = length(when_to_assimilate)
I_F_CI_rho_large = list()
length(I_F_CI_rho_large) = length(when_to_assimilate)

for(tt in 1:length(when_to_assimilate)){
  
  I_F_CI_rho[[tt]] = array(NA, dim = c(length(dept_names), 3, length(when_to_assimilate[tt]:
                                                                       max(df_dyn$week))))
  I_F_CI_rho_large[[tt]] = array(NA, dim = c(length(dept_names), 3, length(when_to_assimilate[tt]:
                                                                             max(df_dyn$week))))
  
  for(dd in 1:length(dept_names)){
    if(dd != 22){
      I_F_CI_rho[[tt]][dd,,] = write_CI(I_F[[tt]][dd,,])
      I_F_CI_rho_large[[tt]][dd,,] = write_CI(I_F[[tt]][dd,,], lwr = 0.125, upr = 0.875)
    }
  }
}

# assigning new names
I_F_CI_rho_equal = I_F_CI_rho
I_F_CI_rho_large_equal = I_F_CI_rho_large

I_A_CI_rho = array(NA, dim = c(length(dept_names), 3,  length(week_first_case:max(df_dyn$week))))
I_A_CI_rho_large = array(NA, dim = c(length(dept_names), 3, length(week_first_case:max(df_dyn$week))))
for(dd in 1:length(dept_names)){
  
  if(dd != 22){
    I_A_CI_rho[dd,,] = write_CI(I_A[dd,,])
    I_A_CI_rho_large[dd,,] = write_CI(I_A[dd,,], lwr = 0.125, upr = 0.875)
  }
}


# assigning new names
I_A_CI_rho_equal = I_A_CI_rho
I_A_CI_rho_large_equal = I_A_CI_rho_large


#=============================================================================#
# oberved versus predicted forecast
#=============================================================================#
library(viridis)
# library(wesanderson)
library(yarrr)

times_to_plot = seq(1, length(when_to_assimilate), by = 3)
# cols = colorRampPalette(c('navy', 'darkolivegreen3'))(length(times_to_plot))
# cols = wes_palette('GrandBudapest1')
# cols = colorRampPalette(c(cols[1], cols[2], cols[3]))(length(times_to_plot))
# cols = viridis(5)
cols = c(piratepal(palette = "espresso")[5], piratepal(palette = "espresso")[1:4])
pts = c(15:18)

pdf('../output/main_text_figures/supplement/observed_vs_forecasted_ensembles.pdf',
    height = 3, width = 6)

depts_to_plot = c(9, 8, 26, 32)
layout(matrix(1:2, nrow = 1, byrow = T))
par(mar = c(1,1,1,1), oma = c(4,5,1,1))

for(runs_now in c('EM', 'equal')){
  
  plot(-100, -100, xlim = c(0, 4),
       ylim = c(0, 4), xlab = '', ylab = '',
       xaxt = 'n', yaxt = 'n')
  
  
    axis(side = 1, at = c(0:4), labels =c('0','10', '100', '1000', '10000'), las = 1)

  
  if(runs_now %in% 'EM'){
    axis(side = 2, at = c(0:4), labels =c('0','10', '100', '1000', '10000'), las = 1)
    mtext(side = 3, text = 'EM-weighted')
  }else{
    axis(side = 2, at = c(0:4), labels = NA, las = 1)
    mtext(side = 3, text = 'Equally-weighted')
  }
  abline(0, 1, lwd = 2, col = 'grey')
  
  for(dd in depts_to_plot){
    for(tt in times_to_plot){
      
      eval(parse(text = paste0(
        'I_F_CI_rho = I_F_CI_rho_', runs_now
      )))
      
      if(runs_now %in% 5:8){
        
        # non-spatial models
        if(tt <= which(when_to_assimilate == when_to_assimilate_NS[[dd]][2]) - 1){
          
          I_F_to_plot = log10(round(I_F_CI_rho[[1]][dd,,1:time_btwn_assim]))
          I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
          
        }else{
          I_F_to_plot = log10(round(I_F_CI_rho[[tt]][dd,,1:time_btwn_assim]))
          I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
        }
      }else{
        
        # spatial models
        I_F_to_plot = log10(round(I_F_CI_rho[[tt]][dd,,1:time_btwn_assim]))
        I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
      }
      
      dat_current = log10(df_dyn$cases_new[which(
        df_dyn$week %in% (wks_assimilated_at_when_to_assimilate[tt,] + week_first_case - 1) & 
          df_dyn$dept_name == dept_names[dd])])
      dat_current[which(is.infinite(dat_current))] = 0
      
      segments(x0 = dat_current, x1 = dat_current,
               y0 = I_F_to_plot[1,], y1 = I_F_to_plot[2,], 
               col = adjustcolor(cols[which(times_to_plot == tt)], alpha.f = 0.5), lwd = 2)
      
      points(x = dat_current, y = I_F_to_plot[3,], 
             col = adjustcolor(cols[which(times_to_plot == tt)], alpha.f = 0.7), 
             pch = pts[which(depts_to_plot == dd)], cex =1.5)
    }
  }
}
mtext(side = 1, text = 'Observed', line = 2, outer = T)
mtext(side = 2, text = 'Forecasted', line = 3, outer = T)

dev.off()


# legend
pdf('../output/main_text_figures/supplement/observed_vs_forecasted_ensembles_legend.pdf',
    height = 3, width = 1.5)
par(mar = c(1,1,1,1), oma = c(2,0,1,0))
layout(matrix(c(1,2,2,2), nrow = 4))
plot(-100, -100, xlim = c(0, 60), 
     ylim = c(0, 1200), 
     xlab = '', ylab = '', xaxt = 'n', las = 1, yaxt = 'n')
for(tt in times_to_plot){
  rect(xleft = wks_assimilated_at_when_to_assimilate[tt,1],
       xright = wks_assimilated_at_when_to_assimilate[tt,4],
       ybottom = -100,
       ytop = 10000, 
       border = adjustcolor(cols[which(times_to_plot == tt)], alpha.f = 0.4),
       col = adjustcolor(cols[which(times_to_plot == tt)], alpha.f = 0.5))
}

for(dd in depts_to_plot){
  cases_now = df_dept$cases_new[which(df_dept$dept_name == dept_names[dd] & 
                                        df_dept$week %in% week_first_case: max(df_dyn$week))]
  lines(cases_now, lwd = 3, col = rgb(0,0,0,0.4))
}

par(mar = c(1,1,1,6))
color.bar(adjustcolor(c(cols, 'white'), alpha.f = 0.7), min = -1)
mtext('Wks', side = 3, at = 14, line = -1.5)
mtext('of forecast', side = 3, at = 18, line = -3)
axis(side = 4, at = seq(-0.8, 0.8, length.out = 5), 
     las = 1, labels = c('1-4', '13-16', '25-28', '37-40', '49-52'), 
     cex.axis = 1.5)

dev.off()



#=============================================================================#
# oberved versus predicted forecast
#=============================================================================#

Metrics::rmse(obs_EM, preds_EM)
Metrics::rmse(obs_EW, preds_EW)


#=============================================================================#
# oberved versus predicted forecast
#=============================================================================#
library(ggplot2)
library(magrittr)

uncertainty_bounds = dplyr::tibble(EW_uncertainty_bounds = var_EW, 
                                   EM_uncertainty_bounds = var_EM, 
                                   dept = rep(dds_to_plot, each = 4),
                                   time = rep(times_to_plot, each = 16)) 
uncertainty_bounds = tidyr::pivot_longer(
  uncertainty_bounds,
  cols = c(EW_uncertainty_bounds, 
           EM_uncertainty_bounds), 
  names_to = 'var_type',
  values_to = 'var') %>%
  dplyr::mutate(dept = dept_names[dept])

pdf('~/Dropbox/forecasting_tsir/output/main_text_figures/supplement/uncertainty_bounds_example.pdf', 
    height = 3, width = 6)
ggplot(uncertainty_bounds, 
       aes(y = var, x = time, fill = as.factor(dept), 
           color = as.factor(dept))) + 
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~var_type) + 
  xlab('Data assimilation period') + 
  ylab('Magnitude of uncertainty (log10)') + 
  labs(fill = "Department", color = 'Department') + 
  theme_bw()
dev.off()

