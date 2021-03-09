#=============================================================================#
# Author: Rachel Oidtman

# Supplementary Figure
# Observed verus predicted values for nationally aggregated forecast
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
load('../output/nationally_agg_forecasts.RData')

time_btwn_assim = 4
when_to_assimilate = seq(week_first_case + time_btwn_assim, max(df_dyn$week), by = time_btwn_assim)
when_to_assimilate = c(week_first_case, when_to_assimilate)


when_to_assimilate_NS = list()
length(when_to_assimilate_NS) = length(dept_names)

time_btwn_assim = 4


np = dim(I_A_nationally_aggregated)[2]


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
# oberved versus predicted forecast
#=============================================================================#
library(viridis)
library(yarrr)

times_to_plot = seq(1, length(when_to_assimilate), by = 1)
cols = viridis(length(times_to_plot))


pdf('../output/main_text_figures/observed_vs_forecasted_national.pdf',
    height = 6, width = 7)

# pts = c(0:2, 4)
pts = c(15:18)
layout(matrix(1:16, 4,4, byrow = T))
par(mar = c(1,1,1,1), oma = c(4,5,1,1))

for(runs_now in runs[1:16]){
  
  plot(-100, -100, xlim = c(0, 6),
       ylim = c(0, 6), xlab = '', ylab = '',
       xaxt = 'n', yaxt = 'n')
  
  if(runs_now %in% 13:16){
    axis(side = 1, at = c(0:6), labels =c('0','10', '100', '1000', '10000', '1e5', '1e6'), las = 1)
  }else{
    axis(side = 1, at = c(0:6), labels =NA, las = 1)
  }
  
  if(runs_now %in% c(1,5,9,13)){
    axis(side = 2, at = c(0:6), labels =c('0','10', '100', '1000', '10000', '1e5', '1e6'), las = 1)
  }else{
    axis(side = 2, at = c(0:6), labels = NA, las = 1)
  }
  abline(0, 1, lwd = 2, col = 'grey')
  
  for(tt in times_to_plot){
    
    if(runs_now %in% 5:8){
      
      I_F_to_plot = log10(round(rowSums(I_F_CI_rho_national[[tt]][runs_now,,1:time_btwn_assim])))
      I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
      
    }else{
      
      # spatial models
      I_F_to_plot = log10(round(rowSums(I_F_CI_rho_national[[tt]][runs_now,,1:time_btwn_assim])))
      I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
    }
    
    dat_current = log10(sum(df_dyn$cases_new[which(
      df_dyn$week %in% (wks_assimilated_at_when_to_assimilate[tt,] + week_first_case - 1))]))
    
    dat_current[which(is.infinite(dat_current))] = 0
    
    
    segments(x0 = dat_current, x1 = dat_current,
             y0 = I_F_to_plot[1], y1 = I_F_to_plot[2],
             col = adjustcolor(cols[which(times_to_plot == tt)], alpha.f = 0.5), lwd = 2)
    
    points(x = dat_current, y = I_F_to_plot[3], 
           col = adjustcolor(cols[which(times_to_plot == tt)], alpha.f = 0.7), 
           pch = 16, cex =1.5)
  }
}
mtext(side = 1, text = 'Observed', line = 2, outer = T)
mtext(side = 2, text = 'Forecasted', line = 3, outer = T)

dev.off()


# legend
pdf('../output/main_text_figures/observed_vs_forecasted_legend_national.pdf',
    height = 6, width = 1.75)

# color bar
par(mar = c(0,1,1,6.5))
color.bar(adjustcolor(c(cols, 'white'), alpha.f = 0.7), min = -1)
axis(side = 4, at = seq(-0.95, 0.95, length.out = length(times_to_plot)), 
     las = 1, labels = paste0(wks_assimilated_at_when_to_assimilate[,1], '-', wks_assimilated_at_when_to_assimilate[,4]), 
     cex.axis = 1.5)

dev.off()

