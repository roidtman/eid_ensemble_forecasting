#=============================================================================#
# Author: Rachel Oidtman

# Supplementary Figure
# Observed verus predicted values for individual-level models 
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


dept_names_to_plot = c('La Guajira', 'Magdalena', 'Atlántico', 'Cesar', 'Bolívar',
                       'Sucre', 'Córdoba', 'Norte de Santander', 'Antioquia', 
                       'Chocó', 'Santander', 'Arauca', 'Boyacá', 'Vichada', 
                       'Casanare', 'Cudinamarca', 'Caldas', 'Risaralda', 'Tolima', 
                       'Valle del Cauca', 'Meta', 'Bogota', 'Quindío', 'Guainía', 'Huila', 
                       'Cauca', 'Caquetá', 'Guaviare', 'Nariño', 'Vaupés', 
                       'Putomayo', 'Amazonas')

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
# library(wesanderson)

times_to_plot = seq(1, length(when_to_assimilate), by = 3)
# cols = colorRampPalette(c('navy', 'darkolivegreen3'))(length(times_to_plot))
# cols = wes_palette('GrandBudapest1')
# cols = colorRampPalette(c(cols[1], cols[2], cols[3]))(length(times_to_plot))
# cols = viridis(5)
cols = c(piratepal(palette = "espresso")[5], piratepal(palette = "espresso")[1:4])

pdf('../output/main_text_figures/observed_vs_forecasted.pdf',
    height = 6, width = 7)

# pts = c(0:2, 4)
pts = c(15:18)
depts_to_plot = c(9, 8, 26, 32)
layout(matrix(1:16, 4,4, byrow = T))
par(mar = c(1,1,1,1), oma = c(4,5,1,1))

for(runs_now in runs[1:16]){
  
  plot(-100, -100, xlim = c(0, 4),
       ylim = c(0, 4), xlab = '', ylab = '',
       xaxt = 'n', yaxt = 'n')
  
  if(runs_now %in% 13:16){
    axis(side = 1, at = c(0:4), labels =c('0','10', '100', '1000', '10000'), las = 1)
  }else{
    axis(side = 1, at = c(0:4), labels =NA, las = 1)
  }
  
  if(runs_now %in% c(1,5,9,13)){
    axis(side = 2, at = c(0:4), labels =c('0','10', '100', '1000', '10000'), las = 1)
  }else{
    axis(side = 2, at = c(0:4), labels = NA, las = 1)
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
          
          I_F_to_plot = log10(round(rowSums(I_F_CI_rho[[1]][dd,,1:time_btwn_assim])))
          # I_F_to_plot = log10(round(I_F_CI_rho[[1]][dd,,1:time_btwn_assim]))
          I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
          
        }else{
          I_F_to_plot = log10(round(rowSums(I_F_CI_rho[[tt]][dd,,1:time_btwn_assim])))
          # I_F_to_plot = log10(round(I_F_CI_rho[[tt]][dd,,1:time_btwn_assim]))
          I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
        }
      }else{
        
        # spatial models
        I_F_to_plot = log10(round(rowSums(I_F_CI_rho[[tt]][dd,,1:time_btwn_assim])))
        # I_F_to_plot = log10(round(I_F_CI_rho[[tt]][dd,,1:time_btwn_assim]))
        I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
      }
      
      dat_current = log10(sum(df_dyn$cases_new[which(
        df_dyn$week %in% (wks_assimilated_at_when_to_assimilate[tt,] + week_first_case - 1) & 
          df_dyn$dept_name == dept_names[dd])]))
    
      dat_current[which(is.infinite(dat_current))] = 0
      
      
      segments(x0 = dat_current, x1 = dat_current,
               y0 = I_F_to_plot[1], y1 = I_F_to_plot[2],
               col = adjustcolor(cols[which(times_to_plot == tt)], alpha.f = 0.5), lwd = 2)
      
      points(x = dat_current, y = I_F_to_plot[3], 
             col = adjustcolor(cols[which(times_to_plot == tt)], alpha.f = 0.7), 
             pch = pts[which(depts_to_plot == dd)], cex =1.5)
    }
  }
}
mtext(side = 1, text = 'Observed', line = 2, outer = T)
mtext(side = 2, text = 'Forecasted', line = 3, outer = T)

dev.off()


# legend
pdf('../output/main_text_figures/observed_vs_forecasted_legend.pdf',
    height = 6, width = 1.75)
par(mar = c(1,1,1,1), oma = c(2,0,1,0))
layout(matrix(c(1,2, 3,3,3), nrow = 5))

# legend
plot(-100, -100, xlim = c(0, 1), ylim = c(0, 1), xaxt = 'n', yaxt = 'n', bty = 'n',
     xlab = '', ylab = '')
legend('top', pch = pts, legend = c('Antioquia', 'Norte de Santander', 'Cauca', 'Amazonas'), pt.cex = 1.5,
       cex = 1.25, bty = 'n')

# time series
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

# color bar
par(mar = c(0,1,1,6.5))
color.bar(adjustcolor(c(cols, 'white'), alpha.f = 0.7), min = -1)
mtext('Wks', side = 3, at = 13, line = -1)
mtext('of forecast', side = 3, at = 16, line = -2.5)
axis(side = 4, at = seq(-0.8, 0.8, length.out = 5), 
     las = 1, labels = c('1-4', '13-16', '25-28', '37-40', '49-52'), 
     cex.axis = 1.5)

dev.off()


