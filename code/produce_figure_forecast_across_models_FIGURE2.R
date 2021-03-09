#=============================================================================#
# Author: Rachel Oidtman

# FIGURE 2
# Median forecasts

# Supplementary Figure 
# Median forecasts on the department-level
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
# MAIN TEXT FIGURE
#=============================================================================#

add_ensemble = T

if(add_ensemble){
  f = '../output/main_text_figures/median_forecasts_with_ensemble_CI.pdf'
  pdf(f, 
      height = 4.5 , width = 7)
}else{
  f = '../output/main_text_figures/median_forecasts.pdf'
  pdf(f, 
      height = 4.5 , width = 7)
}
depts_to_plot = c(9, 8, 26, 32)
# depts_to_plot = c(30:32, 1)
times_to_plot = seq(1, length(when_to_assimilate), by = 3)
cols_to_plot = 'darkolivegreen1'

layout(matrix(1:(length(times_to_plot)* length(depts_to_plot)) , 
              nrow = length(times_to_plot), ncol = 4, byrow = F))
par(mar = c(0.5,1,0,1), oma = c(3,5,2,3))

for(dd in depts_to_plot){
  for(tt in times_to_plot){
    plot(-100, -100, xlim = c(0, 60), 
         ylim = c(0, 3.75), 
         xlab = '', ylab = '', xaxt = 'n', las = 1, yaxt = 'n')
    if(dd == depts_to_plot[1]){
      axis(side = 2, at = c(1:4), labels =c('10', '100', '1000', '10000'), las = 1)
    }
    cases_now = log10(df_dept$cases_new[which(df_dept$dept_name == dept_names[dd] & 
                                                df_dept$week %in% week_first_case: max(df_dyn$week))])
    # lines(cases_now, cex = 1.2, pch = 1, type = 'h', lwd = 2,
    #       col = adjustcolor('navy', alpha.f = 0.4))
    
    if(dd == depts_to_plot[4]){
      axis(side = 4, at = 2, labels = when_to_assimilate[tt] - week_first_case, las = 1,
           lwd.ticks = 0)
      # mtext(side = 4, line = 1, text = when_to_assimilate[tt] - week_first_case)
    }
    
    # add ensemble
    if(add_ensemble){
      I_A_to_plot = log10(round(I_A_CI_rho_ensemble[dd,,]))
      I_A_to_plot[which(is.infinite(I_A_to_plot))] = 0
      I_A_to_plot[which(is.na(I_A_to_plot))] = 0
      
      I_A_to_plot_large = log10(round(I_A_CI_rho_large_ensemble[dd,,]))
      I_A_to_plot_large[which(is.infinite(I_A_to_plot_large))] = 0
      I_A_to_plot_large[which(is.na(I_A_to_plot_large))] = 0
      
      
      I_F_to_plot = log10(round(I_F_CI_rho_ensemble[[tt]][dd,,]))
      I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
      
      I_F_to_plot_large = log10(round(I_F_CI_rho_large_ensemble[[tt]][dd,,]))
      I_F_to_plot_large[which(is.infinite(I_F_to_plot_large))] = 0

      if(tt > 1){
        # 75% CI I_A
        # polygon(x = c(1:(when_to_assimilate[tt] - week_first_case),
        #               rev(1:(when_to_assimilate[tt] - week_first_case))),
        #         y = c(I_A_to_plot_large[1, 1:(when_to_assimilate[tt] - week_first_case)],
        #               rev(I_A_to_plot_large[2, 1:(when_to_assimilate[tt] - week_first_case)])),
        #         col = adjustcolor(cols_to_plot[2], alpha.f = 0.3), 
        #         border = adjustcolor(cols_to_plot[2], alpha.f = 0.3))
        
        # 50% CI I_A
        polygon(x = c(1:(when_to_assimilate[tt] - week_first_case),
                      rev(1:(when_to_assimilate[tt] - week_first_case))),
                y = c(I_A_to_plot[1, 1:(when_to_assimilate[tt] - week_first_case)],
                      rev(I_A_to_plot[2, 1:(when_to_assimilate[tt] - week_first_case)])),
                col = adjustcolor(cols_to_plot[1], alpha.f = 0.5), 
                border = adjustcolor(cols_to_plot[1], alpha.f = 0.5))
      }
      
      # 75% CI I_F
      # polygon(x = c((when_to_assimilate[tt] - week_first_case) : 
      #                 (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)),
      #               rev((when_to_assimilate[tt] - week_first_case) : 
      #                     (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)))),
      #         y = c(I_F_to_plot_large[1,],
      #               rev(I_F_to_plot_large[2,])), 
      #         col = adjustcolor(cols_to_plot[1], alpha.f = 0.1), 
      #         border = adjustcolor(cols_to_plot[1], alpha.f = 0.1))
      
      # 50% CI I_F
      polygon(x = c((when_to_assimilate[tt] - week_first_case) : 
                      (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)),
                    rev((when_to_assimilate[tt] - week_first_case) : 
                          (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)))),
              y = c(I_F_to_plot[1,],
                    rev(I_F_to_plot[2,])), 
              col = adjustcolor(cols_to_plot[1], alpha.f = 0.5), 
              border = adjustcolor(cols_to_plot[1], alpha.f = 0.5))
    }
    
    # add medians for runs
    for(runs_now in runs[1:16]){
      
      eval(parse(text = paste0(
        'I_A_CI_rho = I_A_CI_rho_', runs_now
      )))
      eval(parse(text = paste0(
        'I_F_CI_rho = I_F_CI_rho_', runs_now
      )))
      
      I_A_to_plot = log10(round(I_A_CI_rho[dd,,]))
      I_A_to_plot[which(is.infinite(I_A_to_plot))] = 0
      I_A_to_plot[which(is.na(I_A_to_plot))] = 0
      
      if(runs_now %in% 5:8){
        
        # non-spatial models
        
        if(tt <= which(when_to_assimilate == when_to_assimilate_NS[[dd]][2]) - 1){
          
          I_F_to_plot = log10(round(I_F_CI_rho[[1]][dd,,]))
          I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
          
          lines(x = (1:ncol(I_F_to_plot)), 
                y = I_F_to_plot[3,], col = rgb(0,0,0,0.3), lwd = 1)
        }else{
          I_F_to_plot = log10(round(I_F_CI_rho[[tt]][dd,,]))
          I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
          
          lines(I_A_to_plot[3, 1:(when_to_assimilate[tt] - week_first_case)], 
                col = rgb(0,0,0,0.5), lwd = 1)
          lines(x = (when_to_assimilate[tt] - week_first_case) : 
                  (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)), 
                y = I_F_to_plot[3,], col = rgb(0,0,0,0.5), lwd = 1)
        }
      }else{
        
        # spatial models
        
        I_F_to_plot = log10(round(I_F_CI_rho[[tt]][dd,,]))
        I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
        
        lines(I_A_to_plot[3, 1:(when_to_assimilate[tt] - week_first_case)], 
              col = rgb(0,0,0,0.5), lwd = 1)
        lines(x = (when_to_assimilate[tt] - week_first_case) : 
                (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)), 
              y = I_F_to_plot[3,], col = rgb(0,0,0,0.5), lwd = 1)
      }
    }

    abline(v = when_to_assimilate[tt] - week_first_case, 
           col = adjustcolor('deeppink3', alpha.f = 0.8), lwd = 3)
    
    points(cases_now, cex = 1.25, pch = 16, lwd = 2,
          col = adjustcolor('navy', alpha.f = 0.5))
    
    if(tt == times_to_plot[length(times_to_plot)]){
      axis(side = 1)
    }
    if(tt == times_to_plot[1]){
      mtext(side = 3, text = dept_names_to_plot[dd], cex = 0.75)
    }
  }
}
mtext(side = 1, line = 2, text = 'Number of weeks since first case reported in Colombia', outer = T)
mtext(side = 2, line = 3, text = 'New weekly cases (log10)', outer = T)
mtext(side = 4, line = 2, outer = T, text = 'Week of forecast')

dev.off()


#=============================================================================#
# SUPPLEMENTAL FIGURE
#=============================================================================#

times_to_plot = seq(1, length(when_to_assimilate), by = 3)

for(tt in times_to_plot){
  f_out = paste0('../output/main_text_figures/supplement/median_forecasts_', 
                 when_to_assimilate[tt] - week_first_case,'_wks.pdf')
  pdf(f_out, height = 8, width = 7)

  depts_to_plot = c(1:21, 23:32)
  layout(matrix(1:32, 8,4))
  par(mar = c(1,1,1,0.5), oma = c(3,4.5,3,0.5))
  
  for(dd in depts_to_plot){
    plot(-100, -100, xlim = c(0, 60), 
         ylim = c(0, 3.75), 
         xlab = '', ylab = '', xaxt = 'n', las = 1, yaxt = 'n') 
    
    cases_now = log10(df_dept$cases_new[which(df_dept$dept_name == dept_names[dd] & 
                                                df_dept$week %in% week_first_case: max(df_dyn$week))])
    lines(cases_now, cex = 1.2, pch = 1, type = 'h', lwd = 2, 
          col = adjustcolor('navy', alpha.f = 0.4))
    
    for(runs_now in runs[1:16]){
      
      eval(parse(text = paste0(
        'I_A_CI_rho = I_A_CI_rho_', runs_now
      )))
      eval(parse(text = paste0(
        'I_F_CI_rho = I_F_CI_rho_', runs_now
      )))
      
      I_A_to_plot = log10(round(I_A_CI_rho[dd,,]))
      I_A_to_plot[which(is.infinite(I_A_to_plot))] = 0
      I_A_to_plot[which(is.na(I_A_to_plot))] = 0
      
      if(runs_now %in% 5:8){
        
        # non-spatial models
        
        if(tt <= which(when_to_assimilate == when_to_assimilate_NS[[dd]][2]) - 1){
          
          I_F_to_plot = log10(round(I_F_CI_rho[[1]][dd,,]))
          I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
          
          lines(x = (1:ncol(I_F_to_plot)), 
                y = I_F_to_plot[3,], col = rgb(0,0,0,0.3), lwd = 1)
        }else{
          I_F_to_plot = log10(round(I_F_CI_rho[[tt]][dd,,]))
          I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
          
          lines(I_A_to_plot[3, 1:(when_to_assimilate[tt] - week_first_case)], 
                col = rgb(0,0,0,0.5), lwd = 1)
          lines(x = (when_to_assimilate[tt] - week_first_case) : 
                  (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)), 
                y = I_F_to_plot[3,], col = rgb(0,0,0,0.5), lwd = 1)
        }
      }else{
        
        # spatial models
        
        I_F_to_plot = log10(round(I_F_CI_rho[[tt]][dd,,]))
        I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
        
        lines(I_A_to_plot[3, 1:(when_to_assimilate[tt] - week_first_case)], 
              col = rgb(0,0,0,0.5), lwd = 1)
        lines(x = (when_to_assimilate[tt] - week_first_case) : 
                (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)), 
              y = I_F_to_plot[3,], col = rgb(0,0,0,0.5), lwd = 1)
      }
    }
    
    abline(v = when_to_assimilate[tt] - week_first_case, 
           col = adjustcolor('deeppink3', alpha.f = 0.8), lwd = 3)
    if(dd %in% c(8, 16, 25, 32)){
      axis(side = 1)
    }else{
      axis(side = 1, labels = NA)
    }
    
    if(dd %in% c(1:8)){
      axis(side = 2, at = c(1:4), labels =c('10', '100', '1000', '10000'), las = 1)
    }
    
    mtext(side = 3, text = dept_names_to_plot[dd], cex = 0.75)
  }
  mtext(side = 1, text = 'Weeks since 1st reported case in Colombia', 
        outer = T, line = 1.5)
  mtext(side = 2, text = 'New weekly cases (log10)', outer = T, line = 3)
  mtext(side = 3, text = paste0(when_to_assimilate[tt] - week_first_case, ' weeks of data assimilated into model'), 
        outer = T, line = 1)
  dev.off()
}


