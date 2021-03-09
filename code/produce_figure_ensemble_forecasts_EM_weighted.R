#=============================================================================#
# Author: Rachel Oidtman

# Supplementary Figure 
# Ensemble-weighted forecasts on the department level 
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


#=============================================================================#
# produce forecast through time with weights
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

I_A_CI_rho = array(NA, dim = c(length(dept_names), 3,  length(week_first_case:max(df_dyn$week))))
I_A_CI_rho_large = array(NA, dim = c(length(dept_names), 3, length(week_first_case:max(df_dyn$week))))
for(dd in 1:length(dept_names)){
  
  if(dd != 22){
    I_A_CI_rho[dd,,] = write_CI(I_A[dd,,])
    I_A_CI_rho_large[dd,,] = write_CI(I_A[dd,,], lwr = 0.125, upr = 0.875)
  }
}


#=============================================================================#
# produce forecast through time with weights
#=============================================================================#

dept_names_to_plot = c('La Guajira', 'Magdalena', 'Atlántico', 'Cesar', 'Bolívar',
                       'Sucre', 'Córdoba', 'Norte de Santander', 'Antioquia', 
                       'Chocó', 'Santander', 'Arauca', 'Boyacá', 'Vichada', 
                       'Casanare', 'Cudinamarca', 'Caldas', 'Risaralda', 'Tolima', 
                       'Valle del Cauca', 'Meta', 'Bogota', 'Quindío', 'Guainía', 'Huila', 
                       'Cauca', 'Caquetá', 'Guaviare', 'Nariño', 'Vaupés', 
                       'Putumayo', 'Amazonas')

times_to_plot = seq(1, length(when_to_assimilate), by = 3)

cols_to_plot = c('darkolivegreen1', 'darkolivegreen3', 'darkolivegreen4')
for(tt in times_to_plot){
  f_out = paste0('../output/main_text_figures/supplement/ensemble_forecasts_EM_wts_',
                 when_to_assimilate[tt] - week_first_case,'_wks.pdf')
  pdf(f_out, height = 8, width = 7)
  
  depts_to_plot = c(1:21, 23:32)
  layout(matrix(1:32, 8,4))
  par(mar = c(1,1,1,0.5), oma = c(3,4.5,3,0.5))
  
  for(dd in depts_to_plot){
    
    I_A_to_plot = log10(round(I_A_CI_rho[dd,,]))
    I_A_to_plot[which(is.infinite(I_A_to_plot))] = 0
    I_A_to_plot[which(is.na(I_A_to_plot))] = 0
    
    I_A_to_plot_large = log10(round(I_A_CI_rho_large[dd,,]))
    I_A_to_plot_large[which(is.infinite(I_A_to_plot_large))] = 0
    I_A_to_plot_large[which(is.na(I_A_to_plot_large))] = 0
    
    plot(-100, -100, xlim = c(0, 60), 
         ylim = c(0, 4.5), 
         xlab = '', ylab = '', xaxt = 'n', las = 1, yaxt = 'n') 
    
    cases_now = log10(df_dept$cases_new[which(df_dept$dept_name == dept_names[dd] & 
                                                df_dept$week %in% week_first_case: max(df_dyn$week))])
    lines(cases_now, cex = 1.2, pch = 1, type = 'h', lwd = 2, 
          col = adjustcolor('navy', alpha.f = 0.8))
    
    I_F_to_plot = log10(round(I_F_CI_rho[[tt]][dd,,]))
    I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
    
    I_F_to_plot_large = log10(round(I_F_CI_rho_large[[tt]][dd,,]))
    I_F_to_plot_large[which(is.infinite(I_F_to_plot_large))] = 0
    
    
    abline(v = when_to_assimilate[tt] - week_first_case, 
           col = adjustcolor('deeppink3', alpha.f = 0.8), lwd = 3)
    
    if(tt == 1){
      mtext(side = 3, text = dept_names_to_plot[dd], cex = 0.75)
    }
    
    
    if(tt > 1){
      # 75% CI I_A
      polygon(x = c(1:(when_to_assimilate[tt] - week_first_case),
                    rev(1:(when_to_assimilate[tt] - week_first_case))),
              y = c(I_A_to_plot_large[1, 1:(when_to_assimilate[tt] - week_first_case)],
                    rev(I_A_to_plot_large[2, 1:(when_to_assimilate[tt] - week_first_case)])),
              col = adjustcolor(cols_to_plot[2], alpha.f = 0.3), 
              border = adjustcolor(cols_to_plot[2], alpha.f = 0.3))
      
      # 50% CI I_A
      polygon(x = c(1:(when_to_assimilate[tt] - week_first_case),
                    rev(1:(when_to_assimilate[tt] - week_first_case))),
              y = c(I_A_to_plot[1, 1:(when_to_assimilate[tt] - week_first_case)],
                    rev(I_A_to_plot[2, 1:(when_to_assimilate[tt] - week_first_case)])),
              col = adjustcolor(cols_to_plot[2], alpha.f = 0.3), 
              border = adjustcolor(cols_to_plot[2], alpha.f = 0.3))
    }
    
    # 75% CI I_F
    polygon(x = c((when_to_assimilate[tt] - week_first_case) : 
                    (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)),
                  rev((when_to_assimilate[tt] - week_first_case) : 
                        (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)))),
            y = c(I_F_to_plot_large[1,],
                  rev(I_F_to_plot_large[2,])), 
            col = adjustcolor(cols_to_plot[1], alpha.f = 0.1), 
            border = adjustcolor(cols_to_plot[1], alpha.f = 0.1))
    
    # 50% CI I_F
    polygon(x = c((when_to_assimilate[tt] - week_first_case) : 
                    (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)),
                  rev((when_to_assimilate[tt] - week_first_case) : 
                        (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)))),
            y = c(I_F_to_plot[1,],
                  rev(I_F_to_plot[2,])), 
            col = adjustcolor(cols_to_plot[1], alpha.f = 0.2), 
            border = adjustcolor(cols_to_plot[1], alpha.f = 0.2))
    
    lines(I_A_to_plot[3, 1:(when_to_assimilate[tt] - week_first_case)], col = cols_to_plot[2], lwd = 4)
    lines(x = (when_to_assimilate[tt] - week_first_case) : 
            (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)), 
          y = I_F_to_plot[3,], col = cols_to_plot[2], lwd = 4)
    
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
