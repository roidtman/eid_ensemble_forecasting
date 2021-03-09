#=============================================================================#
# Author: Rachel Oidtman

# Supplementary Figure 
# Comparing non-spatial to spatial models
#=============================================================================#


#=============================================================================#
# load in processed epi data
#=============================================================================#

load('../data/processed/processed_epi_data.RData')

# assign file path
runs = c(2,5)
file_paths = paste0('no_tm_ini_',runs)

dept_names_to_plot = c('La Guajira', 'Magdalena', 'Atlántico', 'Cesar', 'Bolívar',
                       'Sucre', 'Córdoba', 'Norte de Santander', 'Antioquia', 
                       'Chocó', 'Santander', 'Arauca', 'Boyacá', 'Vichada', 
                       'Casanare', 'Cudinamarca', 'Caldas', 'Risaralda', 'Tolima', 
                       'Valle del Cauca', 'Meta', 'Bogota', 'Quindío', 'Guainía', 'Huila', 
                       'Cauca', 'Caquetá', 'Guaviare', 'Nariño', 'Vaupés', 
                       'Putomayo', 'Amazonas')

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

np = dim(I_A_over_time_rho_sample_2)[2]


wks_assimilated_at_when_to_assimilate = matrix(NA, nrow = length(when_to_assimilate),
                                               ncol = time_btwn_assim)
for(ii in 1:length(when_to_assimilate)){
  wks_assimilated_at_when_to_assimilate[ii, ] = (when_to_assimilate[ii] - week_first_case + 1) : (when_to_assimilate[ii] - week_first_case + time_btwn_assim)
}



#=============================================================================#
# looking at forecast across best fit Rt and R models
#=============================================================================#

f = paste0('../output/main_text_figures/supplement/comparing_forecasts_across_spatial_non-spatial.pdf')
pdf(file = f, height = 3, width = 7.5)

depts_to_plot = c(19, 13,8)
# depts_to_plot = c(30:32)
times_to_plot = seq(1, length(when_to_assimilate), by = 6)

layout(matrix(1:(length(times_to_plot)* length(depts_to_plot) * length(runs)) , 
              nrow = length(times_to_plot), ncol = 6, byrow = F))


cols_to_plot_static = c('darkseagreen2', 'darkseagreen3', 'darkseagreen4')
cols_to_plot_dynamic = c('gold1', 'gold3', 'gold4')


par(mar = c(0.5,3.5,0,1), oma = c(3.5,3,3,0))
for(dd in depts_to_plot){
  par(mar = c(0.5,3.5,0,1))
  
  for(rr in 1:length(runs)){
    if(rr == 1){
      par(mar = c(0.5,2.5,0,0))
    }else{
      par(mar = c(0.5,0.5,0,1.5))
    }
    
    
    if(runs[rr] == 2){
      cols = cols_to_plot_static
    }else{
      cols = cols_to_plot_dynamic
    }
    
    eval(parse(text = paste0(
      'I_A_CI_rho = I_A_CI_rho_', runs[rr]
    )))
    eval(parse(text = paste0(
      'I_A_CI_rho_large = I_A_CI_rho_large_', runs[rr]
    )))
    eval(parse(text = paste0(
      'I_F_CI_rho = I_F_CI_rho_', runs[rr]
    )))
    eval(parse(text = paste0(
      'I_F_CI_rho_large = I_F_CI_rho_large_', runs[rr]
    )))
    
    
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
        lines(cases_now, cex = 1.2, pch = 1, type = 'h', lwd = 2, 
              col = adjustcolor('navy', alpha.f = 0.4))
        
        if(rr == 1){
          axis(side = 2, las = 1, at = 1:4, labels= c(10, 100, 1000, 10000))
        }
        
        
        abline(v = when_to_assimilate[tt] - week_first_case, col = 'navy', lwd = 2)
        
        if(tt == 1){
          mtext(side = 3, text = dept_names_to_plot[dd], cex = 0.6)
          
          if(rr == 1){
            mtext(side = 3, text = 'Spatial', line = 1, col = cols_to_plot_static[3], cex = 0.7)
          }else{
            mtext(side = 3, text = 'Non-spatial', line = 1, col = cols_to_plot_dynamic[2], cex = 0.7)
          }
        }
        
        if(tt > 1){
          # 75% CI I_A
          polygon(x = c(1:(when_to_assimilate[tt] - week_first_case),
                        rev(1:(when_to_assimilate[tt] - week_first_case))),
                  y = c(I_A_to_plot_large[1, 1:(when_to_assimilate[tt] - week_first_case)],
                        rev(I_A_to_plot_large[2, 1:(when_to_assimilate[tt] - week_first_case)])),
                  col = adjustcolor(cols[2], alpha.f = 0.3), 
                  border = adjustcolor(cols[2], alpha.f = 0.3))
          
          # 50% CI I_A
          polygon(x = c(1:(when_to_assimilate[tt] - week_first_case),
                        rev(1:(when_to_assimilate[tt] - week_first_case))),
                  y = c(I_A_to_plot[1, 1:(when_to_assimilate[tt] - week_first_case)],
                        rev(I_A_to_plot[2, 1:(when_to_assimilate[tt] - week_first_case)])),
                  col = adjustcolor(cols[2], alpha.f = 0.3), 
                  border = adjustcolor(cols[2], alpha.f = 0.3))
        }
        
        # 75% CI I_F
        polygon(x = c((when_to_assimilate[tt] - week_first_case) : 
                        (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)),
                      rev((when_to_assimilate[tt] - week_first_case) : 
                            (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)))),
                y = c(I_F_to_plot_large[1,],
                      rev(I_F_to_plot_large[2,])), 
                col = adjustcolor(cols[1], alpha.f = 0.1), 
                border = adjustcolor(cols[1], alpha.f = 0.1))
        
        # 50% CI I_F
        polygon(x = c((when_to_assimilate[tt] - week_first_case) : 
                        (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)),
                      rev((when_to_assimilate[tt] - week_first_case) : 
                            (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)))),
                y = c(I_F_to_plot[1,],
                      rev(I_F_to_plot[2,])), 
                col = adjustcolor(cols[1], alpha.f = 0.3), 
                border = adjustcolor(cols[1], alpha.f = 0.3))
        
        lines(I_A_to_plot[3, 1:(when_to_assimilate[tt] - week_first_case)], 
              col = cols[2], lwd = 3)
        lines(x = (when_to_assimilate[tt] - week_first_case) : 
                (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)), 
              y = I_F_to_plot[3,], col = cols[2], lwd = 3)
        
        # points(cases_now, cex = 1, pch = 4, col = rgb(0,0,0,0.8))
      }
    }
    axis(side = 1)
    
  }
}

mtext(side = 1, line = 2, text = 'Number of weeks since first case reported in Colombia', outer = T)
mtext(side = 2, line = 1, text = 'New weekly cases', outer = T)

dev.off()







