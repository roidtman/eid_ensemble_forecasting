#=============================================================================#
# Author: Rachel Oidtman

# Supplementary Figure 
# Comparing forecasts from dynamic and static R models with Rt and R values 
#=============================================================================#


#=============================================================================#
# load in processed epi data
#=============================================================================#

load('../data/processed/processed_epi_data.RData')
load('../data/processed/temp_mosq_data_frame.RData')

load('../output/main_text_figures/R_output_for_RvsRt_figure.RData')


# assign file path
runs = c(2,4)
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
# figuring out months in time period
#=============================================================================#

library(lubridate)
date_first_case = ymd('2015-08-09')

dates_in_time_frame = date_first_case + (seq(7, (59 * 7), by = 7))
months_in_time_frame = month(dates_in_time_frame)
singular_months_in_time_frame = c(8:12, 1:9)


#=============================================================================#
# making main figure
#=============================================================================#


pdf(file = '../output/main_text_figures/static_vs_dynamic_R.pdf', 
    height = 4.5, width = 7.5)

depts_to_plot = c(18, 7)

# define peak week
pks_to_plot = sapply(depts_to_plot, function(ff) which.max(
  df_dept$cases_new[which(df_dept$dept_name == dept_names[ff])][week_first_case: max(df_dyn$week)]))

pp_c = 5

R0_by_month_by_dept_CI = array(NA, dim = c(3, length(singular_months_in_time_frame), 32))
R0_static_by_dept_CI = array(NA, dim = c(3, length(singular_months_in_time_frame), 32))

for(mm in 1:length(singular_months_in_time_frame)){
  print(mm)
  for(dd in depts_to_plot){
    R0_by_month_by_dept_CI[,mm,dd] = quantile(
      R0_by_month_by_dept[[5]][,singular_months_in_time_frame[mm],dd],
      probs = c(0.025, 0.975, 0.5))
  }
}

for(mm in 1:length(singular_months_in_time_frame)){
  for(dd in depts_to_plot){
    R0_static_by_dept_CI[,mm,dd] = quantile(R0_static_by_dept[[5]][,dd],
                                         probs = c(0.025, 0.975, 0.5))
  }
}

R0_by_month_by_dept_samples = array(NA, dim = c(np, length(singular_months_in_time_frame), 32))
R0_static_by_dept_samples = array(NA, dim = c(np, length(singular_months_in_time_frame), 32))

for(mm in 1:length(singular_months_in_time_frame)){
  print(mm)
  for(dd in depts_to_plot){
    R0_by_month_by_dept_samples[,mm,dd] = R0_by_month_by_dept[[5]][ ,singular_months_in_time_frame[mm], dd]
  }
}

for(mm in 1:length(singular_months_in_time_frame)){
  for(dd in depts_to_plot){
    R0_static_by_dept_samples[,mm,dd] = R0_static_by_dept[[5]][,dd]
  }
}


cols_to_plot_static = c('slategray1', 'slategray3', 'slategray4', 'navy')
cols_to_plot_dynamic = c('indianred1', 'indianred3', 'indianred4', 'firebrick4')

layout(matrix(c(1,1,1,1,5,5,5,6,6,6,
                2,2,2,2,5,5,5,6,6,6,
                3,3,3,3,7,7,7,8,8,8,
                4,4,4,4,7,7,7,8,8,8), byrow = T, nrow = 4))


par(mar = c(1,2,0,1), oma = c(3, 3, 1, 0))
for(dd in depts_to_plot){
  
  par(mar = c(0.5,2,2,2.5))
  plot(-100, -100, xlim = c(1, length(singular_months_in_time_frame)), 
       ylim = c(0, max(R0_by_month_by_dept_CI[,,dd])), las = 1, xaxt = 'n')
  axis(side = 1, at = 1:length(singular_months_in_time_frame), 
       labels = NA)
  mtext(side = 2, text = 'R', line = 1, cex = 0.8, outer = T)
  mtext(side = 3, text = dept_names[dd], cex = 0.75)
  
  for(pp in 1:np){
    lines(R0_by_month_by_dept_samples[pp,,dd], 
          col = adjustcolor(cols_to_plot_dynamic[1], alpha.f = 0.05))
  }
  lines(R0_by_month_by_dept_CI[3,,dd], col = cols_to_plot_dynamic[3], lwd = 3)
  abline(h = 1, lwd = 2, col = 'black')
  
  
  par(mar = c(2,2,0.5,2.5))
  plot(-100, -100, xlim = c(1, length(singular_months_in_time_frame)), 
       ylim = c(0, max(R0_by_month_by_dept_CI[,,dd])), las = 1, xaxt = 'n')
  axis(side = 1, at = 1:length(singular_months_in_time_frame), 
       labels = month.abb[singular_months_in_time_frame])
  # mtext(side = 2, text = 'R', line = 2.5, cex = 0.8)
  
  for(pp in 1:np){
    lines(R0_static_by_dept_samples[pp,,dd], 
          col = adjustcolor(cols_to_plot_static[1], alpha.f = 0.05))
  }
  lines(R0_static_by_dept_CI[3,,dd], col = cols_to_plot_static[3], lwd = 3)
  abline(h = 1, lwd = 2, col = 'black')
  
  lines(R0_static_by_dept_CI[3,,dd], col = cols_to_plot_static[2], lwd = 3)
}
mtext(side = 1, text = 'Months (ordered by epidemic timeline)', line = 2.5, cex = 0.8)

mtext(side = 2, text = 'New weekly cases (log10)', line = -20.5, cex = 0.8, outer = T)


# par(mar = c(1.5,2,2,1), oma = c(3, 3, 1, 0))
times_to_plot = c(7, 7)
for(dd in depts_to_plot){
  
  for(rr in 1:length(runs)){
    if(rr == 1){
      par(mar = c(1.5,2.5,2,0))
    }else{
      par(mar = c(1.5,0.5,2,1.5))
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
      
      if(tt %in% c(times_to_plot[which(depts_to_plot == dd)])){
        I_F_to_plot = log10(round(I_F_CI_rho[[tt]][dd,,]))
        I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
        
        I_F_to_plot_large = log10(round(I_F_CI_rho_large[[tt]][dd,,]))
        I_F_to_plot_large[which(is.infinite(I_F_to_plot_large))] = 0
        
        plot(-100, -100, xlim = c(0, dim(I_A_to_plot)[2]), 
             ylim = c(0, 4.5),
             xlab = '', ylab = '', xaxt = 'n', las = 1, yaxt= 'n')
        
        if(rr == 1){
          axis(side = 2, las = 1, at = 1:4, labels= c(10, 100, 1000, 10000))
          mtext(side = 3, text = dept_names[dd], at = 65, cex = 0.75)
          
          if(dd == depts_to_plot[1]){
            mtext(side = 3, line = 1.25, text = 'Static R', col = 'slategray4')
          }
        }
        
        if(rr == 2){
          if(dd == depts_to_plot[1]){
            mtext(side = 3, line = 1.25, text = 'Dynamic R', col = 'indianred4')
          }
        }
        
        
        abline(v = when_to_assimilate[tt] - week_first_case, col = 'navy', lwd = 2)
        
        if(tt == 1){
          mtext(side = 3, text = dept_names[dd], cex = 0.75)
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
                  col = adjustcolor(cols[2], alpha.f = 0.6), 
                  border = adjustcolor(cols[2], alpha.f = 0.6))
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
        
        lines(I_A_to_plot[3, 1:(when_to_assimilate[tt] - week_first_case)], col = cols[3], lwd = 4)
        lines(x = (when_to_assimilate[tt] - week_first_case) : 
                (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)), 
              y = I_F_to_plot[3,], col = cols[2], lwd = 4)
        
        points(cases_now, cex = 1, pch = 4, col = rgb(0,0,0,0.8))
      }
    }
    axis(side = 1)
    
  }
}

mtext(side = 1, text = 'Weeks since 1st reported case in Colombia', line = 2.5, cex = 0.8, at = -5)


dev.off()




#=============================================================================#
# making main figure
# 2 time points 
#=============================================================================#

dept_names_to_plot = c('La Guajira', 'Magdalena', 'Atlántico', 'Cesar', 'Bolívar',
                       'Sucre', 'Córdoba', 'Norte de Santander', 'Antioquia', 
                       'Chocó', 'Santander', 'Arauca', 'Boyacá', 'Vichada', 
                       'Casanare', 'Cudinamarca', 'Caldas', 'Risaralda', 'Tolima', 
                       'Valle del Cauca', 'Meta', 'Quindío', 'Guainía', 'Huila', 
                       'Cauca', 'Caquetá', 'Guaviare', 'Nariño', 'Vaupés', 
                       'Putomayo', 'Amazonas')

pdf(file = '../output/main_text_figures/static_vs_dynamic_R_2_time_points.pdf', 
    height = 4.5, width = 7.5)

layout(matrix(c(1,1,1,1,5,5,5,6,6,6,
                1,1,1,1,5,5,5,6,6,6,
                2,2,2,2,7,7,7,8,8,8,
                2,2,2,2,7,7,7,8,8,8,
                3,3,3,3,9,9,9,10,10,10,
                3,3,3,3,9,9,9,10,10,10,
                4,4,4,4,11,11,11,12,12,12,
                4,4,4,4,11,11,11,12,12,12), byrow = T, nrow = 8))


par(mar = c(1,2,0,1), oma = c(3, 3, 1, 0))
for(dd in depts_to_plot){
  
  par(mar = c(0.5,2,2,2.5))
  plot(-100, -100, xlim = c(1, length(singular_months_in_time_frame)), 
       ylim = c(0, max(R0_by_month_by_dept_CI[,,dd])), las = 1, xaxt = 'n')
  axis(side = 1, at = 1:length(singular_months_in_time_frame), 
       labels = NA)
  mtext(side = 2, text = 'R', line = 1, cex = 0.8, outer = T)
  mtext(side = 3, text = dept_names_to_plot[dd], cex = 0.75)
  
  for(pp in 1:np){
    lines(R0_by_month_by_dept_samples[pp,,dd], 
          col = adjustcolor(cols_to_plot_dynamic[1], alpha.f = 0.05))
  }
  lines(R0_by_month_by_dept_CI[3,,dd], col = cols_to_plot_dynamic[3], lwd = 3)
  abline(h = 1, lwd = 2, col = 'black')
  
  
  par(mar = c(2,2,0.5,2.5))
  plot(-100, -100, xlim = c(1, length(singular_months_in_time_frame)), 
       ylim = c(0, max(R0_by_month_by_dept_CI[,,dd])), las = 1, xaxt = 'n')
  axis(side = 1, at = 1:length(singular_months_in_time_frame), 
       labels = month.abb[singular_months_in_time_frame])
  # mtext(side = 2, text = 'R', line = 2.5, cex = 0.8)
  
  for(pp in 1:np){
    lines(R0_static_by_dept_samples[pp,,dd], 
          col = adjustcolor(cols_to_plot_static[1], alpha.f = 0.05))
  }
  lines(R0_static_by_dept_CI[3,,dd], col = cols_to_plot_static[3], lwd = 3)
  abline(h = 1, lwd = 2, col = 'black')
  
  lines(R0_static_by_dept_CI[3,,dd], col = cols_to_plot_static[2], lwd = 3)
}
mtext(side = 1, text = 'Months (ordered by epidemic timeline)', line = 2.5, cex = 0.8)

mtext(side = 2, text = 'New weekly cases (log10)', line = -20.5, cex = 0.8, outer = T)


# par(mar = c(1.5,2,2,1), oma = c(3, 3, 1, 0))
times_to_plot = c(7, 7)
times_to_plot_2 = c(12,12)
runs = c(4, 2)
for(dd in depts_to_plot){
  
  for(rr in 1:length(runs)){
    
    
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
      
      if(tt == times_to_plot[1] & rr == 1){
        par(mar = c(0.5,2.5,2,0))
      }else if(tt == times_to_plot[1] & rr == 2){
        par(mar = c(2,2.5,0.5,0))
      }else if(tt == times_to_plot_2[1] & rr == 2){
        par(mar = c(2,0.5,0.5,1.5))
      }else{
        par(mar = c(0.5,0.5,2,1.5))
      }
      
      # if(tt == times_to_plot[2] & rr == 1){
      #   par(mar = c(0.5,0.5,2,0))
      # }
      
      
      if(tt %in% c(times_to_plot[which(depts_to_plot == dd)], times_to_plot_2[which(depts_to_plot == dd)])){
        I_F_to_plot = log10(round(I_F_CI_rho[[tt]][dd,,]))
        I_F_to_plot[which(is.infinite(I_F_to_plot))] = 0
        
        I_F_to_plot_large = log10(round(I_F_CI_rho_large[[tt]][dd,,]))
        I_F_to_plot_large[which(is.infinite(I_F_to_plot_large))] = 0
        
        plot(-100, -100, xlim = c(0, dim(I_A_to_plot)[2]), 
             ylim = c(0, 4.5),
             xlab = '', ylab = '', xaxt = 'n', las = 1, yaxt= 'n')
        lines(cases_now, cex = 1.2, pch = 1, type = 'h', lwd = 2, 
              col = adjustcolor('navy', alpha.f = 0.4))
        
        if(tt == times_to_plot[which(depts_to_plot == dd)]){
          axis(side = 2, las = 1, at = 1:4, labels= c(10, 100, 1000, 10000))
          # mtext(side = 3, text = dept_names[dd], at = 65, cex = 0.75)
        }
        
        # abline(v = when_to_assimilate[tt] - week_first_case, col = 'navy', lwd = 2)
        rect(xleft = when_to_assimilate[tt] - week_first_case,
             xright = when_to_assimilate[tt + 1] - week_first_case,
             ybottom = -100, ytop = 100, col = rgb(0,0,0,0.2), border = NA)
      
        if(rr == 1 & tt == times_to_plot[1]){
          mtext(side = 3, text = dept_names_to_plot[dd], cex = 0.75, at = 65)
        }
        
        if(rr == 2){
          axis(side = 1)
        }else{
          axis(side = 1, labels = NA)
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
                  col = adjustcolor(cols[2], alpha.f = 0.6), 
                  border = adjustcolor(cols[2], alpha.f = 0.6))
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
        
        lines(I_A_to_plot[3, 1:(when_to_assimilate[tt] - week_first_case)], col = cols[3], lwd = 4)
        lines(x = (when_to_assimilate[tt] - week_first_case) : 
                (ncol(I_F_to_plot) + (when_to_assimilate[tt] - week_first_case - 1)), 
              y = I_F_to_plot[3,], col = cols[3], lwd = 4)
        
        # points(cases_now, cex = 1, pch = 4, col = rgb(0,0,0,0.8))
      }
    }
  }
}

mtext(side = 1, text = 'Weeks since 1st reported case in Colombia', line = 1, cex = 0.8, at = -5)

dev.off()


