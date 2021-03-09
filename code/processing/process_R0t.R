#=============================================================================#
# Author: Rachel Oidtman

# processing Rt
#=============================================================================#


#=============================================================================#
# load in processed epi data
#=============================================================================#

library(viridis)
library(lubridate)
library(dplyr)
library(sn)


load('../data/processed/processed_epi_data.RData')
load('../data/processed/temp_mosq_data_frame.RData')
date_first_case = ymd('2015-08-09')

# runs
runs = c(2, 4)

# load in temp and mosqsuito occurrence probabiity data
temp_mean = read.csv('../data/dept_tmean_weekly_weighted.csv', stringsAsFactors=F)
aedes_op = read.csv('../data/aegypti_pop_monthly_col_dept.csv', stringsAsFactors= F)


time_btwn_assim = 4
when_to_assimilate = seq(week_first_case + time_btwn_assim, max(df_dyn$week), by = time_btwn_assim)
when_to_assimilate = c(week_first_case, when_to_assimilate)


# when_to_assimilate_NS = list()
# length(when_to_assimilate_NS) = length(dept_names)
# 
# time_btwn_assim = 4
# 
# for(dd_curr in 1:length(dept_names)){
#   if(dd_curr != 22){
#     
#     df_stat_tmp = df_stat[which(df_stat$dept_name == as.character(dept_names[dd_curr])),]
#     df_dyn_tmp = df_dyn[which(df_dyn$dept_name == as.character(dept_names[dd_curr])),]
#     week_first_case_tmp = df_stat_tmp$wk_first_case[which(df_stat_tmp$dept_name == as.character(dept_names[dd_curr]))]
#     
#     when_to_assimilate_tmp = seq(week_first_case_tmp + time_btwn_assim, max(df_dyn$week), by = time_btwn_assim)
#     when_to_assimilate_tmp = c(week_first_case_tmp, when_to_assimilate_tmp)
#     
#     when_to_assimilate_dept_specific = when_to_assimilate[which((when_to_assimilate - week_first_case_tmp) > 0)]
#     when_to_assimilate_dept_specific = c(week_first_case_tmp, when_to_assimilate_dept_specific)
#     
#     when_to_assimilate_NS[[dd_curr]] = when_to_assimilate_dept_specific
#   }
# }
# rm(df_stat_tmp, df_dyn_tmp, week_first_case_tmp)

#=============================================================================#
# load in particles
#=============================================================================#

for(rr in 1:length(runs)){
  
  if(runs[rr] %in% c(2, 4)){
    for(tt in when_to_assimilate){
      
      if(tt == when_to_assimilate[1]){
        f = paste0('../output/no_tm_ini_', runs[rr], '/particles/particles_current_original.csv')
      }else{
        f = paste0('../output/no_tm_ini_', runs[rr], '/particles/particles_current_', tt, '.csv')
      }
      particles = read.csv(f)
      particles = particles[sample(1:nrow(particles), 2000), ]
      
      name_tmp = paste0('particles_', tt, '_', runs[rr])
      assign(name_tmp, particles)
      rm(particles)
      
    }
  }else{
    for(dd in 1:length(dept_names)){
      for(tt in when_to_assimilate_NS[[dd]]){
        
        if(tt == when_to_assimilate_NS[[dd]][1]){
          f = paste0('../output/no_tm_ini_', runs[rr], '/', dd, '/particles/particles_current_original.csv')
        }else{
          f = paste0('../output/no_tm_ini_', runs[rr], '/', dd, '/particles/particles_current_', tt, '.csv')
        }
        
        particles = read.csv(f)
        particles = particles[sample(1:nrow(particles), 2000), ]
        
        name_tmp = paste0('particles_', dd, '_', tt, '_', runs[rr])
        assign(name_tmp, particles)
        rm(particles)
      }
    }
  }
}

# 
np = dim(particles_100_2)[1]

#=============================================================================#
# necessary functions
#=============================================================================#

# calc R0
calc_R0t = function(temp_in, mosq_abund_in, k_in, sn1_in, sn2_in, sn3_in){
  R0_tmp = k_in * dsn(temp_in, sn1_in, sn2_in, sn3_in) * mosq_abund_in
  return(R0_tmp)
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


# get sunday of each week
# for labeling
yearSunday <- function(year) {
  dates <- as.Date(paste(year, "-01-01", sep = "")) + 0:6
  days <- weekdays(dates) == "Sunday"
  seq(dates[days], as.Date(paste(year, "-12-31", sep = "")),  by = "+7 day")
}


#=============================================================================#
# R0 curves
#=============================================================================#

R0_curves_through_time = array(NA, dim = c(length(when_to_assimilate), np, length(seq(15, 32, by = 0.1))))

for(tt in 1:length(when_to_assimilate)){
  eval(parse(text = paste0(
    'particles = particles_', when_to_assimilate[tt], '_4'
  )))
  
  for(pp in 1:np){
    R0_curves_through_time[tt, pp, ] = particles$k_scale[pp] * 
      calc_R0t(temp_in = seq(15, 32, by = 0.1), 
               mosq_abund_in = mean(
                 df_dept_temp_mosq$mosq_abund[which(df_dept_temp_mosq$week %in% week_first_case : max(df_dyn$week))]
                 ), 
               k_in = particles$k[pp],
               sn1_in = particles$sn1[pp], 
               sn2 = particles$sn2[pp], 
               sn3 = particles$sn3[pp])
  }
}


# 95% CrI
R0_curves_through_time_CrI = array(NA, dim = c(length(when_to_assimilate), 3, 
                                               length(seq(15, 32, by = 0.1))))
for(tt in 1:length(when_to_assimilate)){
  R0_curves_through_time_CrI[tt, , ] = write_CI(R0_curves_through_time[tt,,])
}

# pdf('../output/no_tm_ini_4/R_versus_temperature.pdf', height = 4.5, width = 7)
pdf('../output/main_text_figures/R_versus_temperature.pdf', height = 4.5, width = 7)
par(mar = c(2,2,2,1), oma = c(2,2,0,1))
# cols = viridis::plasma(length(when_to_assimilate))
layout(matrix(1:15, nrow = 3, byrow = T))
for(tt in 1:length(when_to_assimilate)){
  
  plot(-100, -100, xlim = c(1, length(seq(15, 32, by = 0.1))), 
       ylim = c(0, 6), las = 1, xlab = '', ylab = '', xaxs = 'i', xaxt = 'n')
  axis(side = 1, at = seq(1, length(seq(15, 32, by = 0.1)), length.out = 5), 
       labels = round(seq(15, 32, length.out = 5)))
  
  abline(v = which(seq(15, 32, by = 0.1) == 17.5), col = rgb(0,0,0,0.3), lwd = 0.5)
  abline(v = which(seq(15, 32, by = 0.1) == 22.5), col = rgb(0,0,0,0.3), lwd = 0.5)
  abline(v = which(seq(15, 32, by = 0.1) == 27.5), col = rgb(0,0,0,0.3), lwd = 0.5)
  
  abline(v = which(seq(15, 32, by = 0.1) == 20), col = rgb(0,0,0,0.5), lwd = 0.75)
  abline(v = which(seq(15, 32, by = 0.1) == 25), col = rgb(0,0,0,0.5), lwd = 0.75)
  abline(v = which(seq(15, 32, by = 0.1) == 30), col = rgb(0,0,0,0.5), lwd = 0.75)
  
  abline(h = c(2, 3, 4, 5, 6), col = rgb(0,0,0,0.3), lwd = 0.5)
  abline(h = 1, lwd = 3, col= 'indianred1')
  
  polygon(x = c(1:length(seq(15, 32, by = 0.1)), rev(1:length(seq(15, 32, by = 0.1)))),
          y = c(R0_curves_through_time_CrI[tt,1,],
                rev(R0_curves_through_time_CrI[tt,2,])),
          col = adjustcolor('navy', alpha.f = 0.2),
          border = adjustcolor('navy', alpha.f = 0.5))
  
  lines(x = 1:length(seq(15, 32, by = 0.1)), 
    y = R0_curves_through_time_CrI[tt, 3, ], col = 'navy', lwd = 3)
  
  mtext(side = 3, text = paste0((when_to_assimilate[tt] - week_first_case), ' weeks assimilated'),
        cex = 0.75)
}
mtext(side = 2, outer = T, text = 'Dynamic R')
mtext(side = 1, outer = T, text = 'Temperature (c)', line = 0.5)
dev.off()


#=============================================================================#
# R0 by department by month
#=============================================================================#

dates_in_df = seq(as.Date("2015-08-09"), as.Date("2016-10-02"), by = 'week')[-61]

# do it for one time
particles_to_consider = array(NA, dim = c(6, np, ncol(particles_100_4)))
particles_to_consider[1,,] = as.matrix(particles_84_4)
particles_to_consider[2,,] = as.matrix(particles_100_4)
particles_to_consider[3,,] = as.matrix(particles_116_4)
particles_to_consider[4,,] = as.matrix(particles_132_4)
particles_to_consider[5,,] = as.matrix(particles_108_4)
particles_to_consider[6,,] = as.matrix(particles_120_4)

R0_by_month_by_dept = list()
length(R0_by_month_by_dept) = dim(particles_to_consider)[1]

for(pp_c in 1:dim(particles_to_consider)[1]){
  R0_by_month_by_dept[[pp_c]] = array(NA, dim = c(np, 12, length(dept_names)))
  
  particles = particles_to_consider[pp_c,,]
  colnames(particles) = colnames(particles_100_4)
  particles = as.data.frame(particles)
  
  for(mm in 1:12){
    for(dd in 1:length(dept_names)){
      
      wk_ind = which(month(dates_in_df) == mm)
      month_temp = mean(df_dept_temp_mosq$temp_mean[which(df_dept_temp_mosq$dept_name == dept_names[dd] & 
                                                  df_dept_temp_mosq$week %in% wk_ind)])
      month_mosq = mean(df_dept_temp_mosq$mosq_abund[which(df_dept_temp_mosq$dept_name == dept_names[dd] & 
                                                   df_dept_temp_mosq$week %in% wk_ind)])
      
      for(pp in 1:np){
        R0_by_month_by_dept[[pp_c]][pp, mm, dd] = particles$k_scale[pp] * 
          calc_R0t(temp_in = month_temp, 
                   mosq_abund_in = month_mosq,
                   k_in = particles$k[pp], sn1_in = particles$sn1[pp],
                   sn2_in = particles$sn2[pp], sn3_in = particles$sn3[pp])
      }
    }
  }
}



# STATIC R0 by department
particles_to_consider_static = array(NA, dim = c(6, np, ncol(particles_100_2)))
particles_to_consider_static[1,,] = as.matrix(particles_84_2)
particles_to_consider_static[2,,] = as.matrix(particles_100_2)
particles_to_consider_static[3,,] = as.matrix(particles_116_2)
particles_to_consider_static[4,,] = as.matrix(particles_132_2)
particles_to_consider_static[5,,] = as.matrix(particles_108_2)
particles_to_consider_static[6,,] = as.matrix(particles_120_2)

R0_static_by_dept = list()
length(R0_static_by_dept) = dim(particles_to_consider_static)[1]

for(pp_c in 1:dim(particles_to_consider_static)[1]){
  R0_static_by_dept[[pp_c]] = array(NA, dim = c(np, length(dept_names)))
  
  particles = particles_to_consider_static[pp_c,,]
  colnames(particles) = colnames(particles_100_2)
  particles = as.data.frame(particles)
  
  for(pp in 1:np){
    R0_static_by_dept[[pp_c]][pp, ] = particles$k[pp] * df_stat$R0[1:32]
  }
}

save(R0_static_by_dept, R0_by_month_by_dept, R0_curves_through_time, R0_curves_through_time_CrI,
     file = '../output/main_text_figures/R_output_for_RvsRt_figure.RData')



depts_to_plot = c(18, 7, 11)

cols_to_plot_static = c('slategray1', 'slategray3', 'slategray4', 'navy')
cols_to_plot_dynamic = c('indianred1', 'indianred3', 'indianred4', 'firebrick4')

maxs_by_dept = c(7, 30, 15)
times_assimilated_to_plot = c('Initial forecast', '16 weeks assimilated',
                              '32 weeks assimilated', '48 weeks assimilated')

pdf('../output/ensemble/monthly_R_through_time_select_depts.pdf',
    height = 7, width = 6)
par(mar = c(2,2,2,1), oma = c(3, 3, 3, 0))
for(dd in depts_to_plot){
  layout(matrix(1:4, nrow = 4))
  
  if(dd != 22){
    
    
    for(pp_c in 1:dim(particles_to_consider)[1]){
      plot(-100, -100, xlim = c(1, 24), ylim = c(0, maxs_by_dept[which(depts_to_plot == dd)]),
           xlab = '', ylab = '', xaxt = 'n', las = 1)
      axis(side = 1, at = seq(2, 24, length.out = 12) - 0.5, labels = month.abb)
      abline(h = 1, lwd = 2, col = rgb(0,0,0,0.5))
      mtext(side = 3, times_assimilated_to_plot[pp_c])
      mtext(side = 2, text = 'R', line = 2.5)
      
      if(pp_c == 1){
        mtext(side = 3, outer = T, line = 1, text = dept_names[dd])
      }
      
      l = 1
      for(mm in 1:12){
        boxplot(R0_by_month_by_dept[[pp_c]][,mm,dd], add = T, at = l + 0.25, 
                col = adjustcolor(cols_to_plot_dynamic[1], alpha.f = 1)
                ,
                outcol = adjustcolor(cols_to_plot_dynamic[2], alpha.f = 0.1)
                , outpch = 20
                , border = adjustcolor(cols_to_plot_dynamic[3], alpha.f = 0.9)
                , medcol = adjustcolor(cols_to_plot_dynamic[4], alpha.f = 1)
                , pars = list(xaxt = 'n', yaxt = 'n', axes = F))
        l = l + 1
        
        boxplot(R0_static_by_dept[[pp_c]][,dd], add = T, at = l - 0.25, 
                col = adjustcolor(cols_to_plot_static[1], alpha.f = 1)
                ,
                outcol = adjustcolor(cols_to_plot_static[2], alpha.f = 0.1)
                , outpch = 20
                , border = adjustcolor(cols_to_plot_static[3], alpha.f = 0.9)
                , medcol = adjustcolor(cols_to_plot_static[4], alpha.f = 1)
                , pars = list(xaxt = 'n', yaxt = 'n', axes = F))
        l = l + 1
      }
    }
    mtext(side = 1, outer = T, text = 'Month', line = 1)
  }
}
dev.off()

