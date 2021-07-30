#=============================================================================#
# Author: Rachel Oidtman

# Supplementary Figure
# parameter correlations and parameter values violin plots through time
# normalized posterior - prior figures
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
runs = c(1, 9, 13)

# load in temp and mosqsuito occurrence probabiity data
temp_mean = read.csv('../data/dept_tmean_weekly_weighted.csv', stringsAsFactors=F)
aedes_op = read.csv('../data/aegypti_pop_monthly_col_dept.csv', stringsAsFactors= F)


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

#=============================================================================#
# necessary functions
#=============================================================================#

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
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
# load in particles
#=============================================================================#


for(rr in 1:length(runs)){
  
  if(runs[rr] %in% c(1:4, 9:16)){
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
np = dim(particles_100_1)[1]


#=============================================================================#
# R MULTIPLIER, REPORTING RATE, OVERDISPERSION
# particles through time for spatial models
# static r: k, rho, dispersion, tm_ini, loc_ini
# dynamic r: c, k, rho, dispersion, tm_ini, loc_ini, sn1, sn2, sn3
#=============================================================================#

library(vioplot)
library(viridis)

pdf(file = '../output/main_text_figures/supplement/params_through_time_R-1-spatial.pdf',
    height = 5, width = 6)

times_to_plot = seq(1, length(when_to_assimilate), by = 3)
model_names = c('R-1-CDRs','R-1-gravity','R-1-radiation')
cols = viridis(5)

layout(matrix(1:12, nrow = 3, ncol = 4, byrow = T))
par(mar = c(2,2,1,1), oma = c(2, 2, 1,1))

for(mm in runs){
  
  if(mm %in% c(1,9,13)){
    particles = array(NA, dim = c(length(when_to_assimilate), np, 5))
    
    for(tt in 1:length(when_to_assimilate)){
      eval(parse(text = paste0(
        'particles_tmp = particles_', when_to_assimilate[tt], '_', mm
      )))
      
      particles[tt,,] = as.matrix(particles_tmp)
    } 
  }
  
  
  # k
  plot(-100, -100, xlim = c(1-0.25, length(times_to_plot)+0.25), ylim = c(0, 15), las = 1, xaxt = 'n',
       xlab = '', ylab = '')
  abline(h = c(2.5, 5, 7.5, 10, 12.5, 15), col = rgb(0,0,0,0.5), lwd = 0.75)
  abline(h = c(2.5-1.25, 2.5 + 1.25, 5 + 1.25, 7.5 + 1.25,
               10 + 1.25, 12.5 + 1.25), col = rgb(0,0,0,0.3), lwd = 0.5)
  axis(side = 1, at = 1:length(times_to_plot), 
       labels = when_to_assimilate[times_to_plot] - week_first_case)
  
  if(mm == 1){
    mtext(side = 3, text = 'R multiplier', cex = 0.8)
  }
  for(tt in times_to_plot){
    vioplot(particles[tt,,1], at = which(tt == times_to_plot), add = T, 
            col = cols[which(tt == times_to_plot)])
  }
  
  
  # rho
  plot(-100, -100, xlim = c(1-0.25, length(times_to_plot)+0.25), ylim = c(0, 1), las = 1, xaxt = 'n',
       xlab = '', ylab = '')
  abline(h = c(0.2, 0.4, 0.6, 0.8, 1), col = rgb(0,0,0,0.5), lwd = 0.75)
  abline(h = c(0.1, 0.3, 0.5, 0.7, 0.9), col = rgb(0,0,0,0.3), lwd = 0.5)
  axis(side = 1, at = 1:length(times_to_plot), 
       labels = when_to_assimilate[times_to_plot] - week_first_case)
  if(mm == 1){
    mtext(side = 3, text = 'Reporting rate', cex = 0.8)
  }
  for(tt in times_to_plot){
    vioplot(particles[tt,,2], at = which(tt == times_to_plot), add = T, 
            col = cols[which(tt == times_to_plot)])
  }
  
  # dispersion
  plot(-100, -100, xlim = c(1-0.25, length(times_to_plot)+0.25), ylim = c(0, 5), las = 1, xaxt = 'n',
       xlab = '', ylab = '')
  abline(h = c(1:5), col = rgb(0,0,0,0.5), lwd = 0.75)
  abline(h = c(0.5, 1.5, 2.5, 3.5, 4.5), col = rgb(0,0,0,0.3), lwd = 0.5)
  axis(side = 1, at = 1:length(times_to_plot), 
       labels = when_to_assimilate[times_to_plot] - week_first_case)
  if(mm == 1){
    mtext(side = 3, text = 'Overdispersion', cex = 0.8)
  }
  for(tt in times_to_plot){
    vioplot(particles[tt,,3], at = which(tt == times_to_plot), add = T, 
            col = cols[which(tt == times_to_plot)])
  }
  
  # tm_ini
  # tm_ini
  plot(-100, -100, xlim = c(1-0.25, length(times_to_plot)+0.25),
       ylim = c(min(particles[,,4]),
                max(particles[,,4])),
       las = 1, xaxt = 'n', yaxt = 'n', ylab = '', xlab = '')
  abline(h = seq(50, 85, by = 5), col = rgb(0,0,0,0.5), lwd = 0.75)
  abline(h = seq(50, 85, by = 5) + 2.5, col = rgb(0,0,0,0.3), lwd = 0.5)
  # mtext(side = 3, text = 'Week of initial importations', cex = 0.8)
  axis(side = 2, at = seq(50, 85, by = 5), labels = seq(50, 85, by = 5) - week_first_case - 1, las = 1)
  axis(side = 1, at = 1:length(times_to_plot), 
       labels = when_to_assimilate[times_to_plot] - week_first_case)
  if(mm == 1){
    mtext(side = 3, text = 'Wk. of 1st infections', cex = 0.8)
  }
  for(tt in times_to_plot){
    vioplot(particles[tt,,4], at = which(tt == times_to_plot), add = T, 
            col = cols[which(tt == times_to_plot)])
  }
  mtext(side = 4, text = model_names[which(mm == runs)], line = 0.5)
}

mtext(side = 1, text = 'Weeks since 1st reported case in Colombia',outer = T, line = 0.25)
mtext(side = 2, text = 'Parameter value',outer = T, line = 0)
dev.off()



#=============================================================================#
# R MULTIPLIER, REPORTING RATE, OVERDISPERSION - logged
# particles through time for spatial models
# static r: k, rho, dispersion, tm_ini, loc_ini
#=============================================================================#

library(vioplot)
library(viridis)

pdf(file = '../output/main_text_figures/supplement/params_through_time_R-1-spatial_logged.pdf',
    height = 5, width = 6)

times_to_plot = seq(1, length(when_to_assimilate), by = 3)
model_names = c('R-1-CDRs','R-1-gravity','R-1-radiation')
cols = viridis(5)

layout(matrix(1:12, nrow = 3, ncol = 4, byrow = T))
par(mar = c(2,2,1,1), oma = c(2, 3, 1,1))

for(mm in runs){
  
  if(mm %in% c(1,9,13)){
    particles = array(NA, dim = c(length(when_to_assimilate), np, 5))
    
    for(tt in 1:length(when_to_assimilate)){
      eval(parse(text = paste0(
        'particles_tmp = particles_', when_to_assimilate[tt], '_', mm
      )))
      
      particles[tt,,] = as.matrix(particles_tmp)
    } 
  }
  
  
  # k
  plot(-100, -100, xlim = c(1-0.25, length(times_to_plot)+0.25), ylim = c(log10(0.01), log10(30)), las = 1, xaxt = 'n',
       xlab = '', ylab = '')
  # abline(h = log10(c(2.5, 5, 7.5, 10, 12.5, 15)), col = rgb(0,0,0,0.5), lwd = 0.75)
  # abline(h = log10(c(2.5-1.25, 2.5 + 1.25, 5 + 1.25, 7.5 + 1.25,
  #              10 + 1.25, 12.5 + 1.25)), col = rgb(0,0,0,0.3), lwd = 0.5)
  axis(side = 1, at = 1:length(times_to_plot), 
       labels = when_to_assimilate[times_to_plot] - week_first_case)
  
  if(mm == 1){
    mtext(side = 3, text = 'R multiplier', cex = 0.8)
  }
  for(tt in times_to_plot){
    vioplot(log10(particles[tt,,1]), at = which(tt == times_to_plot), add = T, 
            col = cols[which(tt == times_to_plot)])
  }
  
  
  # rho
  plot(-100, -100, xlim = c(1-0.25, length(times_to_plot)+0.25), 
       ylim = c(log10(0.0000001), 0), las = 1, xaxt = 'n',
       xlab = '', ylab = '')
  # abline(h = c(0.2, 0.4, 0.6, 0.8, 1), col = rgb(0,0,0,0.5), lwd = 0.75)
  # abline(h = c(0.1, 0.3, 0.5, 0.7, 0.9), col = rgb(0,0,0,0.3), lwd = 0.5)
  axis(side = 1, at = 1:length(times_to_plot), 
       labels = when_to_assimilate[times_to_plot] - week_first_case)
  if(mm == 1){
    mtext(side = 3, text = 'Reporting rate', cex = 0.8)
  }
  for(tt in times_to_plot){
    vioplot(log10(particles[tt,,2]), at = which(tt == times_to_plot), add = T, 
            col = cols[which(tt == times_to_plot)])
  }
  
  
  # dispersion
  plot(-100, -100, xlim = c(1-0.25, length(times_to_plot)+0.25), 
       ylim = c(log10(0.001), log10(5)), las = 1, xaxt = 'n',
       xlab = '', ylab = '')
  # abline(h = c(1:5), col = rgb(0,0,0,0.5), lwd = 0.75)
  # abline(h = c(0.5, 1.5, 2.5, 3.5, 4.5), col = rgb(0,0,0,0.3), lwd = 0.5)
  axis(side = 1, at = 1:length(times_to_plot), 
       labels = when_to_assimilate[times_to_plot] - week_first_case)
  if(mm == 1){
    mtext(side = 3, text = 'Overdispersion', cex = 0.8)
  }
  for(tt in times_to_plot){
    vioplot(log10(particles[tt,,3]), at = which(tt == times_to_plot), add = T, 
            col = cols[which(tt == times_to_plot)])
  }
  
  # tm_ini
  # tm_ini
  plot(-100, -100, xlim = c(1-0.25, length(times_to_plot)+0.25),
       ylim = c(min(particles[,,4]),
                max(particles[,,4])),
       las = 1, xaxt = 'n', yaxt = 'n', ylab = '', xlab = '')
  # abline(h = seq(50, 85, by = 5), col = rgb(0,0,0,0.5), lwd = 0.75)
  # abline(h = seq(50, 85, by = 5) + 2.5, col = rgb(0,0,0,0.3), lwd = 0.5)
  # mtext(side = 3, text = 'Week of initial importations', cex = 0.8)
  axis(side = 2, at = seq(50, 85, by = 5), labels = seq(50, 85, by = 5) - week_first_case - 1, las = 1)
  axis(side = 1, at = 1:length(times_to_plot), 
       labels = when_to_assimilate[times_to_plot] - week_first_case)
  if(mm == 1){
    mtext(side = 3, text = 'Wk. of 1st infections', cex = 0.8)
  }
  for(tt in times_to_plot){
    vioplot(particles[tt,,4], at = which(tt == times_to_plot), add = T, 
            col = cols[which(tt == times_to_plot)])
  }
  mtext(side = 4, text = model_names[which(mm == runs)], line = 0.5)
}

mtext(side = 1, text = 'Weeks since 1st reported case in Colombia',outer = T, line = 0.25)
mtext(side = 2, text = 'Parameter value (log10)',outer = T, line = 1)
dev.off()




#=============================================================================#
# correlation between paramters
#=============================================================================#
library(corrplot)
library(RColorBrewer)

pdf('../output/main_text_figures/supplement/params_correlation_though_time_R-1-spatial.pdf',
    height = 7, width = 9)


times_to_plot = seq(1, length(when_to_assimilate), by = 3)
model_names = c('R-1-CDRs','R-1-gravity','R-1-radiation')
par(mar = c(2,4,4,1))
layout(matrix(1:15, nrow = 3, ncol = 5, byrow = T))
for(mm in runs){
  
  if(mm %in% c(1,9,13)){
    particles = array(NA, dim = c(length(when_to_assimilate), np, 5))
    
    for(tt in 1:length(when_to_assimilate)){
      eval(parse(text = paste0(
        'particles_tmp = particles_', when_to_assimilate[tt], '_', mm
      )))
      
      particles[tt,,] = as.matrix(particles_tmp)
    } 
  }
  
  for(tt in times_to_plot){
    C = cor(particles[tt,,1:4])
    row.names(C) = c('R mulitplier', 'reporting rate', 'overdispersion', 'wk. of 1st inf.')
    colnames(C) = c('R mulitplier', 'reporting rate', 'overdispersion', 'wk. of 1st inf.')
    corrplot(C, type = 'upper', method = 'color', tl.col = 'black', cl.pos = 'n'
             # ,
             # col = brewer.pal(n = 11, name = 'RdYlBu')
             # col = cols
             )
  }
}

dev.off()


#=============================================================================#
# Image plot with posterior versus prior
# prior - posterior
# for nonspatial models, look at median value for parameters
#=============================================================================#

prior_df = data.frame(k = rep(NA, length(runs)),
                      rho = rep(NA, length(runs)),
                      dispersion = rep(NA, length(runs)),
                      tm_ini_1 = rep(NA, length(runs)),
                      loc_ini_1 = rep(NA, length(runs)),
                      tm_ini_2 = rep(NA, length(runs)),
                      loc_ini_2 = rep(NA, length(runs)),
                      c = rep(NA, length(runs)),
                      sn1 = rep(NA, length(runs)),
                      sn2 = rep(NA, length(runs)),
                      sn3 = rep(NA, length(runs)))

posterior_df = data.frame(k = rep(NA, length(runs)),
                          rho = rep(NA, length(runs)),
                          dispersion = rep(NA, length(runs)),
                          tm_ini_1 = rep(NA, length(runs)),
                          loc_ini_1 = rep(NA, length(runs)),
                          tm_ini_2 = rep(NA, length(runs)),
                          loc_ini_2 = rep(NA, length(runs)),
                          c = rep(NA, length(runs)),
                          sn1 = rep(NA, length(runs)),
                          sn2 = rep(NA, length(runs)),
                          sn3 = rep(NA, length(runs)))

for(mm in c(1:16)){
  
  if(mm %in% c(1,9,13)){
    particles = array(NA, dim = c(length(when_to_assimilate), np, 5))
    
    for(tt in 1:length(when_to_assimilate)){
      eval(parse(text = paste0(
        'particles_tmp = particles_', when_to_assimilate[tt], '_', mm
      )))
      
      particles[tt,,] = as.matrix(particles_tmp)
    } 
    
    # data frame wiht priors and posteriors
    prior_df$k[mm] = median(particles[1,,1])
    posterior_df$k[mm] = median(particles[length(when_to_assimilate),,1])
    
    prior_df$rho[mm] = median(particles[1,,2])
    posterior_df$rho[mm] = median(particles[length(when_to_assimilate),,2])
    
    prior_df$dispersion[mm] = median(particles[1,,3])
    posterior_df$dispersion[mm] = median(particles[length(when_to_assimilate),,3])
    
    prior_df$tm_ini_1[mm] = getmode(particles[1,,4])
    posterior_df$tm_ini_1[mm] = getmode(particles[length(when_to_assimilate),,4])
    
    prior_df$loc_ini_2[mm] = getmode(particles[1,,5])
    posterior_df$loc_ini_2[mm] = getmode(particles[length(when_to_assimilate),,5])
    
  }else if(mm %in% c(2, 10, 14)){
    particles = array(NA, dim = c(length(when_to_assimilate), np, 7))
    
    for(tt in 1:length(when_to_assimilate)){
      eval(parse(text = paste0(
        'particles_tmp = particles_', when_to_assimilate[tt], '_', mm
      )))
      
      particles[tt,,] = as.matrix(particles_tmp)
    }
    
    # data frame wiht priors and posteriors
    prior_df$k[mm] = median(particles[1,,1])
    posterior_df$k[mm] = median(particles[length(when_to_assimilate),,1])
    
    prior_df$rho[mm] = median(particles[1,,2])
    posterior_df$rho[mm] = median(particles[length(when_to_assimilate),,2])
    
    prior_df$dispersion[mm] = median(particles[1,,3])
    posterior_df$dispersion[mm] = median(particles[length(when_to_assimilate),,3])
    
    prior_df$tm_ini_1[mm] = getmode(particles[1,,4])
    posterior_df$tm_ini_1[mm] = getmode(particles[length(when_to_assimilate),,4])
    
    prior_df$loc_ini_1[mm] = getmode(particles[1,,5])
    posterior_df$loc_ini_1[mm] = getmode(particles[length(when_to_assimilate),,5])
    
    prior_df$tm_ini_2[mm] = getmode(particles[1,,6])
    posterior_df$tm_ini_2[mm] = getmode(particles[length(when_to_assimilate),,6])
    
    prior_df$loc_ini_2[mm] = getmode(particles[1,,7])
    posterior_df$loc_ini_2[mm] = getmode(particles[length(when_to_assimilate),,7])
    
  }else if(mm %in% c(3, 11, 15)){
    particles = array(NA, dim = c(length(when_to_assimilate), np, 9))
    
    for(tt in 1:length(when_to_assimilate)){
      eval(parse(text = paste0(
        'particles_tmp = particles_', when_to_assimilate[tt], '_', mm
      )))
      
      particles[tt,,] = as.matrix(particles_tmp)
    } 
    
    # data frame wiht priors and posteriors
    prior_df$k[mm] = median(particles[1,,2])
    posterior_df$k[mm] = median(particles[length(when_to_assimilate),,2])
    
    prior_df$rho[mm] = median(particles[1,,3])
    posterior_df$rho[mm] = median(particles[length(when_to_assimilate),,3])
    
    prior_df$dispersion[mm] = median(particles[1,,4])
    posterior_df$dispersion[mm] = median(particles[length(when_to_assimilate),,4])
    
    prior_df$tm_ini_1[mm] = getmode(particles[1,,5])
    posterior_df$tm_ini_1[mm] = getmode(particles[length(when_to_assimilate),,5])
    
    prior_df$loc_ini_1[mm] = getmode(particles[1,,6])
    posterior_df$loc_ini_1[mm] = getmode(particles[length(when_to_assimilate),,6])
    
    prior_df$c[mm] = median(particles[1,, 1])
    posterior_df$c[mm] = median(particles[length(when_to_assimilate),, 1])
    
    prior_df$sn1[mm] = median(particles[1,,7])
    posterior_df$sn1[mm] = median(particles[length(when_to_assimilate),,7])
    
    prior_df$sn2[mm] = median(particles[1,,8])
    posterior_df$sn2[mm] = median(particles[length(when_to_assimilate),,8])
    
    prior_df$sn3[mm] = median(particles[1,,9])
    posterior_df$sn3[mm] = median(particles[length(when_to_assimilate),,9])
    
  }else if(mm %in% c(4, 12, 16)){
    particles = array(NA, dim = c(length(when_to_assimilate), np, 11))
    
    for(tt in 1:length(when_to_assimilate)){
      eval(parse(text = paste0(
        'particles_tmp = particles_', when_to_assimilate[tt], '_', mm
      )))
      
      particles[tt,,] = as.matrix(particles_tmp)
    }
    
    # data frame wiht priors and posteriors
    prior_df$k[mm] = median(particles[1,,2])
    posterior_df$k[mm] = median(particles[length(when_to_assimilate),,2])
    
    prior_df$rho[mm] = median(particles[1,,3])
    posterior_df$rho[mm] = median(particles[length(when_to_assimilate),,3])
    
    prior_df$dispersion[mm] = median(particles[1,,4])
    posterior_df$dispersion[mm] = median(particles[length(when_to_assimilate),,4])
    
    prior_df$tm_ini_1[mm] = getmode(particles[1,,5])
    posterior_df$tm_ini_1[mm] = getmode(particles[length(when_to_assimilate),,5])
    
    prior_df$loc_ini_1[mm] = getmode(particles[1,,6])
    posterior_df$loc_ini_1[mm] = getmode(particles[length(when_to_assimilate),,6])
    
    prior_df$tm_ini_2[mm] = getmode(particles[1,,7])
    posterior_df$tm_ini_2[mm] = getmode(particles[length(when_to_assimilate),,7])
    
    prior_df$loc_ini_2[mm] = getmode(particles[1,,8])
    posterior_df$loc_ini_2[mm] = getmode(particles[length(when_to_assimilate),,8])
    
    prior_df$c[mm] = median(particles[1,, 1])
    posterior_df$c[mm] = median(particles[length(when_to_assimilate),, 1])
    
    prior_df$sn1[mm] = median(particles[1,,9])
    posterior_df$sn1[mm] = median(particles[length(when_to_assimilate),,9])
    
    prior_df$sn2[mm] = median(particles[1,,10])
    posterior_df$sn2[mm] = median(particles[length(when_to_assimilate),,10])
    
    prior_df$sn3[mm] = median(particles[1,,11])
    posterior_df$sn3[mm] = median(particles[length(when_to_assimilate),,11])
    
  }else if(mm == 5){
    particles = list()
    length(particles) = length(dept_names)
    
    for(dd in 1:length(dept_names)){
      if(dd != 22){
        particles[[dd]] = array(NA, dim = c(length(when_to_assimilate_NS[[dd]]),
                                            np, 4))
        for(tt in 1:length(when_to_assimilate_NS[[dd]])){
          
          eval(parse(text = paste0(
            'particles_tmp = particles_', dd, '_', when_to_assimilate_NS[[dd]][tt], '_', runs[mm]
          )))
          
          particles[[dd]][tt,,] = as.matrix(particles_tmp)
        }
      }
    }
    
    # median particles
    depts_to_consider = c(1:21, 23:32)
    
    prior_df$k[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][1,,1]))
    posterior_df$k[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,1]))
    
    prior_df$rho[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][1,,2]))
    posterior_df$rho[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,2]))
    
    prior_df$dispersion[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][1,,3]))
    posterior_df$dispersion[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,3]))
    
    prior_df$tm_ini_1[mm] = getmode(sapply(depts_to_consider, function(ff) particles[[ff]][1,,4]))
    posterior_df$tm_ini_1[mm] = getmode(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,4]))
    
  }else if(mm == 6){
    particles = list()
    length(particles) = length(dept_names)
    
    for(dd in 1:length(dept_names)){
      if(dd != 22){
        particles[[dd]] = array(NA, dim = c(length(when_to_assimilate_NS[[dd]]),
                                            np, 5))
        for(tt in 1:length(when_to_assimilate_NS[[dd]])){
          
          eval(parse(text = paste0(
            'particles_tmp = particles_', dd, '_', when_to_assimilate_NS[[dd]][tt], '_', runs[mm]
          )))
          
          particles[[dd]][tt,,] = as.matrix(particles_tmp)
        }
      }
    }
    
    # median particles
    depts_to_consider = c(1:21, 23:32)
    
    prior_df$k[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][1,,1]))
    posterior_df$k[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,1]))
    
    prior_df$rho[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][1,,2]))
    posterior_df$rho[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,2]))
    
    prior_df$dispersion[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][1,,3]))
    posterior_df$dispersion[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,3]))
    
    prior_df$tm_ini_1[mm] = getmode(sapply(depts_to_consider, function(ff) particles[[ff]][1,,4]))
    posterior_df$tm_ini_1[mm] = getmode(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,4]))
    
    prior_df$tm_ini_2[mm] = getmode(sapply(depts_to_consider, function(ff) particles[[ff]][1,,5]))
    posterior_df$tm_ini_2[mm] = getmode(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,5]))
    
  }else if(mm == 7){
    particles = list()
    length(particles) = length(dept_names)
    
    for(dd in 1:length(dept_names)){
      if(dd != 22){
        particles[[dd]] = array(NA, dim = c(length(when_to_assimilate_NS[[dd]]),
                                            np, 8))
        for(tt in 1:length(when_to_assimilate_NS[[dd]])){
          
          eval(parse(text = paste0(
            'particles_tmp = particles_', dd, '_', when_to_assimilate_NS[[dd]][tt], '_', runs[mm]
          )))
          
          particles[[dd]][tt,,] = as.matrix(particles_tmp)
        }
      }
    }
    
    # median particles
    depts_to_consider = c(1:21, 23:32)
    
    prior_df$k[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][1,,2]))
    posterior_df$k[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,2]))
    
    prior_df$rho[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][1,,3]))
    posterior_df$rho[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,3]))
    
    prior_df$dispersion[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][1,,4]))
    posterior_df$dispersion[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,4]))
    
    prior_df$tm_ini_1[mm] = getmode(sapply(depts_to_consider, function(ff) particles[[ff]][1,,5]))
    posterior_df$tm_ini_1[mm] = getmode(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,5]))
    
    prior_df$c[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][1,,1]))
    posterior_df$c[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,1]))
    
    prior_df$sn1[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][1,,6]))
    posterior_df$sn1[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,6]))
    
    prior_df$sn2[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][1,,7]))
    posterior_df$sn2[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,7]))
    
    prior_df$sn3[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][1,,8]))
    posterior_df$sn3[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,8]))
    
  }else if(mm == 8){
    particles = list()
    length(particles) = length(dept_names)
    
    for(dd in 1:length(dept_names)){
      if(dd != 22){
        particles[[dd]] = array(NA, dim = c(length(when_to_assimilate_NS[[dd]]),
                                            np, 9))
        for(tt in 1:length(when_to_assimilate_NS[[dd]])){
          
          eval(parse(text = paste0(
            'particles_tmp = particles_', dd, '_', when_to_assimilate_NS[[dd]][tt], '_', runs[mm]
          )))
          
          particles[[dd]][tt,,] = as.matrix(particles_tmp)
        }
      }
    }
    
    # median particles
    depts_to_consider = c(1:21, 23:32)
    
    prior_df$k[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][1,,2]))
    posterior_df$k[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,2]))
    
    prior_df$rho[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][1,,3]))
    posterior_df$rho[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,3]))
    
    prior_df$dispersion[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][1,,4]))
    posterior_df$dispersion[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,4]))
    
    prior_df$tm_ini_1[mm] = getmode(sapply(depts_to_consider, function(ff) particles[[ff]][1,,5]))
    posterior_df$tm_ini_1[mm] = getmode(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,5]))
    
    prior_df$tm_ini_2[mm] = getmode(sapply(depts_to_consider, function(ff) particles[[ff]][1,,6]))
    posterior_df$tm_ini_2[mm] = getmode(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,6]))
    
    prior_df$c[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][1,,1]))
    posterior_df$c[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,1]))
    
    prior_df$sn1[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][1,,7]))
    posterior_df$sn1[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,7]))
    
    prior_df$sn2[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][1,,8]))
    posterior_df$sn2[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,8]))
    
    prior_df$sn3[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][1,,9]))
    posterior_df$sn3[mm] = median(sapply(depts_to_consider, function(ff) particles[[ff]][length(when_to_assimilate_NS[[ff]]),,9]))
  }
}


prior_post_diff = matrix(NA, nrow = length(runs), ncol = 11)
for(ii in 1:nrow(prior_post_diff)){
  prior_post_diff[ii, ] = as.numeric(posterior_df[ii, ]) - as.numeric(prior_df[ii, ])
}

prior_post_diff_to_plot = prior_post_diff[,-c(4:7)]

# DIVIDE BY THE PRIOR
prior_post_diff_to_plot = prior_post_diff[,-c(4:7)]
prior_post_diff_to_plot[,1] = prior_post_diff_to_plot[,1] / prior_df$k
prior_post_diff_to_plot[,2] = prior_post_diff_to_plot[,2] / prior_df$rho
prior_post_diff_to_plot[,3] = prior_post_diff_to_plot[,3] / prior_df$dispersion
prior_post_diff_to_plot[,4] = prior_post_diff_to_plot[,4] / prior_df$c
prior_post_diff_to_plot[,5] = prior_post_diff_to_plot[,5] / prior_df$sn1
prior_post_diff_to_plot[,6] = prior_post_diff_to_plot[,6] / prior_df$sn2
prior_post_diff_to_plot[,7] = prior_post_diff_to_plot[,7] / prior_df$sn3


model_names = c('R-1-CDRs', 'R-2-CDRs', 'Rt-1-CDRs', 'Rt-2-CDRs',
                'R-1-nonspatial', 'R-2-nonspatial', 'Rt-1-nonspatial', 'Rt-2-nonspatial',
                'R-1-gravity', 'R-2-gravity', 'Rt-1-gravity', 'Rt-2-gravity', 
                'R-1-radiation', 'R-2-radiation', 'Rt-1-radiation', 'Rt-2-radiation')
param_names = c('R multiplier', 'Reporting rate', 'Overdispersion', 'Rt scalar', 
                'location (skew normal)', 'scale (skew normal)', 'shape (skew normal)')


pdf('../output/main_text_figures/supplement/posterior-prior_parameters_normalized.pdf', height = 5, width = 5.5)
layout(matrix(c(1,1,1,1,1,1,2), nrow = 1))
par(mar = c(1,1,1,0), oma = c(10,8.5,1,3))
cols = center_data_fxn(vec_in = as.vector(prior_post_diff_to_plot), browse = F,
                       col_less_in = 'brown1', col_greater_in = 'cornflowerblue', by_in = 0.01, add_in = 0.01, col_middle_in = 'gray95')
image(t(prior_post_diff_to_plot), col = cols, xaxt = 'n', yaxt = 'n')
box()
axis(side = 1, at = seq(0, 1, length.out = ncol(prior_post_diff_to_plot)), labels = param_names, las = 2)
axis(side = 2, at = seq(0, 1, length.out = length(runs)), labels = model_names, las = 1)
mtext(side = 2, line = 8, text = 'Models', cex = 0.75)
mtext(side = 1, line = 9.5, text = 'Parameters', cex = 0.75)
par(mar = c(0,1,0,0))
color.bar(cols, min = -1)
axis(side = 4, at = c(-1, 1), labels = signif(c(min(prior_post_diff_to_plot, na.rm = T),
                                                max(prior_post_diff_to_plot, na.rm = T)), 2), las = 1)
mtext(side = 3, 'Posterior - Prior', line = -0.5, cex = 0.75)
dev.off()
