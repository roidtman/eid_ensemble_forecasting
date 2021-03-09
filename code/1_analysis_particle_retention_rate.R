#=============================================================================#
# Author: Rachel Oidtman
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
# you can update this to change the run number to consider
runs = 1

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
      # particles = particles[sample(1:nrow(particles), 2000), ]
      
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
        # particles = particles[sample(1:nrow(particles), 2000), ]
        
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
# number unique particles through time
#=============================================================================#

unique_particles_through_time = rep(NA, length(when_to_assimilate))
mm = 1
for(tt in 1:length(when_to_assimilate)){
  eval(parse(text = paste0(
    'particles_tmp = particles_', when_to_assimilate[tt], '_', mm
  )))
  print(when_to_assimilate[tt])
  unique_particles_through_time[tt] = nrow(unique(particles_tmp))
}

plot(unique_particles_through_time[2:length(when_to_assimilate)])



unique_retained_particles_through_time = rep(NA, length(when_to_assimilate))
mm = 1
for(tt in 1:length(when_to_assimilate)){
  eval(parse(text = paste0(
    'particles_tmp = particles_', when_to_assimilate[tt], '_', mm
  )))
  print(when_to_assimilate[tt])
  if(tt > 2){
    unique_retained_particles_through_time[tt] = nrow(unique(particles_tmp[1:18000,]))
  }else{
    unique_retained_particles_through_time[tt] = nrow(unique(particles_tmp))
  }
}

plot(unique_retained_particles_through_time[2:length(when_to_assimilate)])

#=============================================================================#
# particle retention
#=============================================================================#

# 2000 particles are generated 

unique_retention_rate = rep(NA, length(when_to_assimilate))
for(tt in 1:length(when_to_assimilate)){
  unique_retention_rate[tt] = unique_retained_particles_through_time[tt + 1] / unique_retained_particles_through_time[tt]
}

# layout(matrix(1:2, nrow = 2))
pdf('../output/main_text_figures/particle_retention/prop_original_retained_m1.pdf',
    height = 3.5, width = 5)
par(mar = c(5,5,1,1))
layout(1)
plot(unique_retained_particles_through_time[2:length(when_to_assimilate)], las = 1, ylim = c(0, 20000), 
     ylab = '', xlab = '', type = 'b', col = 'steelblue4', lwd = 3, yaxt = 'n')
# axis(side = 2, col = 'yellow3', las = 1, lwd = 3)
axis(side = 2, las = 1)
mtext(side = 2, line = 4, text = '# of unique, retained particles')
mtext(side = 1, line = 2.5, text = 'Assimilation period')

# par(new = T)
# plot(y = unique_retention_rate[2:(length(when_to_assimilate)-1)], type = 'b',
#      x = 2:(length(when_to_assimilate)-1), xlim = c(1, (length(when_to_assimilate)-1)), xlab = '', xaxt = 'n',
#      yaxt = 'n', bty = 'n', ylim = c(0, 1), ylab = '', col = 'steelblue4', lwd = 3)
# axis(side = 4, las = 1, col = 'steelblue4', lwd = 3)
# mtext(side = 4, line = 3, text = '# unique particles in t+1 / # unique particles in t')
dev.off()

#=============================================================================#
# how many original particles are left through time
#=============================================================================#

num_original_retained = rep(NA, length(3:length(when_to_assimilate)))

particles_original = particles_88_1
original_retained = list()
length(original_retained) = length(3:length(when_to_assimilate))

f = '../output/main_text_figures/particle_retention/original_retained_particles_m1.RData'
if(file.exists(f)){
  load(f)
}else{
  for(tt in 1:length(when_to_assimilate)){
    eval(parse(text = paste0(
      'particles_tmp = particles_', when_to_assimilate[tt], '_', mm
    )))
    
    original_retained[[which(tt == (3:length(when_to_assimilate)))]] = c()
    
    for(ii in 1:nrow(unique(particles_tmp[1:18000,]))){
      print(paste(tt, ii)) 
      if(sum(which(particles_original[,1] == unique(particles_tmp[1:18000,])[ii,1])) == 0){
        original_retained[[which(tt == 3:length(when_to_assimilate))]][ii] = NA
      }else{
        original_retained[[which(tt == 3:length(when_to_assimilate))]][ii] = 
          which(particles_original[,1] == unique(particles_tmp[1:18000,])[ii,1]) 
      }
    }
  }
  save(original_retained, file = f)
}


# 
pdf('../output/main_text_figures/particle_retention/prop_original_particles_m1.pdf',
    height = 3.5, width = 5)
layout(1)
par(mar = c(5,5,3,1))
plot(y = c(NA, sapply(1:length(original_retained), function(ff) length(which(!is.na(original_retained[[ff]])))) / 
       unique_retained_particles_through_time[3:length(when_to_assimilate)]),
     x = 1:14,
     type = 'b', lwd =3, col = 'steelblue4', las = 1, ylab = '', xlab = '',
     ylim = c(0,1))
mtext(side = 1, line = 2.5, text = 'Assimilation period')
mtext(side = 2, line = 2.75, text = 'Prop. of particles from original ensemble')
dev.off()


#=============================================================================#
# how many original particles are left through time - hexbin plot
#=============================================================================#
library(ggplot2)

tibble_full = tibble(`R multiplier` = numeric(),
                     `Reporting probability` = numeric(), 
                    `Overdispersion` = numeric(),
                     time = numeric())

dat_text = tibble(cor_k_rho = rep(NA, length(when_to_assimilate) - 1),
                  cor_k_disp = rep(NA, length(when_to_assimilate) - 1),
                  cor_disp_rho = rep(NA, length(when_to_assimilate) - 1),
                  time = rep(NA, length(when_to_assimilate) - 1))
l = 1
for(tt in 2:length(when_to_assimilate)){
  
  eval(parse(text = paste0(
    'particles_tmp = particles_', when_to_assimilate[tt], '_', mm
  )))
  
  df_tmp = tibble(`R multiplier` = particles_tmp$k,
                  `Reporting probability` = particles_tmp$rho,
                  `Overdispersion` = particles_tmp$dispersion,
                  time = when_to_assimilate[tt] - week_first_case)
  tibble_full = rbind(tibble_full, df_tmp)
  
  
  dat_text$cor_k_rho[l] = signif(cor(particles_tmp$k, particles_tmp$rho),3)
  dat_text$cor_k_disp[l] = signif(cor(particles_tmp$k, particles_tmp$dispersion),3)
  dat_text$cor_disp_rho[l] = signif(cor(particles_tmp$dispersion, particles_tmp$rho),3)
  dat_text$time[l] = when_to_assimilate[tt] - week_first_case
  l = l + 1
}



# k versus rho
d <- ggplot(tibble_full, aes(`R multiplier`, `Reporting probability`))
d + geom_hex(bins = 20) + 
  xlim(0, 30) + 
  ylim(0, 1) + theme_bw() + 
  facet_wrap(~time, nrow = 5) + 
  scale_fill_gradient(low="lightblue1",high="darkblue",trans="log10") + 
  geom_text(
    data = dat_text,
    mapping = aes(x = 27, y = 0.9, label = cor_k_rho),
    size = 2
  )
ggsave(filename="../output/main_text_figures/particle_retention/k_vs_rho_m1.pdf",
       width=6,height=6)



# k versus dispersion
d <- ggplot(tibble_full, aes(`R multiplier`, `Overdispersion`))
d + geom_hex(bins = 20) + 
  xlim(0, 30) + 
  ylim(0, 3) + theme_bw() + 
  facet_wrap(~time, nrow = 5) + 
  scale_fill_gradient(low="lightblue1",high="darkblue",trans="log10") + 
  geom_text(
    data = dat_text,
    mapping = aes(x = 27, y = 2.8, label = cor_k_disp),
    size =2
  )
ggsave(filename="../output/main_text_figures/particle_retention/k_vs_dispersion_m1.pdf",
       width=6,height=6)



# k versus rho
d <- ggplot(tibble_full, aes(`Overdispersion`, `Reporting probability`))
d + geom_hex(bins = 20) + 
  xlim(0, 3) + 
  ylim(0, 1) + theme_bw() + 
  facet_wrap(~time, nrow = 5) + 
  scale_fill_gradient(low="lightblue1",high="darkblue",trans="log10") + 
  geom_text(
    data = dat_text,
    mapping = aes(x = 2.7, y = 0.9, label = cor_disp_rho),
    size = 2
  )
ggsave(filename="../output/main_text_figures/particle_retention/dispersion_vs_rho_m1.pdf",
       width=6,height=6)



#=============================================================================#
# figure out retained states
#=============================================================================#


load('../output/no_tm_ini_1/I_F_processed.RData')
load('../output/no_tm_ini_1/I_A_processed.RData')

dim(I_F_over_time[[1]])

plot(I_F_over_time[[1]][1,7,])

total_infs_per_depts = matrix(NA, nrow = np, ncol = length(dept_names))
for(dd in 1:32){
  total_infs_per_depts[,dd] = sapply(1:np, function(ff) sum(I_F_over_time[[1]][dd,ff,]))
}



total_infs_per_run_per_time = matrix(NA, nrow = length(when_to_assimilate), ncol = np)
  
for(tt in when_to_assimilate){
  # if(tt > when_to_assimilate[1]){
  #   total_infs_per_run_per_time[which(when_to_assimilate == tt),] = sapply(1:np, function(ff) 
  #     sum(I_A_over_time[,ff,(when_to_assimilate[2] - week_first_case + 1) : (tt - week_first_case + time_btwn_assim)]))
  # }else{
    total_infs_per_run_per_time[which(when_to_assimilate == tt),] = sapply(1:np, function(ff) 
      sum(I_F_over_time[[which(when_to_assimilate == tt)]][,ff,]))
  # }
}

prop_forecasts_zero = sapply(1:length(when_to_assimilate), function(ff) 
  length(which(total_infs_per_run_per_time[ff,] == 0))) / np


pdf('../output/main_text_figures/particle_retention/prop_forecasts_died_out.pdf',
    height = 3.5, width = 5)
layout(1)
par(mar = c(5,5,3,1))
plot(prop_forecasts_zero, ylim = c(0, 1), las = 1, 
     xlab = '', ylab = '', type = 'b', col = 'steelblue4', lwd = 3)
mtext(side = 2, text = 'Prop. forecasts predicting zero infections',
      line = 2.75)
mtext(side = 1, line = 2.5, text = 'Assimilation period')
dev.off()





pdf('../output/main_text_figures/particle_retention/forecasted_infections_m1.pdf',
    height = 5.5, width = 6.5)
layout(matrix(c(1,1,1,1,1,2:16), 4, 5, byrow = T))
par(mar = c(2,1,1,0), oma = c(4.5,4,1,1))
barplot(sapply(unique(df_dyn$week), function(ww) 
  sum(df_dyn$cases_new[which(df_dyn$week == ww)]))[week_first_case:max(df_dyn$week)],
  space = 0, col = 'grey90',
  border = 'grey', las = 1)
mtext(side = 2, text = 'Reported cases', line = 3.5)
abline(v = when_to_assimilate - week_first_case, col = cols, lwd = 2)
par(mar = c(1,1,1,0))
for(tt in 1:length(when_to_assimilate)){
  
  dens = log10(total_infs_per_run_per_time[tt,])
  dens[which(is.infinite(dens))] = 0
  
  hist(dens,
       col = adjustcolor(cols[tt], alpha.f = 0.1),
       breaks = 100, main = '', xlab = '', ylab = '',
       las = 1, 
       ylim = c(0, 6e3),
       xlim = c(0, 8),
       border = adjustcolor(cols[tt], alpha.f = 0.6), 
       xaxt = 'n',
       yaxt = 'n')
  
  if(tt %in% 11:15){
    axis(side = 1, at = c(0:8), 
         labels = c(0, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8),
         las = 2)
  }else{
    axis(side = 1, at = c(0:8), 
         labels = NA,
         las = 2)
  }
  
  if(tt %in% c(1,6,11)){
    axis(side = 2, las = 1)
  }else{
    axis(side = 2, las = 1, labels = NA)
  }
  
  
  if(tt == 6){
    mtext(side = 2, line = 3.5, text = 'Count')
  }
}
mtext(side = 1, text = 'Forecasted infections (log10)', outer = T, line = 3)
dev.off()





