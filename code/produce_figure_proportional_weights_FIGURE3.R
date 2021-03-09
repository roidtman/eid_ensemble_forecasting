#=============================================================================#
# Author: Rachel Oidtman

# FIGURE 3
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
# visualize ensemble weights
#=============================================================================#


# weight of mobility through time
# order: cell phone mobility, gravity, radiation
mobility_weights_through_time = matrix(NA, ncol = length(when_to_assimilate),
                                       nrow = 4)

R_weights_through_time = matrix(NA, ncol = length(when_to_assimilate),
                                nrow = 2)

introduction_weights_through_time = matrix(NA, ncol = length(when_to_assimilate),
                                           nrow = 2)

for(tt in 1:length(when_to_assimilate)){
  mobility_weights_through_time[1, tt] = sum(ensemble_wts_through_time[1:4,tt]) / sum(ensemble_wts_through_time[,tt])
  mobility_weights_through_time[2, tt] = sum(ensemble_wts_through_time[9:12,tt]) / sum(ensemble_wts_through_time[,tt])
  mobility_weights_through_time[3, tt] = sum(ensemble_wts_through_time[13:16,tt]) / sum(ensemble_wts_through_time[,tt])
  mobility_weights_through_time[4, tt] = sum(ensemble_wts_through_time[5:8,tt]) / sum(ensemble_wts_through_time[,tt])
  
  R_weights_through_time[1, tt] = sum(ensemble_wts_through_time[c(1,2,5,6,9,10,13,14), tt]) / sum(ensemble_wts_through_time[, tt])
  R_weights_through_time[2, tt] = sum(ensemble_wts_through_time[c(3,4,7,8,11,12,15,16), tt]) / sum(ensemble_wts_through_time[, tt])
  
  introduction_weights_through_time[1, tt] = sum(ensemble_wts_through_time[c(1,3,5,7,9,11,13,15), tt]) / sum(ensemble_wts_through_time[, tt])
  introduction_weights_through_time[2, tt] = sum(ensemble_wts_through_time[c(2,4,6,8,10,12,14,16), tt]) / sum(ensemble_wts_through_time[, tt])
}



# cumulative number of locations reporting cases
number_of_depts_with_week_first_case = rep(0, length(when_to_assimilate))

for(dd in 1:length(dept_names)){
  
  tmp_wk_first_case = df_stat$wk_first_case[dd] - week_first_case + 1
  tmp = sapply(1:nrow(wks_assimilated_at_when_to_assimilate), function(ff) 
    tmp_wk_first_case %in% wks_assimilated_at_when_to_assimilate[ff,])
  
  number_of_depts_with_week_first_case[which(tmp)] = number_of_depts_with_week_first_case[which(tmp)] + 1
}

pdf('../output/main_text_figures/prop_ensemble_weight_WITH_data.pdf', height = 6, width = 5)
layout(matrix(c(
  1,1,1,1,5,5,
  2,2,2,2,6,6,
  2,2,2,2,6,6,
  3,3,3,3,7,7,
  3,3,3,3,7,7,
  4,4,4,4,8,8,
  4,4,4,4,8,8), byrow = T, nrow = 7))

par(mar = c(3,0,1,0), oma = c(2,4.5,2.5,0))
cols = c('mistyrose', 'darkolivegreen3', 'lightsalmon','lightsteelblue1',
         'thistle', 'khaki1',
         'lightslategrey', 'snow2')

par(mar = c(2,0,0,0))

national_cases = sapply(week_first_case : max(df_dyn$week), function(ff) sum(df_dept$cases_new[which(df_dept$week == ff)]))

plot(national_cases, type = 'l', lwd = 3, col = adjustcolor('grey30', alpha.f = 0.6), bty = 'n', las = 1, xaxt = 'n', xlab = '', ylab = '',
     yaxt = 'n', ylim = c(0, 5800))
axis(side = 2, at = c(0, 1000, 2000, 3000, 4000, 5000, 6000), las = 1)
points(national_cases, pch = 16, cex = 1.5, col = adjustcolor('grey30', alpha.f = 0.6))
axis(side = 1, at = c((when_to_assimilate - week_first_case), 60), las = 2)
mtext(side = 2, line = 3.2, text = 'National incidence', cex = 0.8)


par(mar = c(3,0,2,0))
barplot(mobility_weights_through_time, col = cols[1:4], border = 'white', las = 1,
        space = rep(0, ncol(mobility_weights_through_time))
)
axis(side = 1, at = 0:15, labels = NA)
axis(side = 1, at = seq(0.5, 14.5, by = 1), labels = when_to_assimilate - week_first_case, tick = F, las = 2)
mtext(side = 3, text = 'Assumption: Human mobility', cex = 0.8)


barplot(R_weights_through_time, col = cols[5:6], border = 'white', las = 1,
        space = rep(0, ncol(mobility_weights_through_time))
)
axis(side = 1, at = 0:15, labels = NA)
axis(side = 1, at = seq(0.5, 14.5, by = 1), labels = when_to_assimilate - week_first_case, tick = F, las= 2)
mtext(side = 3, text = 'Assumption: Transmission potential', cex = 0.8)
mtext(side = 2, text = 'Proportion of ensemble weight', cex = 0.8, line = 3.2)


barplot(introduction_weights_through_time, col = cols[7:8], border = 'white', las = 1,
        space = rep(0, ncol(mobility_weights_through_time))
)
axis(side = 1, at = 0:15, labels = NA)
axis(side = 1, at = seq(0.5, 14.5, by = 1), labels = when_to_assimilate - week_first_case, tick = F, las = 2)
mtext(side = 3, text = 'Assumption: Number of introductions', cex = 0.8)

mtext(side = 1, line = 2.5, text = 'Weeks since 1st reported case in Colombia', cex = 0.8)

# ADD LEGENDS
plot(-100, -100, xlim = c(-1, 1), ylim = c(0, 1), xaxt = 'n', yaxt = 'n', bty = 'n')


plot(-100, -100, xlim = c(-1, 1), ylim = c(0, 1), xaxt = 'n', yaxt = 'n', bty = 'n')
legend('topleft', legend = c('CDR-informed', 
                             'Gravity', 
                             'Radiation',
                             'No movement'), bty = 'n', 
       fill = cols[1:4], cex = 1.1, border = cols[1:4])

plot(-100, -100, xlim = c(-1, 1), ylim = c(0, 1), xaxt = 'n', yaxt = 'n', bty = 'n')
legend('topleft', legend = c('Static R', 
                             'Dynamic R'), bty = 'n', 
       fill = cols[5:6], cex = 1.1, border = cols[5:6])

plot(-100, -100, xlim = c(-2, 2), ylim = c(-1, 1), xaxt = 'n', yaxt = 'n', bty = 'n')
legend('topleft', legend = c('1 introduction', 
                             '2 introductions'), bty = 'n', 
       fill = cols[7:8], cex = 1.1, border = cols[7:8])


dev.off()
