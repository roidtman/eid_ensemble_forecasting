#=============================================================================#
# Author: Rachel Oidtman

# FIGURE 1
# DATA
#=============================================================================#


#=============================================================================#
# load in processed epi data
#=============================================================================#

require(stats)
require(viridis)

load('../data/processed/processed_epi_data.RData')
load('../data/processed/temp_mosq_data_frame.RData')

dept_models = c(1:21, 23:32)

time_btwn_assim = 4
when_to_assimilate = seq(week_first_case + time_btwn_assim, max(df_dyn$week), by = time_btwn_assim)
when_to_assimilate = c(week_first_case, when_to_assimilate)

load('../data/mobility_models/processed_gravity_model.RData')
load('../data/mobility_models/processed_radiation_model.RData')
connectivity_CDR_normalized = read.csv('../data/mobility_models/connectivity_normalized_department_min_value_0.5.csv')
connectivity_CDR_normalized = as.matrix(connectivity_CDR_normalized)


dept_names_to_plot = c('La Guajira', 'Magdalena', 'Atlántico', 'Cesar', 'Bolívar',
                       'Sucre', 'Córdoba', 'Norte de Santander', 'Antioquia', 
                       'Chocó', 'Santander', 'Arauca', 'Boyacá', 'Vichada', 
                       'Casanare', 'Cudinamarca', 'Caldas', 'Risaralda', 'Tolima', 
                       'Valle del Cauca', 'Meta', 'Quindío', 'Guainía', 'Huila', 
                       'Cauca', 'Caquetá', 'Guaviare', 'Nariño', 'Vaupés', 
                       'Putumayo', 'Amazonas')

#=============================================================================#
# necessary functions
#=============================================================================#


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


#=============================================================================#
# plots
#=============================================================================#

mat_to_plot = matrix(NA, nrow = length(dept_models), ncol = length(week_first_case : max(df_dyn$week)))

l = 1
for(dd in dept_models){
  mat_to_plot[l,] = df_dept$cases_new[which(df_dept$week %in% week_first_case : max(df_dyn$week) & df_dept$dept_name == dept_names[dd])]
  l = l + 1
}

cum_cases_by_dept = sapply(1:length(dept_names), function(ff) sum(df_dept$cases_new[which(df_dept$dept_name == dept_names[ff])][week_first_case : max(df_dyn$week)]))
cum_cases_by_dept = cum_cases_by_dept[-22]
ys_to_plot = order(cum_cases_by_dept)


pdf('../output/main_text_figures/data_figure.pdf', height = 5.5, width = 8.5)
# layout(matrix(c(1,1,1,1,1,1,1,2,3,3,3,3,3,3,
#                 1,1,1,1,1,1,1,2,3,3,3,3,3,3,
#                 1,1,1,1,1,1,1,2,4,4,4,4,4,4,
#                 1,1,1,1,1,1,1,2,4,4,4,4,4,4,
#                 1,1,1,1,1,1,1,2,5,5,5,5,5,5,
#                 1,1,1,1,1,1,1,2,5,5,5,5,5,5,
#                 1,1,1,1,1,1,1,2,6,6,7,7,8,8,
#                 1,1,1,1,1,1,1,2,6,6,7,7,8,8,
#                 1,1,1,1,1,1,1,2,6,6,7,7,8,8), byrow = T, nrow = 9))

# layout(matrix(c(1,1,1,1,1,1,1,2,3,3,3,3,3,3,
#               1,1,1,1,1,1,1,2,3,3,3,3,3,3,
#               1,1,1,1,1,1,1,2,4,4,4,4,4,4,
#               1,1,1,1,1,1,1,2,4,4,4,4,4,4,
#               1,1,1,1,1,1,1,2,5,5,6,6,7,7,
#               1,1,1,1,1,1,1,2,5,5,6,6,7,7), byrow = T, nrow = 6))

layout(matrix(c(1,1,1,1,1,1,1,2,2,3,3,3,3,3,3,3,3,3,
                1,1,1,1,1,1,1,2,2,3,3,3,3,3,3,3,3,3,
                1,1,1,1,1,1,1,2,2,4,4,4,4,4,4,4,4,4,
                1,1,1,1,1,1,1,2,2,4,4,4,4,4,4,4,4,4,
                1,1,1,1,1,1,1,2,2,5,5,5,6,6,6,7,7,7,
                1,1,1,1,1,1,1,2,2,5,5,5,6,6,6,7,7,7), byrow = T, nrow = 6))

par(mar = c(2,2,0,1), oma = c(3,9,2,1))

# cols = colorRampPalette(c('azure', 'darkslategray4'))(50)

cols = c(colorRampPalette(c('honeydew2', 'lightsteelblue'))(10),
         colorRampPalette(c('lightsteelblue', 'cadetblue3'))(10),
         colorRampPalette(c('cadetblue3', 'midnightblue'))(10),
         colorRampPalette(c('midnightblue', 'magenta4'))(10))

mat_to_plot_rot = apply(mat_to_plot, 2, rev)
# image(t(mat_to_plot_rot))

image(log10(t(mat_to_plot_rot)), xaxt = 'n', yaxt = 'n', col = cols)
box()
axis(side = 2, at = seq(0, 1, length.out = 31), labels = rev(dept_names_to_plot), las = 1)
axis(side = 1, at = seq(0, 1, length.out = length(when_to_assimilate)), labels = when_to_assimilate - week_first_case)
mtext(side = 1, line = 2.25, cex = 0.8, text = 'Weeks since 1st reported case in Colombia')
mtext(side = 2, line = 9, cex = 0.8, text = 'Departments (south to north)')

mtext(side = 3, line = 0, text = 'Weekly', at = 1.16, cex = 0.8)
mtext(side = 3, line = -1, text = 'incidence', at = 1.2, cex = 0.8)

par(mar = c(0.5,0,0,3.5))
color.bar(c(cols, 'white'), min = -1)
axis(side = 4, at = seq(-1, 1, length.out = 1590)[c(1, round(log10(c(25, 50, 100, 300, 500, 750, 1000, 1590)) * (1590 / 3.201397)))], 
     labels = c(1, 25, 50, 100, 300, 500, 750, 1000, 1600), las = 1)

par(mar = c(3,5,0,0))

# national_cases = sapply(week_first_case : max(df_dyn$week), function(ff) sum(df_dept$cases_new[which(df_dept$week == ff)]))
# 
# plot(national_cases, type = 'l', lwd = 3, col = adjustcolor('grey30', alpha.f = 0.6), bty = 'n', las = 1, xaxt = 'n', xlab = '', ylab = '')
# points(national_cases, pch = 16, cex = 1.5, col = adjustcolor('grey30', alpha.f = 0.6))
# axis(side = 1, at = c((when_to_assimilate - week_first_case), 60), las = 2)
# mtext(side = 2, line = 3, text = 'National incidence', cex = 0.8)


mean_temperature_across_depts = sapply(week_first_case : max(df_dyn$week), function(ff) mean(df_dept_temp_mosq$temp_mean[which(df_dept_temp_mosq$week == ff)]))

plot(mean_temperature_across_depts, type = 'l', lwd = 3, col = adjustcolor('grey30', alpha.f = 0), 
     bty = 'n', las = 1, xaxt = 'n', xlab = '', ylab = '',
     ylim = c(min(df_dept_temp_mosq$temp_mean), max(df_dept_temp_mosq$temp_mean)))
axis(side = 1, at = c((when_to_assimilate - week_first_case), 60), las = 2)
for(dd in dept_models){
  lines(df_dept_temp_mosq$temp_mean[which(df_dept_temp_mosq$dept_name == dept_names[dd])][week_first_case:max(df_dyn$week)], col = adjustcolor('grey30', alpha.f = 0.3))
}
lines(mean_temperature_across_depts, lwd = 3, col = adjustcolor('grey30', alpha.f = 0.6))
points(mean_temperature_across_depts, pch = 16, cex = 1.5, col = adjustcolor('grey30', alpha.f = 0.6))
mtext(side = 2, line = 3, text = 'Temperature (c)', cex = 0.8)


mean_mosq_across_depts = sapply(week_first_case : max(df_dyn$week), function(ff) mean(df_dept_temp_mosq$mosq_op[which(df_dept_temp_mosq$week == ff)]))

plot(mean_mosq_across_depts, type = 'l', lwd = 3, col = adjustcolor('grey30', alpha.f = 0), 
     bty = 'n', las = 1, xaxt = 'n', xlab = '', ylab = '',
     ylim = c(min(df_dept_temp_mosq$mosq_op), max(df_dept_temp_mosq$mosq_op)))
axis(side = 1, at = c((when_to_assimilate - week_first_case), 60), las = 2)
for(dd in dept_models){
  lines(df_dept_temp_mosq$mosq_op[which(df_dept_temp_mosq$dept_name == dept_names[dd])][week_first_case:max(df_dyn$week)], 
        col = adjustcolor('grey30', alpha.f = 0.3))
}
lines(mean_mosq_across_depts, lwd = 3, col = adjustcolor('grey30', alpha.f = 0.6))
points(mean_mosq_across_depts, pch = 16, cex = 1.5, col = adjustcolor('grey30', alpha.f = 0.6))
mtext(side = 1, line = 2.5, text = 'Weeks since 1st reported case in Colombia', cex = 0.8)
mtext(side = 2, line = 3, text = 'Pr(mosq occur)', cex = 0.8)


cols_mobility = magma(100)
mobility_models_to_plot = array(NA, dim = c(length(dept_models), length(dept_models), 3))

mobility_models_to_plot[,,1] = log10(connectivity_CDR_normalized[-22, -22])
mobility_models_to_plot[,,2] = log10(connectivity_gravity_normalized[-22, -22])
mobility_models_to_plot[,,3] = log10(connectivity_radiation_normalized[-22, -22])

mobility_models_to_plot[which(is.infinite(mobility_models_to_plot))] = NA

par(mar = c(2,1,3,0))
mat_to_plot_tmp = apply(connectivity_CDR_normalized[-22,-22], 2, rev)
image(log10(t(mat_to_plot_tmp)), zlim = c(min(mobility_models_to_plot, na.rm = T),
                                                   max(mobility_models_to_plot, na.rm = T)), col = cols_mobility,
      xaxt = 'n', yaxt = 'n')
box()
mtext(side = 2, line = 0, text = 'Origin', cex = 0.8)
mtext(side = 1, line = 0, text = 'Destination', cex = 0.8)
mtext(side = 3, line = 0, text = 'CDR-informed', cex = 0.8)

mat_to_plot_tmp = apply(connectivity_gravity_normalized[-22,-22], 2, rev)
image(log10(t(mat_to_plot_tmp)), zlim = c(min(mobility_models_to_plot, na.rm = T),
                                                       max(mobility_models_to_plot, na.rm = T)), col = cols_mobility,
      xaxt = 'n', yaxt = 'n')
box()
mtext(side = 1, line = 0, text = 'Destination', cex = 0.8)
mtext(side = 3, line = 0, text = 'Gravity', cex = 0.8)

mat_to_plot_tmp = apply(connectivity_radiation_normalized[-22,-22], 2, rev)
image(log10(t(mat_to_plot_tmp)), zlim = c(min(mobility_models_to_plot, na.rm = T),
                                                         max(mobility_models_to_plot, na.rm = T)), col = cols_mobility,
      xaxt = 'n', yaxt = 'n')
box()
mtext(side = 1, line = 0, text = 'Destination', cex = 0.8)
mtext(side = 3, line = 0, text = 'Radiation', cex = 0.8)

dev.off()









#=============================================================================#
# mobility plots
#=============================================================================#


pdf('~/Dropbox/PhD_general/defense/figs/movement_assumption.pdf',
    height = 1.9, width =7)
layout(matrix(1:4, nrow = 1))
cols_mobility = magma(100)
mobility_models_to_plot = array(NA, dim = c(length(dept_models), length(dept_models), 4))

mobility_models_to_plot[,,1] = log10(connectivity_CDR_normalized[-22, -22])
mobility_models_to_plot[,,2] = log10(connectivity_gravity_normalized[-22, -22])
mobility_models_to_plot[,,3] = log10(connectivity_radiation_normalized[-22, -22])
mobility_models_to_plot[,,4] = log10(diag(rep(1, 31)))

mobility_models_to_plot[which(is.infinite(mobility_models_to_plot))] = NA

par(mar = c(2,1,3,2))
mat_to_plot_tmp = apply(connectivity_CDR_normalized[-22,-22], 2, rev)
image(log10(t(mat_to_plot_tmp)), zlim = c(min(mobility_models_to_plot, na.rm = T),
                                          max(mobility_models_to_plot, na.rm = T)), col = cols_mobility,
      xaxt = 'n', yaxt = 'n')
box()
mtext(side = 2, line = 0, text = 'Origin', cex = 0.8)
mtext(side = 1, line = 0, text = 'Destination', cex = 0.8)
mtext(side = 3, line = 0, text = 'Cell phone mobility', cex = 0.8)

mat_to_plot_tmp = apply(connectivity_gravity_normalized[-22,-22], 2, rev)
image(log10(t(mat_to_plot_tmp)), zlim = c(min(mobility_models_to_plot, na.rm = T),
                                          max(mobility_models_to_plot, na.rm = T)), col = cols_mobility,
      xaxt = 'n', yaxt = 'n')
box()
mtext(side = 1, line = 0, text = 'Destination', cex = 0.8)
mtext(side = 3, line = 0, text = 'Gravity', cex = 0.8)

mat_to_plot_tmp = apply(connectivity_radiation_normalized[-22,-22], 2, rev)
image(log10(t(mat_to_plot_tmp)), zlim = c(min(mobility_models_to_plot, na.rm = T),
                                          max(mobility_models_to_plot, na.rm = T)), col = cols_mobility,
      xaxt = 'n', yaxt = 'n')
box()
mtext(side = 1, line = 0, text = 'Destination', cex = 0.8)
mtext(side = 3, line = 0, text = 'Radiation', cex = 0.8)


mat_to_plot_tmp = apply(log10(diag(rep(1, 31))), 2, rev)
image((t(mat_to_plot_tmp)), zlim = c(min(mobility_models_to_plot, na.rm = T),
                                          max(mobility_models_to_plot, na.rm = T)), col = cols_mobility,
      xaxt = 'n', yaxt = 'n')
box()
mtext(side = 1, line = 0, text = 'Destination', cex = 0.8)
mtext(side = 3, line = 0, text = 'No inter-dept. movement', cex = 0.8)

dev.off()
