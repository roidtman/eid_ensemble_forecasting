#=============================================================================#
# Author: Rachel Oidtman

# Supplementary Figure 
# Comparing mosquito dynamics and temperature to reported cases 
#=============================================================================#


#=============================================================================#
# load in processed epi data
#=============================================================================#
library(lubridate)

load('../data/processed/processed_epi_data.RData')
load('../data/processed/temp_mosq_data_frame.RData')
date_first_case = ymd('2015-08-09')


dept_names_to_plot = c('La Guajira', 'Magdalena', 'Atlántico', 'Cesar', 'Bolívar',
                       'Sucre', 'Córdoba', 'Norte de Santander', 'Antioquia', 
                       'Chocó', 'Santander', 'Arauca', 'Boyacá', 'Vichada', 
                       'Casanare', 'Cudinamarca', 'Caldas', 'Risaralda', 'Tolima', 
                       'Valle del Cauca', 'Meta', 'Bogota', 'Quindío', 'Guainía', 'Huila', 
                       'Cauca', 'Caquetá', 'Guaviare', 'Nariño', 'Vaupés', 
                       'Putumayo', 'Amazonas')

#=============================================================================#
# load in processed epi data
#=============================================================================#

# mosquito occurrence probability versus cases
pdf('../output/main_text_figures/supplement/mosq_vs_cases.pdf',
    height = 8, width = 7)
layout(matrix(1:32, 8,4))
par(mar = c(1.5,2,2,1), oma = c(3,3,1,0.5))
for(dd in dept_names){
  plot(df_dept$cases_new[which(df_dept$week %in% 84:max(df_dept$week) & df_dept$dept_name == dd)], type = 'h', lwd = 3,
       col = adjustcolor('navy', alpha.f = 0.6), las = 1)
  mtext(side = 3, text = dept_names_to_plot[which(dept_names == dd)], cex = 0.75)
  
  par(new = T)
  plot(df_dept_temp_mosq$mosq_op[which(df_dept$week %in% 84:max(df_dept$week) & df_dept$dept_name == dd)], type = 'l',
       xaxt = 'n', yaxt = 'n', bty = 'n', col = 'seagreen3', lwd = 3)
}
mtext(side = 1, text = 'Weeks since 1st reported case in Colombia', outer = T, line = 1)
mtext(side = 2, text = 'New weekly cases', outer = T, line = 1)
dev.off()


# temperature versus cases
pdf('../output/main_text_figures/supplement/temp_vs_cases.pdf',
    height = 8, width = 7)
layout(matrix(1:32, 8,4))
par(mar = c(1.5,2,2,1), oma = c(3,3,1,0.5))
for(dd in dept_names){
  plot(df_dept$cases_new[which(df_dept$week %in% 84:max(df_dept$week) & df_dept$dept_name == dd)], type = 'h', lwd = 3,
       col = adjustcolor('navy', alpha.f = 0.6), las = 1)
  mtext(side = 3, text = dept_names_to_plot[which(dept_names == dd)], cex = 0.75)

  par(new = T)
  plot(df_dept_temp_mosq$temp_mean[which(df_dept$week %in% 84:max(df_dept$week) & df_dept$dept_name == dd)], type = 'l',
       xaxt = 'n', yaxt = 'n', bty = 'n', col = 'firebrick3', lwd = 3)
}
mtext(side = 1, text = 'Weeks since 1st reported case in Colombia', outer = T, line = 1)
mtext(side = 2, text = 'New weekly cases', outer = T, line = 1)
dev.off()


