#=============================================================================#
# Author: Rachel Oidtman

# Process epidemiological data 
#=============================================================================#


#=============================================================================#
# process and save epi data for easy loading
#=============================================================================#

df_full = read.csv('../data_bricks/epi_data.csv', stringsAsFactors = F)
muni_shape = read.csv('../data/shapefiles_muni_id.csv', stringsAsFactors = F)
dept_names = as.character(unique(muni_shape$dept)[-length(unique(muni_shape$dept))])


#=============================================================================#
# creating department level data frame
#=============================================================================#

df_dept = data.frame(dept_name = rep(dept_names, max(df_full$week)),
                     week = rep(1:max(df_full$week), each = length(dept_names)),
                     pop = rep(0, length(dept_names) * max(df_full$week)),
                     R0 = rep(0, length(dept_names) * max(df_full$week)),
                     cases_new = rep(0, length(dept_names) * max(df_full$week)),
                     week_first_case = rep(NA, length(dept_names) * max(df_full$week)),
                     num_munis = rep(0, length(dept_names) * max(df_full$week)))

for(mm in unique(df_full$muni)){
  
  tmp_id = muni_shape[which(muni_shape$mpios == mm),]
  
  df_dept$cases_new[which(df_dept$dept_name == tmp_id$dept)] =  df_dept$cases[which(df_dept$dept_name == tmp_id$dept)] + df_full$cases_new[which(df_full$muni == mm)]
  df_dept$pop[which(df_dept$dept_name == tmp_id$dept)] =  df_dept$pop[which(df_dept$dept_name == tmp_id$dept)] + df_full$pop[which(df_full$muni == mm)]
  df_dept$R0[which(df_dept$dept_name == tmp_id$dept)] =  df_dept$R0[which(df_dept$dept_name == tmp_id$dept)] + df_full$R0[which(df_full$muni == mm)]
  df_dept$num_munis[which(df_dept$dept_name == tmp_id$dept)] = df_dept$num_munis[which(df_dept$dept_name == tmp_id$dept)] + 1
  
}


df_dept$R0 = df_dept$R0 / df_dept$num_munis

for(dd in 1:length(unique(dept_names))){
  df_tmp = df_dept[which(df_dept$dept_name == dept_names[dd]),]
  wk_min_tmp = min(df_tmp$week[which(df_tmp$cases != 0)])
  
  df_dept$week_first_case[which(df_dept$dept_name == dept_names[dd])] = wk_min_tmp
}

df_dept$dept_name = as.character(df_dept$dept_name)


#=============================================================================#
# important global variables
#=============================================================================#

loc_ini = df_dept$dept_name[which(df_dept$week_first_case == min(df_dept$week_first_case, na.rm =T))][1]
week_first_case = min(df_dept$week_first_case, na.rm = T)
week_first_case_min = 50


#=============================================================================#
# static and dynamic dataframes
#=============================================================================#

df_dyn = df_dept[,c(1:2,5)]
df_dyn$inf_new = rep(0, nrow(df_dyn))
df_dyn$inf_cum = rep(0, nrow(df_dyn))
df_dyn$susc_recon = rep(NA, nrow(df_dyn))
df_dyn$inf_old_move = rep(NA, nrow(df_dyn))

df_stat = data.frame(dept_name = unique(df_dept$dept_name),
                     wk_first_case = df_dept$week_first_case[1:length(unique(df_dept$dept_name))],
                     pop = df_dept$pop[1:length(unique(df_dept$dept_name))],
                     R0 = df_dept$R0[1:length(unique(df_dept$dept_name))])


#=============================================================================#
# save output
#=============================================================================#

save(loc_ini, week_first_case, week_first_case_min, 
     df_dept, dept_names,
     df_dyn, df_stat,
     file = '../data/processed/processed_epi_data.RData')


