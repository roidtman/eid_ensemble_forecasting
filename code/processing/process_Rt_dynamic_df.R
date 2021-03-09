#=============================================================================#
# Author: Rachel Oidtman

# Process file with temperature and mosquito occurrence probability at the department level
#=============================================================================#


#=============================================================================#
# load in processed epi data
#=============================================================================#

library(viridis)
library(lubridate)
library(dplyr)
library(sn)


load('../data/processed/processed_epi_data.RData')
date_first_case = ymd('2015-08-09')

# load in temp and mosqsuito occurrence probabiity data
temp_mean = read.csv('../data/dept_tmean_weekly_weighted.csv', stringsAsFactors=F)
aedes_op = read.csv('../data/aegypti_pop_monthly_col_dept.csv', stringsAsFactors= F)


#=============================================================================#
# processing data
#=============================================================================#

# find weekly mean temperature per department
temp_14 = select(temp_mean, ends_with('14'))
temp_15 = select(temp_mean, ends_with('15'))
temp_16 = select(temp_mean, ends_with('16'))


temp_mean_across_years = matrix(NA, nrow = nrow(temp_14), ncol = ncol(temp_14))


for(ii in 1:ncol(temp_14)){
  
  if(ii <= ncol(temp_16)){
    temp_mean_across_years[,ii] = sapply(1:nrow(temp_14), function(ff) mean(c(temp_14[ff,ii], temp_15[ff,ii], temp_16[ff,ii])))
  }else{
    temp_mean_across_years[,ii] = sapply(1:nrow(temp_14), function(ff) mean(c(temp_14[ff,ii], temp_15[ff,ii])))
  }
}; rm(ii)


# updating department names in temp_mean so that they match those in the df_dept
temp_mean$DEPARTMENT[1] = 'SANTAFE DE BOGOTA D.C'
temp_mean$DEPARTMENT[2] = 'BOLIVAR'
temp_mean$DEPARTMENT[3] = 'BOYACA'
temp_mean$DEPARTMENT[5] = 'CAQUETA'
temp_mean$DEPARTMENT[8] = 'CORDOBA'
temp_mean$DEPARTMENT[10] = 'CHOCO'


# adding temperature data to department dataframe
l = 1
for(tt in 1:max(df_dept$week)){
  if(tt == 53 | tt == 105){
    l = 1
  }
  for(dd in 1:length(unique(df_dept$dept_name))){
    df_dept$temp_mean[which(df_dept$dept_name == unique(df_dept$dept_name)[dd] & df_dept$week == tt)] = temp_mean_across_years[which(temp_mean$DEPARTMENT == unique(df_dept$dept_name)[dd]), l]
  }
  l = l + 1
}
rm(tt, l, dd)


# updating departments in aedes_op so they're the same as df_dept
aedes_op$NOM_DEPART[13] = 'BOYACA'
aedes_op$NOM_DEPART[14] = 'BOLIVAR'
aedes_op$NOM_DEPART[15] = 'SANTAFE DE BOGOTA D.C'
aedes_op$NOM_DEPART[16] = 'CORDOBA'
aedes_op$NOM_DEPART[19] = 'CAQUETA'
aedes_op$NOM_DEPART[24] = 'CHOCO'


# adding months
wks_14 = c(4,4,5,4,4,5,4,4,5,4,4,5)
wks_15 = c(4,4,5,4,4,5,4,4,5,4,4,5)
wks_16 = c(4,5,4,4,5,4,4,5,4)

months = c(rep(1:12, wks_14), rep(1:12, wks_15), rep(1:9, wks_16))


# adding mosquito occurrence probability to department data frame
for(tt in 1:max(df_dept$week)){
  for(dd in 1:length(unique(df_dept$dept_name))){
    df_dept$mosq_op[which(df_dept$dept_name == unique(df_dept$dept_name)[dd] & df_dept$week == tt)] = 
      aedes_op[which(aedes_op$NOM_DEPART == unique(df_dept$dept_name)[dd]), 2 + months[tt]]
  }
}

# converting mosquito occurrence probability to mosquito abundance
occur_abund = function(op){
  -log(1 - op)
}

df_dept$mosq_abund = rep(NA, nrow(df_dept))
df_dept$mosq_abund = occur_abund(df_dept$mosq_op)


head(df_dept)

df_dept_temp_mosq = df_dept

save(df_dept_temp_mosq, file = '../data/processed/temp_mosq_data_frame.RData')




