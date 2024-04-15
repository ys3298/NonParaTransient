
setwd('/storage/work/yqs5519/SurvivalPlot/application')
source('/storage/work/yqs5519/SurvivalPlot/algorithm/calibrateAUC.R')
source('/storage/work/yqs5519/SurvivalPlot/algorithm/CI_functions.R')
source('/storage/work/yqs5519/SurvivalPlot/Mantel_Byar/helper_functions.R')

## transform the data
convert_chimerism_totalcell = function(d1)
{
  n = dim(d1)[1]
  dd = NULL
  for(i in 1:n)
  {
    one_row = d1[i,]
    #################################  
    if(d1$chimerism_drop_group[i]==2) 
    {
      #subject, time1, time2, DFS_STATUS, chimerism_drop_group_dynamic
      temp1 = c(one_row$subject, 0, one_row$chimerism_first_drop_day, 0, 0)
      dd = rbind(dd, temp1)
      
      temp2 = c(one_row$subject, one_row$chimerism_first_drop_day, one_row$chimerism_second_drop_day,  0, 1)
      dd = rbind(dd, temp2)
      
      temp3 = c(one_row$subject, one_row$chimerism_second_drop_day, one_row$relapse_free_days, one_row$relapse_free_status, 2)
      dd = rbind(dd, temp3)
      
    }
    
    if(d1$chimerism_drop_group[i]==1) 
    {
      temp1 = c(one_row$subject, 0, one_row$chimerism_first_drop_day, 0, 0)
      dd = rbind(dd, temp1)
      
      temp2 = c(one_row$subject, one_row$chimerism_first_drop_day, one_row$relapse_free_days,  one_row$relapse_free_status, 1)
      dd = rbind(dd, temp2)
    }
    ####################################
    
    if(d1$chimerism_drop_group[i]==0) 
    {
      temp2 = c(one_row$subject, 0, one_row$relapse_free_days, one_row$relapse_free_status, 0)
      dd = rbind(dd, temp2)
    }
    
  }
  
  dd = as.data.frame(dd)
  names(dd) = c("subject", "time1", "time2", "relapse_free_status", "chimerism_drop_group_dynamic")
  dd = dd %>% mutate(chimerism_drop_group_dynamic = as.factor(chimerism_drop_group_dynamic))
  return(dd)
}

## import data
raw_data = read.csv('./2022-05-20 IMC final.csv')
d1 = raw_data %>% 
  mutate(Total_Drop_group = as.factor(Total_Drop_group)) %>% 
  mutate(relapse_free_status = ifelse(DFS_STATUS == 2, 1, DFS_STATUS),
         relapse_free_days = DFS_days) %>% 
  dplyr::select(Liao.s.., overall_survival_days, Censor_overall_survival, 
                Total_Drop_group, Total_First_drop_days, Total_Second_drop_days,
                relapse_free_status, relapse_free_days)
colnames(d1)[1] = "subject"
colnames(d1)[4] = "chimerism_drop_group"
colnames(d1)[5] = "chimerism_first_drop_day"
colnames(d1)[6] = "chimerism_second_drop_day"

total_cell_time = convert_chimerism_totalcell(d1) %>% na.omit() 
colnames(total_cell_time) = c('id', 'start_td', 'stop_td', 'endpoint_td', "covariate_td")
TempTD = total_cell_time
out1 = Mantel.Byar.test(0, 1, Group = "covariate_td", Event = TempTD$endpoint_td, 
                        StartTime = TempTD$start_td, StopTime = TempTD$stop_td, 
                        method = c("SAS"), plot=0, landmark=0)
out1

out2 = Mantel.Byar.test(0, 2, Group = "covariate_td", Event = TempTD$endpoint_td, 
                        StartTime = TempTD$start_td, StopTime = TempTD$stop_td, 
                        method = c("SAS"), plot=0, landmark=0)
out2

out3 = Mantel.Byar.test(1, 2, Group = "covariate_td", Event = TempTD$endpoint_td, 
                        StartTime = TempTD$start_td, StopTime = TempTD$stop_td, 
                        method = c("SAS"), plot=0, landmark=0)
out3
