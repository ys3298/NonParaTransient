source('/storage/work/yqs5519/SurvivalPlot/Mantel_Byar/helper_functions.R')
source('/storage/work/yqs5519/SurvivalPlot/algorithm/generate_data.R')
library(foreach, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
library(doParallel, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
library(tidyverse, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')

# file_num = readr::parse_number(list.files('/storage/work/yqs5519/SurvivalPlot/Mantel_Byar', 'MB_test_results_power'))
# location_index = ifelse((is_empty(file_num)), 1,  max(file_num, na.rm = T) + 1) ## set location to save results
# location = paste('/storage/work/yqs5519/SurvivalPlot/Mantel_Byar/MB_test_results_power', location_index, '/', sep = '')
# dir.create(location)

## settings
N_ite = 10
corr = 1

Lambda1 = 1
Rho1 = 1
Lambda2 = 1
Rho2 = 1
# scenario = 'test power (time dependent transition)'

n_sample = 200
scenario = 3
Trans_time = 0.1
nocensor = T
time_end = 15
# F = have censor
if (scenario == 1) {
  # scenario 1
  Beta = c(1,1,1)
  Gamma =  c(1.5,1.5,NA)
  CensorC = c(2.7,5,3)
} else if (scenario == 2) {
  # scenario 2
  Beta = c(1,1,1)
  Gamma =  c(1,1.5,NA)
  CensorC = c(2,4.5,4)
} else if (scenario == 3) {
  Beta = c(0.5,1,1.5)
  Gamma = c(1,1,NA)
  CensorC = c(1,5,5.5)
} else if (scenario == 4) {
  Beta = c(0.5,1,1.5)
  Gamma = c(0.5,1.5,NA)
  CensorC = c(1,5,5.5)
} else if (scenario == 5) {
  Beta = c(0.5,1.5,2.5)
  Gamma = c(0.5,1.5,NA)
  CensorC = c(1,3,2)
} else if (scenario == 6) {
  Beta = c(0.5,1,2)
  Gamma = c(0.5,1.5,NA)
  CensorC = c(1,3,2)
}

# fileConn<-file(paste(location, 'settings.txt', sep = ''))
# writeLines(c("sample size:", n_sample, '\n',
#              'Iteraton:', N_ite, '\n',
#              'scenario', scenario, '\n',
#              "Gamma (transition model):", Gamma1, '\n',
#              # "Gamma (transition model):", Gamma1, '\n',
#              # "Gamma (transition model):", Gamma2, '\n',
#              'Beta (event model):', Beta, '\n',
#              'Censor:', CensorC, '\n',
#              'Lambda1 (event)', Lambda1, '\n',
#              'Rho1 (event):', Rho1, '\n',
#              'Lambda2 (transition)', Lambda2, '\n',
#              'Rho2(transition):', Rho2, '\n',
#              'results folder location', location,'\n',
#              'No Censor?', nocensor,'\n',
#              'Time end (a very large number = not truncate):', time_end,'\n',
#              'correlation(smaller -> stronger, 1 = independent):', corr,'\n'), fileConn)
# close(fileConn)

p_value_01 = c()
p_value_12 = c()
p_value_02 = c()
for (i in 1:N_ite) {
  ############ Simulate Original Data
  set.seed(i)
  test = simulWeib_3stage(n_sample,
                          lambda1 = Lambda1, rho1 = Rho1,
                          lambda2 = Lambda2, rho2 = Rho2,
                          beta = Beta,
                          gamma = Gamma,
                          rateC = CensorC,
                          NoCensor = nocensor,
                          Trans_time)
  # test = simulWeib_3stage_td(n_sample,
  #                            lambda1 = Lambda1, rho1 = Rho1,
  #                            lambda2 = Lambda2, rho2 = Rho2,
  #                            beta = Beta,
  #                            gamma1 = Gamma1,
  #                            gamma2 = Gamma2,
  #                            rateC = CensorC,
  #                            NoCensor = nocensor)
  
  test2 = combine_stages(test[[1]], test[[2]], test[[3]])

  
  ## 
  TempTD = reformat_df_pairwise_end(test2, 0, 1, time_end)
  out = Mantel.Byar.test(0, 1, Group = "covariate_td", Event = TempTD$endpoint_td, 
                         StartTime = TempTD$start_td, StopTime = TempTD$stop_td, 
                         method = c('my method'), plot=0, landmark=0)
  table_out = out[["result"]]
  table = data.frame(n1 = table_out$na, da = table_out$da, n2 = table_out$nb, db = table_out$db)
  p_value_01[i] = out$p.value
  
  ## 
  
  TempTD = reformat_df_pairwise_end(test2, 1, 2, time_end)
  out = Mantel.Byar.test(1, 2, Group = "covariate_td", Event = TempTD$endpoint_td,
                         StartTime = TempTD$start_td, StopTime = TempTD$stop_td,
                         method = c('my method'), plot=0, landmark=0)
  p_value_12[i] = out$p.value

  ##
  ## modify one data plot to get dramatic drop
  # TempTD = reformat_df_pairwise_end(test2, 1, 2, time_end)
  # TempTD[which(TempTD$start_td == min(TempTD$start_td[which(TempTD$covariate_td == 2)]) & TempTD$covariate_td == 2),4] = 1
  # TempTD[which(TempTD$start_td == min(TempTD$start_td[which(TempTD$covariate_td == 2)]) & TempTD$covariate_td == 2),3] = TempTD[which(TempTD$start_td == min(TempTD$start_td[which(TempTD$covariate_td == 2)]) & TempTD$covariate_td == 2),2] + 0.0005
  # imc2 = TempTD %>% filter(covariate_td == 2)
  # model = survfit(Surv(time=start_td, time2=stop_td, event=endpoint_td) ~ covariate_td, TempTD)
  # a = ggsurvplot(model, data=TempTD, risk.table = T, tables.height = 0.3, conf.int = T, axes.offset =F)
  ##

  TempTD = reformat_df_pairwise_end(test2, 0, 2, time_end)
  TempTD[which(TempTD$id==47),4] = 1
  out = Mantel.Byar.test(0, 2, Group = "covariate_td", Event = TempTD$endpoint_td,
                         StartTime = TempTD$start_td, StopTime = TempTD$stop_td,
                         method = c('my method'), plot=0, landmark=0)

  p_value_02[i] = out$p.value
  
  ## 
  
  # TempTD = reformat_df_all(test2)
  # out = Mantel.Byar.overall(Group = "covariate_td", Event = TempTD$endpoint_td, 
  #                           StartTime = TempTD$start_td, StopTime = TempTD$stop_td, 
  #                           method = c("Tominaga"), plot=0, landmark=0)
  # p_value_all[i] = out$p.value
}

## type1 error / power
# mean(p_value_01<=0.05, na.rm = T)
# mean(p_value_02<=0.05, na.rm = T)
# mean(p_value_12<=0.05, na.rm = T)

mean(p_value_01<=0.05)
mean(p_value_02<=0.05)
mean(p_value_12<=0.05)

# pvalue = data.frame(group01 = p_value_01,
#                     group12 = p_value_12,
#                     group02 = p_value_02)

# write.csv(pvalue, paste(location, 'pvalue.csv', sep = ''))



# data <- data.frame(
#   n1 = c(200, 199, 198, 197, 196, 195, 194, 193, 192, 191, 190, 189, 188, 187, 186, 185, 184, 183, 182, 181, 180, 179, 178, 0),
#   da = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0),
#   n2 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 178),
#   db = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
# )
# 
# # Summarize the data to create a two by two table
# table_data <- data.frame(
#   da = c(sum(data$n1[data$da == 1]), sum(data$n1[data$da == 0])),
#   db = c(sum(data$n2[data$db == 1]), sum(data$n2[data$db == 0]))
# )
# 
# rownames(table_data) <- c("da", "not da")
# colnames(table_data) <- c("db", "not db")
# 
# table_data
