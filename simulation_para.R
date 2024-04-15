assign(".lib.loc", "/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2", envir = environment(.libPaths))

rm(list = ls(all.names = TRUE))
start_time = Sys.time()
library(rlang, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
library(tidyverse, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
source('/storage/work/yqs5519/SurvivalPlot/algorithm/generate_data.R')
source('/storage/work/yqs5519/SurvivalPlot/algorithm/calibrateAUC.R')
source('/storage/work/yqs5519/SurvivalPlot/algorithm/generate_correlated_data.R')
library(foreach, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
library(doParallel, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')

file_num = readr::parse_number(list.files('/storage/work/yqs5519/SurvivalPlot/', 'formal_sim_results'))
location_index = ifelse((is_empty(file_num)), 1,  max(file_num, na.rm = T) + 1) ## set location to save results
location = paste('/storage/work/yqs5519/SurvivalPlot/formal_sim_results', location_index, '/', sep = '')
dir.create(location) 

## settings
N_ite = 1000
B = 1000
corr = 1

n_sample = 500
scenario = 3
Time_End = 0.25
nocensor = F
# T = no censor
# F = have censor
noCalibration = F
# T = cal_auc_raw
# F = cal_auc_calibrate



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
  Gamma = c(1.3,1.3,NA)
  CensorC = c(1,5,5.5)
} else if (scenario == 4) {
  Beta = c(0.5,1,1.5)
  Gamma = c(0.5,1.5,NA)
  CensorC = c(1,5,5.5)
} 

Lambda1 = 1
Rho1 = 1
Lambda2 = 1
Rho2 = 1

fileConn<-file(paste(location, 'settings.txt', sep = ''))
writeLines(c("sample size:", n_sample, '\n',
             'last time:', Time_End, '\n',
             'boostrap B:', B, '\n',
             'Iteraton:', N_ite, '\n',
             "Gamma (transition model):", Gamma, '\n',
             'Beta (event model):', Beta, '\n',
             'Censor:', CensorC, '\n',
             'Lambda1 (event)', Lambda1, '\n',
             'Rho1 (event):', Rho1, '\n',
             'Lambda2 (transition)', Lambda2, '\n',
             'Rho2(transition):', Rho2, '\n',
             'results folder location', location,'\n',
             'No Censor?', nocensor,'\n',
             'correlation(smaller -> stronger, 1 = independent):', corr,'\n',
             'No calibration?', noCalibration,'\n'), fileConn)
close(fileConn)

## simulate
multiResultClass <- function(result1=NULL,result2=NULL,result3=NULL,result4=NULL,result5=NULL)
{
  me <- list(
    orig_stat = result1,
    stat_list = result2,
    AUCt_list = result3,
    CIFt_list = result4,
    quant80 = result5
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}

stat_list = list()
AUCt_list = list()
CIFt_list = list()
orig_stat = data.frame(matrix(ncol = 3, nrow = 0))
quant80 = data.frame(matrix(ncol = 0, nrow = 3))
num_cores <- detectCores()-1
cl <- makeCluster(num_cores)
registerDoParallel(cl)


out = foreach(i = 1:N_ite) %dopar%  {

  library(rsample, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
  library(simsurv, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
  library(dplyr, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
  library(cmprsk, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
  library(tidyverse, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
  library(survsim, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
  library(survival, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
  library(survminer, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
  library(MASS, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
  
  ############ Simulate Original Data
  set.seed(i)
  result <- multiResultClass()

  test = simulWeib_3stage(n_sample,
                          lambda1 = Lambda1, rho1 = Rho1,
                          lambda2 = Lambda2, rho2 = Rho2,
                          beta = Beta,
                          gamma = Gamma,
                          rateC = CensorC,
                          NoCensor = nocensor,
                          trans_time = 'no')

  
  test2 = combine_stages(test[[1]], test[[2]], test[[3]])
  data = reformat_df(test2)

  
  CR0 <- cuminc(ftime = data$time_diff,
                fstatus = data$status,
                cencode = 0,
                group = data$IMC)

  result$orig_stat = rbind(result$orig_stat, cal_auc_calibrate(CR0, data, time_end = Time_End)[[1]])
  
  
  #### Bootstrap 
  D <- nest(data,-id)
  bs <- bootstraps(D, times = B)
  stat = data.frame(matrix(ncol = 3, nrow = 0))
  for (b in 1:B) {
    b_data = as.tibble(bs$splits[[b]]) %>% arrange(id) %>% unnest(cols = c(data))

    CR <- cuminc(ftime = b_data$time_diff,
                 fstatus = b_data$status,
                 cencode = 0,
                 group = b_data$IMC)

    # some function to calculate the statistic
    stat[b,] = cal_auc_calibrate(CR, b_data, time_end = Time_End)[[1]]
    colnames(stat) = names(orig_stat)
  }
  result$stat_list = c(result$stat_list, list(stat))

  return(result)
}



# Stop the parallel processing
stopCluster(cl)

orig_stat = do.call(rbind, lapply(out, function(x) x$orig_stat))
write.csv(orig_stat, paste(location, 'original_statistic.csv', sep = ''), col.names = F)

stat_list = do.call(c, lapply(out, function(x) x$stat_list))
saveRDS(stat_list, paste(location, 'stat_list.rds', sep = ''))

time = Sys.time() - start_time
