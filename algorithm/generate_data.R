# simulate three stage data
library(simsurv, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
library(tidyverse, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
library(dplyr, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
library("cmprsk", lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
library(survsim, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
library(survival, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
library("survminer", lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
library("cmprsk", lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
library("MASS", lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
round_func = function(value_vec, mean_min = 'min') {
  if (mean_min == 'min') {
    original_value <- min(value_vec)  # The original value
  }
  if (mean_min == 'median') {
    original_value <- median(value_vec)  # The original value
  }
  rounded_value <- round(original_value / 0.05) * 0.05  # Round to the nearest multiple of 0.05
  
  # Ensure the rounded value is less than the original value
  # if (rounded_value >= original_value) {
  #   rounded_value <- rounded_value - 0.05
  # }
  return(rounded_value)
}

###################################################### Data generation
# baseline hazard: Weibull
# N = sample size
# covariate X (for example IMC) with level 0,1,2 with probability p
# lambda1, rho1 = scale and shape parameters in h0() in the weibull event process
# lambda2, rho2 = scale and shape parameters in h0() in the transition event process
# beta = fixed effect parameter for event
# gamma = fixed effect parameter for transition
# rateC = rate parameter of the exponential censoring distribution of C
# https://stats.stackexchange.com/questions/135124/how-to-create-a-toy-survival-time-to-event-data-with-right-censoring


simulWeib_3stage = function(N=500, 
                            lambda1 = 1, rho1 = 1, 
                            lambda2 = 1, rho2 = 1, 
                            beta = c(1,1,1), 
                            gamma = c(0.5,0.5,0.5), 
                            rateC = c(1,2,3),
                            NoCensor = FALSE,
                            trans_time)
{
  df = list()
  # x = IMC status, start from all 0
  x = rep(0,N)
  for (i in 1:3) {
    # if IMC = 0, take the first parameter; 
    # if IMC = 1, take the second, etc..
    beta_x = c()
    gamma_x = c()
    rateC_x = c()
    beta_x = beta[1+x] # get the (1+x)th of beta for N times
    gamma_x = gamma[1+x]
    rateC_x = rateC[1+x]
    
    # Weibull latent event times
    v = runif(n=N)
    T1lat = (- log(v) / (lambda1 * exp(beta_x)))^(1 / rho1)
    
    # Weibull latent transition times
    t = runif(n=N)
    T2lat = (- log(t) / (lambda2 * exp(gamma_x)))^(1 / rho2)
    
    # censoring times
    if (NoCensor) {C = rep(100000, N)} else{C = rexp(n=N, rate=rateC_x)}
    
    
    ## Get the status and time for the status;
    ## When IMC = 2 (i = 3, last stage), it is not possible to transit further. 
    ## So just compare Event time and censor time.
    if (trans_time == 'no') {
      if (i == 3) {
        # follow-up times and event indicators
        time = pmin(T1lat, C)
        status = case_when(time == T1lat ~ 1,
                           time == C ~ 0)

      } else {
        # follow-up times and event indicators
        time = pmin(T1lat, T2lat, C)
        status = case_when(time == T1lat ~ 1,
                           time == T2lat ~ 2,
                           time == C ~ 0)
        ## event = 1; transition = 2;
      }
    } 
    if (trans_time != 'no') {
      if (i == 3) {
        # follow-up times and event indicators
        time = pmin(T1lat, C)
        status = case_when(time == T1lat ~ 1,
                           time == C ~ 0)
        
      } else {
        # follow-up times and event indicators
        time = c()
        status = c()
        for (m in 1:length(T2lat)) {
          time[m] <- min(T1lat[m], C[m], trans_time)
          status[m] = case_when(time[m] == T1lat[m] ~ 1,
                                time[m] == trans_time ~ 2,
                                time[m] == C[m] ~ 0)
        }
        ## event = 1; transition = 2;
      }
    }
    
    # update x base on transition or not
    x = ifelse(status == 2, x+1, x)
    # save the data set for each stage
    df[[i]] = data.frame(id=1:N,time=time,status=status)
  }
  
  return(df)
}


combine_stages = function(df1,df2,df3) {
  n = nrow(df1) 
  df = data.frame(matrix(ncol = 6, nrow = n))
  colnames(df) = c('id', 'IMC', 'drop1', 'drop2', 'time', 'status')
  for (i in 1:n) {
    if (df1$status[i] == 1) {
      df[i,] = c(i, 0, NA, NA, df1$time[i], df1$status[i])
    } else if (df1$status[i] == 0) {
      df[i,] = c(i, 0, NA, NA, df1$time[i], df1$status[i])
    } else if (df1$status[i] == 2 & df2$status[i] != 2) {
      df[i,] = c(i, 1, df1$time[i], NA, df1$time[i]+df2$time[i], df2$status[i])
    } else if (df1$status[i] == 2 & df2$status[i] == 2) {
      df[i,] = c(i, 2, df1$time[i], df1$time[i]+df2$time[i], df1$time[i]+df2$time[i]+df3$time[i], df3$status[i])
    } 
  }
  return(df)
}


## reformat data
reformat_df = function(d1)
{
  n = dim(d1)[1]
  dd = NULL
  for(i in 1:n)
  {
    one_row = d1[i,]
    #################################  
    if(d1$IMC[i]==2) 
    {
      #subject, time1, time2, DFS_STATUS, IMC
      #DFS_STATUS=1 event; DFS_STATUS=2 transit to next; DFS_STATUS=0 censor;
      
      temp1 = c(one_row$id, 0, one_row$drop1, 2, 0)
      dd = rbind(dd, temp1)
      
      temp2 = c(one_row$id, one_row$drop1, one_row$drop2,  2, 1)
      dd = rbind(dd, temp2)
      
      temp3 = c(one_row$id, one_row$drop2, one_row$time, one_row$status, 2)
      dd = rbind(dd, temp3)
      
    }
    ####################################
    if(d1$IMC[i]==1) 
    {
      temp1 = c(one_row$id, 0, one_row$drop1, 2, 0)
      dd = rbind(dd, temp1)
      
      temp2 = c(one_row$id, one_row$drop1, one_row$time,  one_row$status, 1)
      dd = rbind(dd, temp2)
    }
    ####################################
    
    if(d1$IMC[i]==0) 
    {
      temp2 = c(one_row$id, 0, one_row$time, one_row$status, 0)
      dd = rbind(dd, temp2)
    }
    
  }
  
  colnames(dd) = c('id', 'time1', 'time2', 'status', 'IMC')
  dd = as.data.frame(dd) %>% mutate(time_diff = time2 - time1)
  return(dd)
}

# test2 = combine_stages(test[[1]], test[[2]], test[[3]])
# data = reformat_df(test2)
