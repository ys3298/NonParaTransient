library(tidyverse, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
library("survminer", lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
library("cmprsk", lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
library(rsample, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')

## In CR, there are some combinations (0,1), (1,1), (2,1). 
## The first number is the level of grouping variable, and the second are the outcome of interests

# ### PART1: raw ROC and AUC
# cal_auc = function(CR_out) {
#   interests = names(CR_out)[which(substr(names(CR_out), 3, 3) == "1")]
#   ind = which(substr(names(CR_out), 3, 3) == "1")
#   t_largest = 0
#   for (i in 1:length(ind)) {
#     if (max(CR_out[[i]]$time) > t_largest) {
#       t_largest = max(CR_out[[i]]$time)
#     }
#   }
#   for (m in 1:length(ind)) {
#     CR_out[[m]]$time[length(CR_out[[m]]$time)] = t_largest
#   }
# 
#   AUC = c()
#   for (i in 1:length(ind)) {
#     index = ind[i]
#     df = data.frame(time = CR_out[[index]]$time, est = CR_out[[index]]$est)
#     AUC[i] = sum(diff(df$time, lag = 1) * head(df$est, -1))
#   }
#   names(AUC) = paste('group', substr(names(CR_out)[ind], 1, 1), sep='')
#   return(AUC)
# }
 
ggcompetingrisks = function(fit, gnames = NULL, gsep=" ",
                            multiple_panels = TRUE,
                            ggtheme = theme_survminer(),
                            coef = 1.96, conf.int = FALSE, ...) {
  stopifnot(any(class(fit) %in% c("cuminc", "survfitms")))

  if (any(class(fit) == "cuminc")) {
    pl = ggcompetingrisks.cuminc(fit = fit, gnames=gnames,
                                 gsep=gsep, multiple_panels=multiple_panels,
                                 coef = coef, conf.int = conf.int)
  }
  if (any(class(fit) == "survfitms")) {
    pl = ggcompetingrisks.survfitms(fit = fit)
  }

  pl = pl + ggtheme +
    ylab("Probability of an event") + xlab("Time") +
    ggtitle("Cumulative incidence functions")
  ggpubr::ggpar(pl, ...)
}

ggcompetingrisks.cuminc = function(fit, gnames = NULL, gsep=" ",
                                   multiple_panels = TRUE, coef = 1.96, conf.int = FALSE) {
  if (!is.null(fit$Tests))
    fit <- fit[names(fit) != "Tests"]
  fit2 <- lapply(fit, `[`, 1:3)
  if (is.null(gnames)) gnames <- names(fit2)
  fit2_list <- lapply(seq_along(gnames), function(ind) {
    df <- as.data.frame(fit2[[ind]])
    df$name <- gnames[ind]
    df
  })
  time <- est <- event <- group <- NULL
  df <- do.call(rbind, fit2_list)
  df$event <- sapply(strsplit(df$name, split=gsep), `[`, 2)
  df$group <- sapply(strsplit(df$name, split=gsep), `[`, 1)
  df$std <- std <- sqrt(df$var)
  # pl <- ggplot(df, aes(time, est, color=event))
  df2 = df %>% filter(event == 1) %>%
    rename('IMC' = 'group')

  if (multiple_panels) {
    pl <- ggplot(df, aes(time, est)) + facet_wrap(~IMC)
  } else {
    df2 = df %>% filter(event == 1) %>%
      rename('IMC' = 'group')
    pl <- ggplot(df2, aes(time, est, color=event, linetype = IMC))
  }
  if (conf.int) {
    pl <- pl + geom_ribbon(aes(ymin = est - coef*std, ymax=est + coef*std), alpha = 0.2, linetype=0)
  }
  pl +
    geom_line() +scale_color_manual(name="", values=c("blue","red"), labels=c("Event", "Transit"))
}

# CR <- cuminc(ftime = data$time_diff,
#              fstatus = data$status,
#              cencode = 0,
#              group = data$IMC)
# ggcompetingrisks(fit = CR, multiple_panels = F, 
#                  xlab = "Days", 
#                  ylab = "Cumulative incidence of Event", 
#                  title = "Competing Risks Analysis")


### PART2: calibrated ROC and AUC

auc = function(df) {
  time = df$time
  est = df$est
  area = 0
  # if (sum(!is.na(est)) != 0) {
    for (i in 1:(nrow(df)-1)) {
      est_now = est[i]
      est_next = est[i+1]
      if (est_now == est_next) {
        area = area + est_now*(time[i+1] - time[i])
      } else if (est_now != est_next) {
        area = area + (est_now + est_next)*(time[i+1] - time[i])/2
      }
    }
  # }
  return(area)
}

# version 2 by 08/01/2021
# cal_auc_calibrate = function(CR_out, raw_df, time_end = 100000) {
#   # time_end could be manually chosen, calculate the AUC up to that time point,
#   # if max time < time end, it will stretch the last observation to the time end
#   interests = names(CR_out)[which(substr(names(CR_out), 3, 3) == "1")]
#   ind = which(substr(names(CR_out), 3, 3) == "1")
#   t_largest = 0
#   for (i in 1:length(ind)) {
#     if (max(CR_out[[i]]$time) > t_largest) {
#       t_largest = max(CR_out[[i]]$time)
#     }
#   }
#   for (m in 1:length(ind)) {
#     CR_out[[m]]$time[length(CR_out[[m]]$time)] = t_largest
#   }
#   
#   # print(t_largest)
#   AUC = c()
#   all_df = data.frame(time = double(),
#                       est = double(),
#                       group = double())
#   for (i in 1:length(ind)) {
#     index = ind[i]
#     if (index != length(ind)) { # calibrate the survival funtion for stages before the last
#       df = data.frame(time = CR_out[[index]]$time, est = CR_out[[index]]$est, group = index)
#       # get trans_prob
#       data_subgroup = raw_df %>% filter(IMC == index-1)
#       prob_trans = sum(data_subgroup$status == 2 & data_subgroup$time_diff <= max(data_subgroup$time_diff)) / sum(data_subgroup$status != 0 & data_subgroup$time_diff <= max(data_subgroup$time_diff))
# 
#       data_km = raw_df %>% mutate(status = ifelse(status == 2, 1, 0)) # transition = event (1); censor + event of interest = censor (0)
#       data_km_subgroup = data_km %>% filter(IMC == index-1) # for each group separately
#       km_fit0 <- survfit(Surv(time_diff, status) ~ 1, data=data_km_subgroup)
#       table1 = summary(km_fit0, times = CR_out[[index]]$time)
#       factor = data.frame(time = table1$time,
#                           survival = table1$surv,
#                           one_minus_sur = 1-table1$surv,
#                           one_minus_sur_fac = (1-table1$surv)*prob_trans,
#                           sur_fac = 1 - (1-table1$surv)*prob_trans,
#                           factor = 1/(1 - (1-table1$surv)*prob_trans))
#       if (max(factor$time) == t_largest) {
#         factor_vec = factor$factor
#       } else {
#         factor_vec = c(factor$factor, factor$factor[length(factor$factor)])
#       }
#       
#       df$est = df$est*factor_vec
#       
#       ##
#       if(df$time[length(df$time)-1] < time_end & t_largest > time_end) {
#         df$time[[length(df$time)]] = time_end
#       } else if(df$time[length(df$time)-1] > time_end) {
#         df$est[length(df$est)] = df$est[max(which(df$time <= time_end))]
#         df$time[[length(df$time)]] = time_end
#       } else if(df$time[length(df$time)] < time_end) {
#         df$time[[length(df$time)]] = time_end
#       }
#       df = df %>% filter(time <= time_end)
#       ##
#       
#       all_df = bind_rows(all_df, df)
#     } else {
#       df = data.frame(time = CR_out[[index]]$time, est = CR_out[[index]]$est, group = index)
#       
#       ##
#       if(df$time[length(df$time)-1] < time_end & t_largest > time_end) {
#         df$time[[length(df$time)]] = time_end
#       } else if(df$time[length(df$time)-1] > time_end) {
#         df$est[length(df$est)] = df$est[max(which(df$time <= time_end))]
#         df$time[[length(df$time)]] = time_end
#       } else if(df$time[length(df$time)] < time_end) {
#         df$time[[length(df$time)]] = time_end
#       }
#       df = df %>% filter(time <= time_end)
#       ##
#       
#       all_df = bind_rows(all_df, df)
#     }
#     
#     # AUC[i] = sum(diff(df$time, lag = 1) * head(df$est, -1))
#     AUC[i] = auc(df)
#     
#   }
#   names(AUC) = paste('group', substr(names(CR_out)[ind], 1, 1), sep='')
#   return(list(AUC, all_df))
# }

# version3: 08/01/2023
# cal_auc_calibrate = function(CR_out, raw_df, time_end = 100000) {
#   # time_end could be manually chosen, calculate the AUC up to that time point,
#   # if max time < time end, it will stretch the last observation to the time end
#   interests = names(CR_out)[which(substr(names(CR_out), 3, 3) == "1")]
#   ind = which(substr(names(CR_out), 3, 3) == "1")
#   t_largest = 0
#   for (i in 1:length(ind)) {
#     if (max(CR_out[[i]]$time) > t_largest) {
#       t_largest = max(CR_out[[i]]$time)
#     }
#   }
#   for (m in 1:length(ind)) {
#     CR_out[[m]]$time[length(CR_out[[m]]$time)] = t_largest
#   }
# 
#   # print(t_largest)
#   AUC = c()
#   all_df = data.frame(time = double(),
#                       est = double(),
#                       group = double())
#   for (i in 1:length(ind)) {
#     index = ind[i]
#     if (index != length(ind)) { # calibrate the survival function for stages before the last
#       df = data.frame(time = CR_out[[index]]$time, est = CR_out[[index]]$est, group = index)
#       # get trans_prob
#       data_subgroup = raw_df %>% filter(IMC == index-1)
#       time_prob_trans =  max(data_subgroup$time_diff)
#       prob_trans = sum(data_subgroup$status == 2 & data_subgroup$time_diff <=time_prob_trans) /
#         sum(data_subgroup$status != 0 & data_subgroup$time_diff <=time_prob_trans)
# 
#       # data_km = raw_df %>% mutate(status = ifelse(status == 2, 1, 0)) # transition = event (1); censor + event of interest = censor (0)
#       # data_km_subgroup = data_km %>% filter(IMC == index-1) # for each group separately
#       # km_fit0 <- survfit(Surv(time_diff, status) ~ 1, data=data_km_subgroup)
#       # table1 = summary(km_fit0, times = CR_out[[index]]$time)
# 
#       km_fit0 = cuminc(ftime = data_subgroup$time_diff,
#                        fstatus = data_subgroup$status,
#                        cencode = 0)
#       table1 = list(time = unique(CR_out[[index]]$time),
#                     est = timepoints(km_fit0, times = CR_out[[index]]$time)$est[2,] ) # transition=2 CIF estimation
#       # table1 = km_fit0$`1 2` # transition
#       # plot(x=table1$time, y=table1$est)
#       # title(paste('func1',i))
#       table1$surv = 1-table1$est
#       factor = data.frame(time = table1$time,
#                           survival = table1$surv,
#                           one_minus_sur = 1-table1$surv,
#                           one_minus_sur_fac = (1-table1$surv)*prob_trans,
#                           sur_fac = 1 - (1-table1$surv)*prob_trans,
#                           factor = 1/(1 - (1-table1$surv)*prob_trans))
#       temp = left_join(df, factor, by = 'time')
#       temp$factor = ifelse(is.na(temp$factor), temp$factor[length(na.omit(temp$factor))], temp$factor)
#       df$est = df$est*temp$factor
# 
#       ##
#       if(df$time[length(df$time)-1] < time_end & t_largest > time_end) {
#         df$time[[length(df$time)]] = time_end
#       } else if(df$time[length(df$time)-1] > time_end) {
#         df$est[length(df$est)] = df$est[max(which(df$time <= time_end))]
#         df$time[[length(df$time)]] = time_end
#       } else if(df$time[length(df$time)] < time_end) {
#         df$time[[length(df$time)]] = time_end
#       }
#       df = df %>% filter(time <= time_end)
#       ##
# 
#       all_df = bind_rows(all_df, df)
#     } else {
#       df = data.frame(time = CR_out[[index]]$time, est = CR_out[[index]]$est, group = index)
# 
#       ##
#       if(df$time[length(df$time)-1] < time_end & t_largest > time_end) {
#         df$time[[length(df$time)]] = time_end
#       } else if(df$time[length(df$time)-1] > time_end) {
#         df$est[length(df$est)] = df$est[max(which(df$time <= time_end))]
#         df$time[[length(df$time)]] = time_end
#       } else if(df$time[length(df$time)] < time_end) {
#         df$time[[length(df$time)]] = time_end
#       }
#       df = df %>% filter(time <= time_end)
#       ##
# 
#       all_df = bind_rows(all_df, df)
#     }
# 
#     # print(df)
#     AUC[i] = auc(df)
# 
#   }
#   names(AUC) = paste('group', substr(names(CR_out)[ind], 1, 1), sep='')
#   return(list(AUC, all_df))
# }


# version4: 08/18/2023
cal_auc_calibrate = function(CR_out, raw_df, time_end = 100000) {
  # time_end could be manually chosen, calculate the AUC up to that time point,
  # if max time < time end, it will stretch the last observation to the time end
  interests = names(CR_out)[which(substr(names(CR_out), 3, 3) == "1")]
  ind = which(substr(names(CR_out), 3, 3) == "1")
  t_largest = 0
  for (i in 1:length(ind)) {
    if (max(CR_out[[i]]$time) > t_largest) {
      t_largest = max(CR_out[[i]]$time)
    }
  }
  for (m in 1:length(ind)) {
    CR_out[[m]]$time[length(CR_out[[m]]$time)] = t_largest
  }
  
  AUC = c()
  all_df = data.frame(time = double(),
                      est = double(),
                      group = double())

  AUCt = data.frame(matrix(NA, nrow = 3, ncol = 5-1))
  CIFt = data.frame(matrix(NA, nrow = 3, ncol = 5-1))
  for (i in 1:length(ind)) {
    index = ind[i]
    if (index != length(ind)) { # calibrate the survival function for stages before the last
      df = data.frame(time = CR_out[[index]]$time, est = CR_out[[index]]$est, group = index)
      # get trans_prob
      data_subgroup = raw_df %>% filter(IMC == index-1)
      time_prob_trans =  min(max(data_subgroup$time_diff), time_end) ## change from last time point to pre-defined time-end
      # time_prob_trans =  (max(data_subgroup$time_diff))
      prob_trans = sum(data_subgroup$status == 2 & data_subgroup$time_diff <= time_prob_trans) /
        sum(data_subgroup$status != 0 & data_subgroup$time_diff <= time_prob_trans) # alpha_z
      # print(prob_trans)
      km_fit0 = cuminc(ftime = data_subgroup$time_diff,
                       fstatus = data_subgroup$status,
                       cencode = 0)

      table1 = list(time = unique(CR_out[[index]]$time),
                    est = timepoints(km_fit0, times = CR_out[[index]]$time)$est["1 2",] ) # transition=2 CIF estimation
      table1$est = ifelse(is.na(table1$est), max(km_fit0[["1 2"]][["est"]]),table1$est)
      # table1 = km_fit0$`1 2` # transition
      # plot(x=table1$time, y=table1$est)
      # title(paste('func1',i))

      table1$surv = 1-table1$est
      factor = data.frame(time = table1$time,
                          survival = table1$surv,
                          one_minus_sur = 1-table1$surv,
                          one_minus_sur_fac = (1-table1$surv)*prob_trans,
                          sur_fac = 1 - (1-table1$surv)*prob_trans,
                          factor = 1/(1 - (1-table1$surv)*prob_trans))
      temp = left_join(df, factor, by = 'time')
      temp$factor = ifelse(is.na(temp$factor), temp$factor[length(na.omit(temp$factor))], temp$factor)
      # when no event (df$est = 0), F_tran = 1, factor = inf, let factor = 1, since est = 0, the value of factor will not influence estimation
      if (all(any(is.infinite(temp$factor)) & (df$est == 0))) {
        temp$factor = ifelse(is.infinite(temp$factor), 1.0, temp$factor)
      }
      if (all(is.na(temp$factor))) {
        temp$factor = ifelse((is.na(temp$factor)), 1.0, temp$factor)
      }
      df$est = df$est*temp$factor
      df_till_end = df
      ##
      if(df$time[length(df$time)-1] < time_end & t_largest > time_end) {
        df$time[[length(df$time)]] = time_end
      } else if(df$time[length(df$time)-1] > time_end) {
        df$est[length(df$est)] = df$est[max(which(df$time <= time_end))]
        df$time[[length(df$time)]] = time_end
      } else if(df$time[length(df$time)] < time_end) {
        df$time[[length(df$time)]] = time_end
      }
      df = df %>% filter(time <= time_end)
      ##
      all_df = bind_rows(all_df, df)
    } else {
      df = data.frame(time = CR_out[[index]]$time, est = CR_out[[index]]$est, group = index)
      df_till_end = df
      ##
      if(df$time[length(df$time)-1] < time_end & t_largest > time_end) {
        df$time[[length(df$time)]] = time_end
      } else if(df$time[length(df$time)-1] > time_end) {
        df$est[length(df$est)] = df$est[max(which(df$time <= time_end))]
        df$time[[length(df$time)]] = time_end
      } else if(df$time[length(df$time)] < time_end) {
        df$time[[length(df$time)]] = time_end
      }
      df = df %>% filter(time <= time_end)
      ##
      
      all_df = bind_rows(all_df, df)
    }
    
    AUC[i] = auc(df)
    
    # ## AUC(t)
    # time_length = 9
    # time_end_plot = 0.4
    # time_cut = seq(0,time_end_plot,length.out = time_length)
    # for (j in 2:time_length) {
    #   df_temp = df_till_end %>% filter(time <= time_cut[j])
    #   df_temp = rbind(df_temp, df_temp[nrow(df_temp),])
    #   df_temp[nrow(df_temp),1] = time_cut[j]
    #   if (nrow(df_temp) <= 1) {
    #     AUCt[i,j-1] = 0
    #   } else {
    #     AUCt[i,j-1] = auc(df_temp)
    #   }
    # }
    # 
    # ## CIF(t)
    # for (j in 2:time_length) {
    #   df_temp = df_till_end %>% filter(time <= time_cut[j])
    #   CIFt[i,j-1] = df_temp$est[length(df_temp$est)]
    # }
  }
  
  
  names(AUC) = paste('group', substr(names(CR_out)[ind], 1, 1), sep='')
  # return(list(AUC, all_df, AUCt, CIFt))
  return(list(AUC, all_df))
}

cal_auc_raw = function(CR_out, raw_df, time_end = 100000) {
  # time_end could be manually chosen, calculate the AUC up to that time point,
  # if max time < time end, it will stretch the last observation to the time end
  interests = names(CR_out)[which(substr(names(CR_out), 3, 3) == "1")]
  ind = which(substr(names(CR_out), 3, 3) == "1")
  t_largest = 0
  for (i in 1:length(ind)) {
    if (max(CR_out[[i]]$time) > t_largest) {
      t_largest = max(CR_out[[i]]$time)
    }
  }
  for (m in 1:length(ind)) {
    CR_out[[m]]$time[length(CR_out[[m]]$time)] = t_largest
  }
  
  AUC = c()
  all_df = data.frame(time = double(),
                      est = double(),
                      group = double())
  
  AUCt = data.frame(matrix(NA, nrow = 3, ncol = 5-1))
  CIFt = data.frame(matrix(NA, nrow = 3, ncol = 5-1))
  for (i in 1:length(ind)) {
    index = ind[i]
    if (index != length(ind)) { # calibrate the survival function for stages before the last
      df = data.frame(time = CR_out[[index]]$time, est = CR_out[[index]]$est, group = index)
      df_till_end = df
      ##
      if(df$time[length(df$time)-1] < time_end & t_largest > time_end) {
        df$time[[length(df$time)]] = time_end
      } else if(df$time[length(df$time)-1] > time_end) {
        df$est[length(df$est)] = df$est[max(which(df$time <= time_end))]
        df$time[[length(df$time)]] = time_end
      } else if(df$time[length(df$time)] < time_end) {
        df$time[[length(df$time)]] = time_end
      }
      df = df %>% filter(time <= time_end)
      ##
      all_df = bind_rows(all_df, df)
    } else {
      df = data.frame(time = CR_out[[index]]$time, est = CR_out[[index]]$est, group = index)
      df_till_end = df
      ##
      if(df$time[length(df$time)-1] < time_end & t_largest > time_end) {
        df$time[[length(df$time)]] = time_end
      } else if(df$time[length(df$time)-1] > time_end) {
        df$est[length(df$est)] = df$est[max(which(df$time <= time_end))]
        df$time[[length(df$time)]] = time_end
      } else if(df$time[length(df$time)] < time_end) {
        df$time[[length(df$time)]] = time_end
      }
      df = df %>% filter(time <= time_end)
      ##
      
      all_df = bind_rows(all_df, df)
    }
    
    AUC[i] = auc(df)
    
    ## AUC(t)
    time_length = 9
    time_end_plot = 0.4
    time_cut = seq(0,time_end_plot,length.out = time_length)
    for (j in 2:time_length) {
      df_temp = df_till_end %>% filter(time <= time_cut[j])
      df_temp = rbind(df_temp, df_temp[nrow(df_temp),])
      df_temp[nrow(df_temp),1] = time_cut[j]
      if (nrow(df_temp) <= 1) {
        AUCt[i,j-1] = 0
      } else {
        AUCt[i,j-1] = auc(df_temp)
      }
    }
    
    ## CIF(t)
    for (j in 2:time_length) {
      df_temp = df_till_end %>% filter(time <= time_cut[j])
      CIFt[i,j-1] = df_temp$est[length(df_temp$est)]
    }
  }
  
  
  names(AUC) = paste('group', substr(names(CR_out)[ind], 1, 1), sep='')
  return(list(AUC, all_df, AUCt, CIFt))
}

plot_calibrate = function(output) {
  df = output[[2]]
  df = df %>% mutate(group = group - 1) %>% mutate(group = factor(group))
  ggplot(df) + geom_line(aes(x = time, y = est, color = group))
}



# source('/storage/work/yqs5519/SurvivalPlot/algorithm/generate_correlated_data.R')
# source('/storage/work/yqs5519/SurvivalPlot/algorithm/generate_data.R')

# set.seed(20)
# test = fun_gendata_modify(n = 1000,
#                           rho = 1,
#                           beta = c(1,1,1),
#                           gamma = c(1.5,1.5,NA),
#                           rateC = c(2.7,5,3.0))
# test2 = combine_stages(test[[1]], test[[2]], test[[3]])
# data = reformat_df(test2)
# CR <- cuminc(ftime = data$time_diff,
#              fstatus = data$status,
#              cencode = 0,
#              group = data$IMC)
# plot(CR)
# a = cal_auc_raw(CR, data, time_end = 0.2)
# plot_calibrate(a)
# a[[1]]
# 
# a = cal_auc_calibrate(CR, data, time_end = 0.2)
# plot_calibrate(a)
# a[[1]]

# C_Survival_weibull1 = function(t) {
#   1-exp(-1*(t^1)*exp(c(1,1,1)))
# }
# 
# time_end = 0.2
# b = cal_auc_calibrate(CR, data, time_end = time_end)

# ### cif
# df = b[[2]]
# true = data.frame(time = seq(0, time_end, 0.01), est = C_Survival_weibull1(seq(0, time_end, 0.01)), group = 'true')
# df = rbind(true,df) %>% mutate(group = factor(group))
# ggplot(df) + geom_line(aes(x = time, y = est, color = group))
# c(b[[1]], true = integrate(C_Survival_weibull1, 0, time_end)$value)
# 
# 
# ### auc
# b[[3]]
# True_AUCt1 = c()
# time_point = seq(0,0.8,length.out = 17)
# for (i in 1:length(time_point)) {
#   True_AUCt1[i] =  integrate(C_Survival_weibull1, 0, time_point[i])$value
# }
# # b[[3]]
# # True_AUCt1[-1]
# bias = rbind(group1 = b[[3]][1,] - True_AUCt1[-1],
#              group2 = b[[3]][2,] - True_AUCt1[-1],
#              group3 = b[[3]][3,] - True_AUCt1[-1])
# bias
# b[[3]]
