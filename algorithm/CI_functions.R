library(tidyverse, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
library("survminer", lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
library("cmprsk", lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')
library(rsample, lib.loc = '/storage/home/yqs5519/R/x86_64-pc-linux-gnu-library/4.2')

################################################################################### the following are CI and p-value
## percentile
percentile_ci = function(statisitcs, alpha = 0.05, side = "two-sided") {
  if (side == "two-sided") {
    CI_l = quantile(statisitcs, alpha/2)
    CI_u = quantile(statisitcs, 1-alpha/2)

    CI = c(CI_l, CI_u)
  }

  return(CI)
}

reverse_percentile_ci = function(original_stat, b_stat, alpha = 0.05, side = "two-sided") {
  if (side == "two-sided") {
    CI_l = 2*original_stat - quantile(b_stat, 1-alpha/2)
    CI_u = 2*original_stat - quantile(b_stat, alpha/2)

    CI = c(CI_l[[1]], CI_u[[1]]) %>% as.numeric()
    names(CI) = c(paste((alpha/2)*100,"%"), paste((1-alpha/2)*100,"%"))
  }
  return(CI)
}

# student_t_ci = function(original_stat, b_stat, alpha = 0.05, side = "two-sided") {
#   n = length(b_stat)
#   # se = sqrt((sum((b_stat - mean(b_stat))^2))/(n-1)) # todo check
#   se = sd(b_stat) / sqrt(length(b_stat))
#   # http://www.comparingpartitions.info/?link=Tut9
#   if (side == "two-sided") {
#     CI_l = original_stat - se*qt(1-alpha/2,n-1)
#     CI_u = original_stat - se*qt(alpha/2,n-1)
# 
#     CI = c(CI_l, CI_u) %>% unlist()
#     names(CI) = c(paste((alpha/2)*100,"%"), paste((1-alpha/2)*100,"%"))
#   }
#   return(CI)
# }

# percentile_ci(stat)
# student_t_ci(original_stat, stat)
# reverse_percentile_ci(original_stat, stat)

pairwise_ci = function(STAT_df, G1, G2, method, ORIG_stat) {
  orig_diff = ORIG_stat[,G1] - ORIG_stat[,G2]
  diff = (STAT_df[G1] - STAT_df[G2]) %>% t()

  if (method == "percentile") {return (percentile_ci(diff))}
  if (method == "reverse") {return(reverse_percentile_ci(orig_diff, b_stat = diff))}
  if (method == "student") {return(student_t_ci(orig_diff, b_stat = diff))}
  else {return ('method not valid')}
}

single_ci = function(STAT_df, G, method, ORIG_stat) {
  orig_diff = ORIG_stat[G]
  diff = (STAT_df[,G]) %>% t()
  
  if (method == "percentile") {return (percentile_ci(diff))}
  if (method == "reverse") {return(reverse_percentile_ci(orig_diff, b_stat = diff))}
  if (method == "student") {return(student_t_ci(orig_diff, b_stat = diff))}
  else {return ('method not valid')}
}

std.error <- function(x) sd(x)/sqrt(length(x))


cover_prob_pair = function(true_stat, origin_stat, stat_list, G1, G2, method) {
  true = true_stat[G1] - true_stat[G2]
  cover = c()
  range = c()
  for (i in 1:length(stat_list)) {
    stat_now = origin_stat[i,]
    boot_stat = stat_list[[i]]
    
    CI = pairwise_ci(boot_stat, G1, G2, method, stat_now)
    # print(CI)
    range[i] = CI[2]-CI[1]
    
    if (true >= CI[1] & true <= CI[2]) {
      cover[i] = 1
    } else {cover[i] = 0}
  }
  return(c(mean(cover), mean(range)))
}

cover_prob_single = function(true_stat, origin_stat, stat_list, G, method) {
  true = true_stat[G]
  cover = c()
  range = c()
  for (i in 1:length(stat_list)) {
    stat_now = origin_stat[i,]
    boot_stat = stat_list[[i]]
    
    CI = single_ci(boot_stat, G, method, stat_now)
    # print(CI)
    range[i] = CI[2]-CI[1]
    
    if (true >= CI[1] & true <= CI[2]) {
      cover[i] = 1
    } else {cover[i] = 0}
  }
  return(c(mean(cover), mean(range)))
}

# True_AUC = c(integrate(Survival_weibull, 0, Time_End)$value,
#              integrate(Survival_weibull, 0, Time_End)$value,
#              integrate(Survival_weibull, 0, Time_End)$value)

# cover_prob_single(true_stat=True_AUC, origin_stat = original_stat, stat_list, G=1, 'reverse')
# cover_prob_single(true_stat=True_AUC, origin_stat = original_stat, stat_list, G=1, 'percentile')
# cover_prob_single(true_stat=True_AUC, origin_stat = original_stat, stat_list, G=1, 'student')


# true_diff = c(0,0,0)
# diff_df = data.frame(diff12 = original_stat[,1] - original_stat[,2] - true_diff12,
#                      diff23 = original_stat[,2] - original_stat[,3] - true_diff23,
#                      diff13 = original_stat[,1] - original_stat[,3] - true_diff13)
# cover_prob_pair(true_stat = true_diff, origin_stat = original_stat, 
#                 stat_list = stat_list, G1 = 1, G2 = 2, 'percentile')