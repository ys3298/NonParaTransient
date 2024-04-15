reformat_df_long = function(d1) {
  n = dim(d1)[1]
  dd = NULL
  for(i in 1:n)
  {
    one_row = d1[i,]
    #################################  
    if(d1$IMC[i]==2) 
    {
      #subject, time1, time2, DFS_STATUS, IMC
      #DFS_STATUS=1 event; DFS_STATUS=0 censor;
      
      temp1 = c(one_row$id, 0, one_row$drop1, 0, 0, i)
      dd = rbind(dd, temp1)
      
      temp2 = c(one_row$id, one_row$drop1, one_row$drop2,  0, 1, i+0.1)
      dd = rbind(dd, temp2)
      
      temp3 = c(one_row$id, one_row$drop2, one_row$time, one_row$status, 2, i+0.2)
      dd = rbind(dd, temp3)
      
    }
    ####################################
    if(d1$IMC[i]==1) 
    {
      temp1 = c(one_row$id, 0, one_row$drop1, 0, 0, i)
      dd = rbind(dd, temp1)
      
      temp2 = c(one_row$id, one_row$drop1, one_row$time,  one_row$status, 1, i+0.1)
      dd = rbind(dd, temp2)
    }
    ####################################
    
    if(d1$IMC[i]==0) 
    {
      temp2 = c(one_row$id, 0, one_row$time, one_row$status, 0, i)
      dd = rbind(dd, temp2)
    }
    
  }
  
  colnames(dd) = c('id', 'time1', 'time2', 'status', 'IMC', 'id2')
  dd = as.data.frame(dd)
  dd = dd %>% dplyr::select('time1', 'time2', 'status', 'IMC', 'id2') 
  colnames(dd) = c('start_td', 'stop_td', 'endpoint_td', "covariate_td", 'patientsnumber_td')
  return(dd)
}

reformat_df_pairwise = function(d1, group1, group2) {
  n = dim(d1)[1]
  dd = NULL
  for (i in 1:n) {
    one_row = d1[i,]
    ############################################################################################################
    if (group1 == 0 & group2 == 1) {
      if(d1$IMC[i]==2) 
      {
        #subject, time1, time2, DFS_STATUS, IMC, ID
        #DFS_STATUS=1 event; DFS_STATUS=0 censor;
        
        temp1 = c(one_row$id, 0, one_row$drop1, 0, 0, i)
        dd = rbind(dd, temp1)
        
        temp2 = c(one_row$id, one_row$drop1, one_row$drop2,  0, 1, i+0.1)
        dd = rbind(dd, temp2)
      }
      ####################################
      if(d1$IMC[i]==1) 
      {
        temp1 = c(one_row$id, 0, one_row$drop1, 0, 0, i)
        dd = rbind(dd, temp1)
        
        temp2 = c(one_row$id, one_row$drop1, one_row$time,  one_row$status, 1, i+0.1)
        dd = rbind(dd, temp2)
      }
      ####################################
      
      if(d1$IMC[i]==0) 
      {
        temp2 = c(one_row$id, 0, one_row$time, one_row$status, 0, i)
        dd = rbind(dd, temp2)
      }
    }
    ############################################################################################################
    if (group1 == 1 & group2 == 2) {
      #################################  
      if(d1$IMC[i]==2) 
      {
        #subject, time1, time2, DFS_STATUS, IMC
        #DFS_STATUS=1 event; DFS_STATUS=0 censor;
        
        # temp2 = c(one_row$id, one_row$drop1,  one_row$drop2,  0, 1, i+0.1)
        temp2 = c(one_row$id, one_row$drop1 - one_row$drop1, one_row$drop2- one_row$drop1,  0, 1, i+0.1)
        dd = rbind(dd, temp2)
        
        # temp3 = c(one_row$id, one_row$drop2, one_row$time, one_row$status, 2, i+0.2)
        temp3 = c(one_row$id, one_row$drop2- one_row$drop1, one_row$time- one_row$drop1, one_row$status, 2, i+0.2)
        dd = rbind(dd, temp3)
        
      }
      ####################################
      if(d1$IMC[i]==1) 
      {
        # temp2 = c(one_row$id, one_row$drop1,  one_row$time,  one_row$status, 1, i+0.1)
        temp2 = c(one_row$id, one_row$drop1- one_row$drop1, one_row$time- one_row$drop1,  one_row$status, 1, i+0.1)
        dd = rbind(dd, temp2)
      }
    }
    ############################################################################################################  
    if (group1 == 0 & group2 == 2) {
      if(d1$IMC[i]==2) 
      {
        #subject, time1, time2, DFS_STATUS, IMC, ID
        #DFS_STATUS=1 event; DFS_STATUS=0 censor;
        
        temp1 = c(one_row$id, 0, one_row$drop1, 0, 0, i)
        dd = rbind(dd, temp1)
        
        # temp3 = c(one_row$id, one_row$drop2-(one_row$drop2-one_row$drop1), one_row$time-(one_row$drop2-one_row$drop1), one_row$status, 2, i+0.2)
        temp3 = c(one_row$id, one_row$drop2, one_row$time, one_row$status, 2, i+0.2)
        dd = rbind(dd, temp3)
        
      }
      ####################################
      if(d1$IMC[i]==1) 
      {
        temp1 = c(one_row$id, 0, one_row$drop1, 0, 0, i)
        dd = rbind(dd, temp1)
      }
      ####################################
      
      if(d1$IMC[i]==0) 
      {
        temp2 = c(one_row$id, 0, one_row$time, one_row$status, 0, i)
        dd = rbind(dd, temp2)
      }
    }
  }
  
  colnames(dd) = c('id', 'time1', 'time2', 'status', 'IMC', 'id2')
  dd = as.data.frame(dd)
  dd = dd %>% dplyr::select('id', 'time1', 'time2', 'status', 'IMC', 'id2') 
  colnames(dd) = c('id', 'start_td', 'stop_td', 'endpoint_td', "covariate_td", 'patientsnumber_td')
  return(dd)
}

reformat_df_pairwise_end = function(d1, group1, group2, time_end = 0.15) {
  n = dim(d1)[1]
  dd = NULL
  for (i in 1:n) {
    one_row = d1[i,]
    ############################################################################################################
    if (group1 == 0 & group2 == 1) {
      if(d1$IMC[i]==2) 
      {
        #subject, time1, time2, DFS_STATUS, IMC, ID
        #DFS_STATUS=1 event; DFS_STATUS=0 censor;
        
        if (one_row$drop1 - 0 > time_end) {
          temp1 = c(one_row$id, 0, time_end, 0, 0, i)
        } else {
          temp1 = c(one_row$id, 0, one_row$drop1, 0, 0, i)
        }
        dd = rbind(dd, temp1)
        
        if (one_row$drop2 - one_row$drop1 > time_end) {
          temp2 = c(one_row$id, one_row$drop1, one_row$drop1 + time_end,  0, 1, i+0.1)
        } else {
          temp2 = c(one_row$id, one_row$drop1, one_row$drop2,  0, 1, i+0.1)
        }
        dd = rbind(dd, temp2)
      }
      ####################################
      if(d1$IMC[i]==1) 
      {
        if (one_row$drop1 - 0 > time_end) {
          temp1 = c(one_row$id, 0, time_end, 0, 0, i)
        } else {
          temp1 = c(one_row$id, 0, one_row$drop1, 0, 0, i)
        }
        dd = rbind(dd, temp1)
        
        if (one_row$time - one_row$drop1 > time_end) {
          temp2 = c(one_row$id, one_row$drop1, one_row$drop1 + time_end,  0, 1, i+0.1)
        } else {
          temp2 = c(one_row$id, one_row$drop1, one_row$time,  one_row$status, 1, i+0.1)
        }
        dd = rbind(dd, temp2)
      }
      ####################################
      
      if(d1$IMC[i]==0) 
      {
        if (one_row$time - 0 > time_end) {
          temp2 = c(one_row$id, 0, 0 + time_end,  0, 0, i)
        } else {
          temp2 = c(one_row$id, 0, one_row$time,  one_row$status, 0, i)
        }
        dd = rbind(dd, temp2)
      }
    }
    ############################################################################################################
    if (group1 == 1 & group2 == 2) {
      #################################  
      if(d1$IMC[i]==2) 
      {
        #subject, time1, time2, DFS_STATUS, IMC
        #DFS_STATUS=1 event; DFS_STATUS=0 censor;
        
        # temp2 = c(one_row$id, one_row$drop1,  one_row$drop2,  0, 1, i+0.1)
        if (one_row$drop2 - one_row$drop1 > time_end) {
          temp2 = c(one_row$id, one_row$drop1 - one_row$drop1, time_end,  0, 1, i+0.1)
        } else {
          temp2 = c(one_row$id, one_row$drop1 - one_row$drop1, one_row$drop2- one_row$drop1,  0, 1, i+0.1)
        }
        dd = rbind(dd, temp2)
        
        # temp3 = c(one_row$id, one_row$drop2, one_row$time, one_row$status, 2, i+0.2)
        if (one_row$time - one_row$drop2 > time_end) {
          temp3 = c(one_row$id, one_row$drop2- one_row$drop1, one_row$drop2- one_row$drop1 + time_end, 0, 2, i+0.2)
        }
        temp3 = c(one_row$id, one_row$drop2- one_row$drop1, one_row$time- one_row$drop1, one_row$status, 2, i+0.2)
        dd = rbind(dd, temp3)
        
      }
      ####################################
      if(d1$IMC[i]==1) 
      {
        # temp2 = c(one_row$id, one_row$drop1,  one_row$time,  one_row$status, 1, i+0.1)
        if (one_row$time - one_row$drop1 > time_end) {
          temp2 = c(one_row$id, one_row$drop1- one_row$drop1, time_end,  0, 1, i+0.1)
        } else {
          temp2 = c(one_row$id, one_row$drop1- one_row$drop1, one_row$time- one_row$drop1,  one_row$status, 1, i+0.1)
        }

        dd = rbind(dd, temp2)
      }
    }
    ############################################################################################################  
    if (group1 == 0 & group2 == 2) {
      if(d1$IMC[i]==2) 
      {
        #subject, time1, time2, DFS_STATUS, IMC, ID
        #DFS_STATUS=1 event; DFS_STATUS=0 censor;
        
        if (one_row$drop1 - 0 > time_end) {
          temp1 = c(one_row$id, 0, 0+time_end, 0, 0, i)
        } else {
          temp1 = c(one_row$id, 0, one_row$drop1, 0, 0, i)
        }
        dd = rbind(dd, temp1)
        
        if (one_row$time - one_row$drop2 > time_end) {
          temp3 = c(one_row$id, one_row$drop2, one_row$drop2+time_end, 0, 2, i+0.2)
        } else {
          temp3 = c(one_row$id, one_row$drop2, one_row$time, one_row$status, 2, i+0.2)
        }
        # temp3 = c(one_row$id, one_row$drop2-(one_row$drop2-one_row$drop1), one_row$time-(one_row$drop2-one_row$drop1), one_row$status, 2, i+0.2)
        dd = rbind(dd, temp3)
        
      }
      ####################################
      if(d1$IMC[i]==1) 
      {
        if (one_row$drop1 - 0 > time_end) {
          temp1 = c(one_row$id, 0, 0+time_end, 0, 0, i)
        } else {
          temp1 = c(one_row$id, 0, one_row$drop1, 0, 0, i)
        }
        dd = rbind(dd, temp1)
      }
      ####################################
      
      if(d1$IMC[i]==0) 
      {
        if (one_row$time - 0 > time_end) {
          temp2 = c(one_row$id, 0, 0+time_end, 0, 0, i)
        } else {
          temp2 = c(one_row$id, 0, one_row$time, one_row$status, 0, i)
        }
        
        dd = rbind(dd, temp2)
      }
    }
  }
  
  colnames(dd) = c('id', 'time1', 'time2', 'status', 'IMC', 'id2')
  dd = as.data.frame(dd)
  dd = dd %>% dplyr::select('id', 'time1', 'time2', 'status', 'IMC', 'id2') 
  colnames(dd) = c('id', 'start_td', 'stop_td', 'endpoint_td', "covariate_td", 'patientsnumber_td')
  return(dd)
}

reformat_df_all = function(d1) {
  n = dim(d1)[1]
  dd = NULL
  for(i in 1:n)
  {
    one_row = d1[i,]
    #################################  
    if(d1$IMC[i]==2) 
    {
      #subject, time1, time2, DFS_STATUS, IMC
      #DFS_STATUS=1 event; DFS_STATUS=0 censor;
      
      temp1 = c(one_row$id, 0, one_row$drop1, 0, 0, i)
      dd = rbind(dd, temp1)
      
      temp2 = c(one_row$id, one_row$drop1, one_row$drop2,  0, 1, i+0.1)
      dd = rbind(dd, temp2)
      
      temp3 = c(one_row$id, one_row$drop2, one_row$time, one_row$status, 2, i+0.2)
      dd = rbind(dd, temp3)
      
    }
    ####################################
    if(d1$IMC[i]==1) 
    {
      temp1 = c(one_row$id, 0, one_row$drop1, 0, 0, i)
      dd = rbind(dd, temp1)
      
      temp2 = c(one_row$id, one_row$drop1, one_row$time,  one_row$status, 1, i+0.1)
      dd = rbind(dd, temp2)
    }
    ####################################
    
    if(d1$IMC[i]==0) 
    {
      temp2 = c(one_row$id, 0, one_row$time, one_row$status, 0, i)
      dd = rbind(dd, temp2)
    }
    
  }
  
  colnames(dd) = c('id', 'time1', 'time2', 'status', 'IMC', 'id2')
  dd = as.data.frame(dd)
  dd = dd %>% dplyr::select('id', 'time1', 'time2', 'status', 'IMC', 'id2') 
  colnames(dd) = c('id', 'start_td', 'stop_td', 'endpoint_td', "covariate_td", 'patientsnumber_td')
  return(dd)
}

Mantel.Byar.test <- function(group1 = 1, group2 = 0, Group=NULL, Event=TempTD$endpoint_td, StartTime=TempTD$start_td,	StopTime=TempTD$stop_td, method=c("SAS", "Tominaga", 'my method'), plot=0, landmark=0) {
  
  #modified from logrank test in http://aoki2.si.gunma-u.ac.jp/R/logrank.html
  #Reuire TempTD dataset created by Cox with TD variable in EZR
  
  Group.name <- Group
  if(!is.null(Group)){
    Group <- eval(parse(text=paste("TempTD$", Group, sep="")))
  } else {
    cn <- colnames(TempTD)
    len.cn <- length(cn)
    if(substring(cn[len.cn], nchar(cn[len.cn])-2, nchar(cn[len.cn]))!="_td") {
      print("Mantel.Byar() function should be done just after Cox proportional hazard modeling with time-deopendent covariate.")
    } else {
      Group.name <- cn[len.cn]
      Group <- eval(parse(text=paste("TempTD$", cn[len.cn], sep="")))
    }
  }
  
  
  method <- match.arg(method)
  data.name <- sprintf("StartTime: %s, StopTime: %s, Event: %s, Group: %s",
                       deparse(substitute(StartTime)),
                       deparse(substitute(StopTime)), deparse(substitute(Event)),
                       paste("TempTD$", Group.name, sep=""))
  OK <- complete.cases(Group, Event, StartTime, StopTime)
  Group <- Group[OK]
  Event <- Event[OK]
  StartTime <- StartTime[OK]
  StopTime <- StopTime[OK]
  
  Start <- pmin(StartTime, StopTime)					#for samples with StartTime>StopTime
  Stop <- pmax(StartTime, StopTime)
  StartTime <- Start
  StopTime <- Stop
  
  len <- length(Group)
  stopifnot(length(Event) == len, length(StopTime) == len)
  
  tg <- table(c(StopTime, rep(NA, 4)),  
              c(Group, group1, group1, group2, group2)*10+c(Event, 1, 0, 1, 0))
  k <- nrow(tg)
  nia <- table(Group)[1]
  nib <- len-nia
  na <- c(nia, (rep(nia, k)-cumsum(tg[,1]+tg[,2]))[-k]) %>% unname()	 ## transition + event loss
  nb <- c(nib, (rep(nib, k)-cumsum(tg[,3]+tg[,4]))[-k])	%>% unname() ## transition + event loss
  # na <- (rep(nia, k)-cumsum(tg[,1]+tg[,2])) %>% unname()	 ## transition + event loss
  # nb <- (rep(nib, k)-cumsum(tg[,3]+tg[,4]))	%>% unname() ## transition + event loss
  #following part is different from log-rank test
  minus <- NULL							
  for (i in 1:length(tg[,1])){
    if(as.numeric(rownames(tg))[i]==0){
      minus[i] <- sum((as.numeric(rownames(tg))[i] < StartTime)) ## remove those not in the stage 2 
    } else {
      minus[i] <- sum((as.numeric(rownames(tg))[i] <= StartTime))
    }
  }
  nb <- nb - minus
  #Following part is same with log-ranktest 
  da <- tg[,2]							
  db <- tg[,4]							
  dt <- da+db								
  nt <- na+nb								
  d <- dt/nt	# dt = A+C nt = N1+N2							
  O <- c(sum(da), sum(db))				
  ea <- na*d								
  eb <- nb*d
  ea[is.infinite(ea)]<-NA
  eb[is.infinite(eb)]<-NA
  E <- c(sum(ea, na.rm=TRUE), sum(eb, na.rm=TRUE))				
  result <- data.frame(da, db, dt, na, nb, nt, d, ea, eb)
  if (method == "Tominaga") {                   
    method <- "Mantel Byar(Tominaga)"
    chi <- sum((O-E)^2/E)
  } else if (method == "SAS") {                                      
    method <- "Mantel Byar test"
    v <- sum(dt*(nt-dt)/(nt-1)*na/nt*(1-na/nt), na.rm=TRUE)
    chi <- (sum(da, na.rm=TRUE)-sum(na*d, na.rm=TRUE))^2/v
  } else{
    method = 'my method'
    v = sum(na*nb*dt*(nt-dt) / (nt^2 *(nt-1)), na.rm=TRUE) # in the last row, when nt = 1, v = NaN
    exp = sum(na*dt/nt)
    obs = sum(da)
    chi = (abs(obs - exp) - 0.5)^2 / v
    
    # v = sum(na*nb*dt*(nt-dt) / (nt^2 *(nt-1)))
    # exp = sum(nb*dt/nt)
    # obs = sum(db)
    # chi = (abs(obs - exp) - 0.5)^2 / v
  }
  P <- pchisq(chi, 1, lower.tail=FALSE)
  if(plot>=1){							#If plot>=1, draw Simon Makuch plot with a landmark as specified.
    StartTime2 <- StartTime[StopTime>=landmark]
    StopTime2 <- StopTime[StopTime>=landmark]
    Event2 <- Event[StopTime>=landmark]
    Group2 <- Group[StopTime>=landmark]
    km <- survfit(Surv(StartTime2,StopTime2,Event2)~Group2, na.action = na.omit, conf.type="log-log")
    print(summary(km))
    #		diff <- survdiff(Surv(StopTime2,Event2)~Group2)
    n.atrisk.G1 <- NULL
    n.atrisk.G2 <- NULL
    #		n.atrisk.G1[1] <- diff$n[1]			#To correct number at risk at zero point in no event group
    #		n.atrisk.G2[1] <- 0
    len <- nchar("Group2")
    legend <- substring(names(km$strata), len+2)
    #		windows(width=7, height=7); par(lwd=1, las=1, family="sans", cex=1)
    #		dev.new()
    if (.Platform$OS.type == 'windows'){
      justDoIt(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))
    } else if (MacOSXP()==TRUE) {
      justDoIt(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))
    } else {
      justDoIt(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))
    }
    mar <- par("mar")
    mar[1] <- mar[1] + length(km$strata) + 0.5
    mar[2] <- mar[2] + 2
    par(mar=mar)
    opar <- par(mar = mar)
    on.exit(par(opar))
    #		plot(km, ylab="Probability", bty="l", col=1:32, lty=1, lwd=1, conf.int=FALSE, mark.time=TRUE)
    if(plot==1) plot(km, ylab="Probability", bty="l", col=1:32, lty=1, lwd=1, conf.int=FALSE, mark.time=TRUE)
    if(plot==2) plot(km, ylab="Probability", bty="l", col=1, lty=1:32, lwd=1, conf.int=FALSE, mark.time=TRUE)
    if(plot>=3) plot(km, ylab="Probability", bty="l", col=1, lty=1, lwd=1:32, conf.int=FALSE, mark.time=TRUE)
    xticks <- axTicks(1)
    #		n.atrisk <- nrisk(km, xticks)		#nrisk does not work properly in Simon-Makuch plot
    
    for(i in 1:length(xticks)){
      n.atrisk.G1[i] <- length(which(Group2==0 & StartTime2<=xticks[i] & xticks[i]<=StopTime2))
      n.atrisk.G2[i] <- length(which(Group2==1 & StartTime2<=xticks[i] & xticks[i]<=StopTime2))
    }
    
    n.atrisk <- rbind(n.atrisk.G1, n.atrisk.G2)
    colnames(n.atrisk) <- xticks
    
    for (i in 1:length(km$strata)){axis(1, at = xticks, labels = n.atrisk[i,], line=3+i, tick = FALSE)}
    for (i in 1:length(km$strata)){mtext(legend[i], at=-(xticks[2]-xticks[1])/2, side=1, line=4+i, cex=1)}
    title(xlab = "Number at risk", line = 3.5, adj = 0)
    #		legend ("topright", legend, col=1:32, lty=1, lwd=1,  box.lty=0, title="Time-dependent covariate")
    #		if(plot==1) legend ("topright", legend, col=1:32, lty=1, lwd=1,  box.lty=0, title="Time-dependent covariate")
    if(plot==1) legend ("topright", legend, col=1:32, lty=1, lwd=1,  box.lty=0, title=Group.name)
    if(plot==2) legend ("topright", legend, col=1, lty=1:32, lwd=1,  box.lty=0, title=Group.name)
    if(plot>=3) legend ("topright", legend, col=1, lty=1, lwd=1:32,  box.lty=0, title=Group.name)
  }
  return(structure(list(statistic=c("X-squared"=chi), parameter=c(df=1), p.value=P, 
                        method=method, data.name=data.name, result=result), class="htest"))
}


Mantel.Byar.overall <- function(Group=NULL, Event=TempTD$endpoint_td, StartTime=TempTD$start_td,	StopTime=TempTD$stop_td, method=c("SAS", "Tominaga"), plot=0, landmark=0) {
  
  #modified from logrank test in http://aoki2.si.gunma-u.ac.jp/R/logrank.html
  #Reuire TempTD dataset created by Cox with TD variable in EZR
  
  Group.name <- Group
  if(!is.null(Group)){
    Group <- eval(parse(text=paste("TempTD$", Group, sep="")))
  } else {
    cn <- colnames(TempTD)
    len.cn <- length(cn)
    if(substring(cn[len.cn], nchar(cn[len.cn])-2, nchar(cn[len.cn]))!="_td") {
      print("Mantel.Byar() function should be done just after Cox proportional hazard modeling with time-deopendent covariate.")
    } else {
      Group.name <- cn[len.cn]
      Group <- eval(parse(text=paste("TempTD$", cn[len.cn], sep="")))
    }
  }
  
  
  method <- match.arg(method)
  data.name <- sprintf("StartTime: %s, StopTime: %s, Event: %s, Group: %s",
                       deparse(substitute(StartTime)),
                       deparse(substitute(StopTime)), deparse(substitute(Event)),
                       paste("TempTD$", Group.name, sep=""))
  OK <- complete.cases(Group, Event, StartTime, StopTime)
  Group <- Group[OK]
  Event <- Event[OK]
  StartTime <- StartTime[OK]
  StopTime <- StopTime[OK]
  
  Start <- pmin(StartTime, StopTime)					#for samples with StartTime>StopTime
  Stop <- pmax(StartTime, StopTime)
  StartTime <- Start
  StopTime <- Stop
  
  StartTime1 <- Start[Group == 0 | Group == 1]
  StartTime2 <- Start[Group == 0 | Group == 2]
  
  len <- length(Group)
  stopifnot(length(Event) == len, length(StopTime) == len)
  
  tg <- table(c(StopTime, rep(NA, 6)),  
              c(Group, 0,0,1,1,2,2)*10+c(Event, 1, 0, 1, 0, 1, 0))
  k <- nrow(tg)
  nia <- table(Group)[1] %>% unname()
  nib <- table(Group)[2] %>% unname()
  nic <- table(Group)[3] %>% unname()
  na <- c(nia, (rep(nia, k)-cumsum(tg[,1]+tg[,2]))[-k])	 ## transition + event loss
  nb <- c(nib, (rep(nib, k)-cumsum(tg[,3]+tg[,4]))[-k])	 ## transition + event loss
  nc <- c(nic, (rep(nic, k)-cumsum(tg[,5]+tg[,6]))[-k])	 ## transition + event loss
  #following part is different from log-rank test
  minus1 <- NULL		
  minus2 <- NULL		
  for (i in 1:length(tg[,1])){
    if(as.numeric(rownames(tg))[i]==0){
      minus1[i] <- sum((as.numeric(rownames(tg))[i] < StartTime1))
      minus2[i] <- sum((as.numeric(rownames(tg))[i] < StartTime2))
    } else {
      minus1[i] <- sum((as.numeric(rownames(tg))[i] <= StartTime1))
      minus2[i] <- sum((as.numeric(rownames(tg))[i] <= StartTime2))
    }
  }
  nb <- nb - minus1
  nc <- nc - minus2
  #Following part is same with log-ranktest 
  da <- tg[,2]							
  db <- tg[,4]
  dc <- tg[,6]						
  dt <- da+db+dc								
  nt <- na+nb+nc
  d <- dt/nt								
  O <- c(sum(da), sum(db), sum(dc))				
  ea <- na*d								
  eb <- nb*d
  ec <- nc*d
  ea[is.infinite(ea)]<-NA
  eb[is.infinite(eb)]<-NA
  eb[is.infinite(ec)]<-NA
  E <- c(sum(ea, na.rm=TRUE), sum(eb, na.rm=TRUE), sum(ec, na.rm=TRUE) )				
  result <- data.frame(da, db, dc, dt, na, nb, nc, nt, d, ea, eb, ec)
  if (method == "Tominaga") {                   
    method <- "Mantel Byar(Tominaga)"
    chi <- sum((O-E)^2/E)
  } else {                                      
    method <- "Mantel Byar test"
    v <- sum(dt*(nt-dt)/(nt-1)*na/nt*(1-na/nt), na.rm=TRUE)
    chi <- (sum(da, na.rm=TRUE)-sum(na*d, na.rm=TRUE))^2/v
  }
  P <- pchisq(chi, 2, lower.tail=FALSE)
  if(plot>=1){							#If plot>=1, draw Simon Makuch plot with a landmark as specified.
    StartTime2 <- StartTime[StopTime>=landmark]
    StopTime2 <- StopTime[StopTime>=landmark]
    Event2 <- Event[StopTime>=landmark]
    Group2 <- Group[StopTime>=landmark]
    km <- survfit(Surv(StartTime2,StopTime2,Event2)~Group2, na.action = na.omit, conf.type="log-log")
    print(summary(km))
    #		diff <- survdiff(Surv(StopTime2,Event2)~Group2)
    n.atrisk.G1 <- NULL
    n.atrisk.G2 <- NULL
    #		n.atrisk.G1[1] <- diff$n[1]			#To correct number at risk at zero point in no event group
    #		n.atrisk.G2[1] <- 0
    len <- nchar("Group2")
    legend <- substring(names(km$strata), len+2)
    #		windows(width=7, height=7); par(lwd=1, las=1, family="sans", cex=1)
    #		dev.new()
    if (.Platform$OS.type == 'windows'){
      justDoIt(paste("windows(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))
    } else if (MacOSXP()==TRUE) {
      justDoIt(paste("quartz(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))
    } else {
      justDoIt(paste("x11(", get("window.type", envir=.GlobalEnv), "); par(", get("par.option", envir=.GlobalEnv), ")", sep=""))
    }
    mar <- par("mar")
    mar[1] <- mar[1] + length(km$strata) + 0.5
    mar[2] <- mar[2] + 2
    par(mar=mar)
    opar <- par(mar = mar)
    on.exit(par(opar))
    #		plot(km, ylab="Probability", bty="l", col=1:32, lty=1, lwd=1, conf.int=FALSE, mark.time=TRUE)
    if(plot==1) plot(km, ylab="Probability", bty="l", col=1:32, lty=1, lwd=1, conf.int=FALSE, mark.time=TRUE)
    if(plot==2) plot(km, ylab="Probability", bty="l", col=1, lty=1:32, lwd=1, conf.int=FALSE, mark.time=TRUE)
    if(plot>=3) plot(km, ylab="Probability", bty="l", col=1, lty=1, lwd=1:32, conf.int=FALSE, mark.time=TRUE)
    xticks <- axTicks(1)
    #		n.atrisk <- nrisk(km, xticks)		#nrisk does not work properly in Simon-Makuch plot
    
    for(i in 1:length(xticks)){
      n.atrisk.G1[i] <- length(which(Group2==0 & StartTime2<=xticks[i] & xticks[i]<=StopTime2))
      n.atrisk.G2[i] <- length(which(Group2==1 & StartTime2<=xticks[i] & xticks[i]<=StopTime2))
    }
    
    n.atrisk <- rbind(n.atrisk.G1, n.atrisk.G2)
    colnames(n.atrisk) <- xticks
    
    for (i in 1:length(km$strata)){axis(1, at = xticks, labels = n.atrisk[i,], line=3+i, tick = FALSE)}
    for (i in 1:length(km$strata)){mtext(legend[i], at=-(xticks[2]-xticks[1])/2, side=1, line=4+i, cex=1)}
    title(xlab = "Number at risk", line = 3.5, adj = 0)
    #		legend ("topright", legend, col=1:32, lty=1, lwd=1,  box.lty=0, title="Time-dependent covariate")
    #		if(plot==1) legend ("topright", legend, col=1:32, lty=1, lwd=1,  box.lty=0, title="Time-dependent covariate")
    if(plot==1) legend ("topright", legend, col=1:32, lty=1, lwd=1,  box.lty=0, title=Group.name)
    if(plot==2) legend ("topright", legend, col=1, lty=1:32, lwd=1,  box.lty=0, title=Group.name)
    if(plot>=3) legend ("topright", legend, col=1, lty=1, lwd=1:32,  box.lty=0, title=Group.name)
  }
  return(structure(list(statistic=c("X-squared"=chi), parameter=c(df=2), p.value=P, 
                        method=method, data.name=data.name, result=result), class="htest"))
}

