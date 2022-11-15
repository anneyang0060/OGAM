setwd('/home/yangy/OGAM_code_and_data')

library(dplyr)

rm(list=ls())
FlightNum_old <- c()
OccurNum_old <- c()
DelayNum_old <- c()
Nt <- 30

for(year in 1987:2005){
  
  if(file.exists(paste0('datasets/flight/',year,'_gam.Rdata'))){next}
  df <- read.csv(paste('datasets/flight/',year,'.csv',sep=''))
  
  df <- df[(!is.na(df$UniqueCarrier) & !is.na(df$Origin)),]
  df$FlightNum <- paste(df$Origin, df$FlightNum, sep='')
  df <- df[(!is.na(df$DepDelay)),]
  df$Delayed <- (df$DepDelay > 15)
  
  df <- df[,c('Year', 'Month', 'DayofMonth', 'DayOfWeek', 'FlightNum',  'CRSDepTime', 'DepDelay', 'Delayed')]
  for(i in 1:ncol(df)){df <- df[(!is.na(df[,i])), ]}
  time <- (df$Month*100 + df$DayofMonth)*1e4 + df$CRSDepTime
  df <- df[order(time),]
  rm(time)
  
  FlightNum_current_year <- unique(df$FlightNum)
  new_idx <- c()
  for(i in 1:length(FlightNum_current_year)){
    if(!(FlightNum_current_year[i]%in%FlightNum_old)){
      new_idx <- c(new_idx, i)
    }
  }
  FlightNum <- c(FlightNum_old, FlightNum_current_year[new_idx])
  OccurNum <- c(OccurNum_old, rep(0, length(new_idx)))
  DelayNum <- c(DelayNum_old, rep(0, length(new_idx)))
  
  ## compute historical delay rate ###
  df$HistDelayRate <- 0
  df$HistOccurNum <- 0

  for(i in 1:length(FlightNum_current_year)){
    
    df_idx <- which(df$FlightNum==FlightNum_current_year[i])
    stat_idx <- which(FlightNum==FlightNum_current_year[i])
    hist_occurnum <- (1:length(df_idx)) + OccurNum[stat_idx] - 1
    hist_delaynum <- c(0,cumsum(df[df_idx, 'Delayed'])[-length(df_idx)]) + DelayNum[stat_idx]
    df[df_idx, 'HistOccurNum'] <- hist_occurnum
    df[df_idx, 'HistDelayRate'] <- hist_delaynum / hist_occurnum
    OccurNum[stat_idx] <- OccurNum[stat_idx] + length(df_idx)
    DelayNum[stat_idx] <- DelayNum[stat_idx] + sum(df[df_idx, 'Delayed'])
    print(paste(year, round(100*i/length(FlightNum_current_year),2), sep=', '))
    
  }

  ####### organize
  df <- df[(!is.na(df$HistDelayRate)),]
  
  df$Origin <- substr(df$FlightNum,start = 1,stop = 2)
  df$subj <- paste(df$Month*100 + df$DayofMonth, df$Origin)
  group <- group_by(df,subj)
  smr <- summarise(group, N_obs = n(), N_delay = sum(Delayed))
  smr$DailyDelayRate <- smr$N_delay/smr$N_obs
  
  df <- df[df$HistOccurNum>Nt,]
  df <- df[df$CRSDepTime>=600 & df$CRSDepTime<=2300, ]
  
  df$CRSDepTime <- (df$CRSDepTime%/%100) * 60 + df$CRSDepTime%%100
  df$CRSDepTime <- (df$CRSDepTime - 6*60)/(23*60 - 6*60)
  
  df<-df[df$HistDelayRate>=0.15,]
  df<-df[df$HistDelayRate<=0.65,]
  df$HistDelayRate=(df$HistDelayRate-0.14)/0.1
  df$HistDelayRate=log(df$HistDelayRate)
  df$HistDelayRate <- (df$HistDelayRate-log((0.15-0.14)/0.1))/(log((0.65-0.14)/0.1)-log((0.15-0.14)/0.1))
  
  df$Date <- df$Month*100+df$DayofMonth

  save(df,file = paste0('datasets/flight/',year,'_gam.Rdata'))
  
  print(year)
  
  FlightNum_old <- FlightNum
  OccurNum_old <- OccurNum
  DelayNum_old <- DelayNum
  
}
