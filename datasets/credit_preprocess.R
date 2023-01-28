setwd(".../OGAM")
library(MASS)

if(!file.exists('datasets/credit/P2P_Macro_Data_processed.Rdata')){
    
    library(haven)
    
    data <- read_dta('datasets/credit/P2P_Macro_Data.dta')
    data <- data[,c('badloan','loan_amnt','int_rate','dti','annual_inc')]
    data <- data[data$loan_status!="Current",]
    data <- data[,c('badloan','loan_amnt','int_rate','dti','annual_inc')]
    
    gc()

    data$loan_amnt <- sqrt(data$loan_amnt)
    data <- data[which(data$loan_amnt<190),]
    data$loan_amnt <- (data$loan_amnt-min(data$loan_amnt))/(max(data$loan_amnt)-min(data$loan_amnt))
    hist(data$loan_amnt[which(data$loan_amnt<0.3)])
    
    data$int_rate <- sqrt(data$int_rate)
    data <- data[which(data$int_rate<0.52),]
    data$int_rate <- (data$int_rate-min(data$int_rate))/(max(data$int_rate)-min(data$int_rate))
    hist(data$int_rate[which(data$int_rate<0.3)])
    
    data$dti <- sqrt(data$dti)
    data <- data[which(data$dti<6.2&data$dti>1.8),]
    data$dti <- (data$dti-min(data$dti))/(max(data$dti)-min(data$dti))
    hist(data$dti[which(data$dti<0.3)])
    
    data$annual_inc <- log(data$annual_inc)
    data <- data[which(data$annual_inc>10&data$annual_inc<12.3),]
    data$annual_inc <- (data$annual_inc-min(data$annual_inc))/(max(data$annual_inc)-min(data$annual_inc))
    hist(data$annual_inc[which(data$annual_inc<0.3)])
    
    for(i in 1:4){
      data[,i+1] <- (data[,i+1]-min(data[,i+1]))/(max(data[,i+1])-min(data[,i+1]))
    }
    
    save(data,file = 'P2P_Macro_Data_processed.Rdata')
  }
