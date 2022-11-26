setwd(".../OGAM")
library(MASS)
source('FNS/FNS_SmoBack_credit.R')

# read in data
{
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
  
  load('datasets/credit/P2P_Macro_Data_processed.Rdata')
  Nfull <- nrow(data)
  Ntrain <- round(Nfull*0.9)
  set.seed(111)
  idx <- sample(1:Nfull,Ntrain)
  ytrain <- data$badloan[idx]
  Xtrain <- as.matrix(data[idx,c('loan_amnt','int_rate','dti','annual_inc')])
  ytest <- data$badloan[(1:Nfull)[-idx]]
  Xtest <- as.matrix(data[(1:Nfull)[-idx],c('loan_amnt','int_rate','dti','annual_inc')])
  rm(idx)
}

#### parameter
{
  d <- 4
  m <- 10 # No evalpoints
  Max_iter <- 50
  K_band <- 200
  L <- 5
  link <- 'logit'
  nn <- rep(1000,822)
  Kmax <- length(nn)
  nn[Kmax] <- Ntrain-sum(nn[1:(Kmax-1)])
}

#### estimate
# online
{
  # stored information
  time <- c()
  beta0_store <-c()
  beta_store <-c()
  band <- c()
  band_select <- FALSE
  load('res/credit/online_constants_for_bandwidths.Rdata')
  NN <- 0

  for (K in 1:Kmax) {
    
    # generate data
    X <- as.matrix(Xtrain[(NN+1):(NN+nn[K]),])
    y <- ytrain[(NN+1):(NN+nn[K])]
    NN <- NN + nn[K]
    
    # delta
    delta_inner <- 0.01*(K<=10)+0.001*(10<K&K<=50)+1e-04*(K>50)
    delta_outer <- delta_inner*10

    t0<-Sys.time()
    ogam(K, X, y, nn[K], m, delta_inner, delta_outer, Max_iter, band_select, K_band, C1[,min(K,K_band)], L)
    t1<-Sys.time()
    
    # store
    band<-cbind(band,h)
    time <-c(time, as.numeric(t1-t0,unit='secs'))
    beta0_store <- c(beta0_store,beta0_est)
    beta_store <- cbind(beta_store,beta_est[,1:d])
    
  }
  
  save(time, beta0_store, beta_store, band,
       file = 'res/credit/credit_online.Rdata')
}

# batch
{
  # stored information
  time <- c()
  beta0_store <-c()
  beta_store <-c()
  band <- c()
  band_select <- FALSE
  load('res/credit/batch_constants_for_bandwidths.Rdata')
  NN <- 0
  X <- c(); y <- c()
  sub_streams <- c(1, seq(20,180,20), seq(200,700,100), 822)
  
  for (K in 1:Kmax) {
    
    # generate data
    X <- rbind(X, as.matrix(Xtrain[(NN+1):(NN+nn[K]),]))
    y <- c(y, ytrain[(NN+1):(NN+nn[K])])
    
    if(K%in%sub_streams){
      
      # delta
      delta_inner <- 0.01*(K<=10)+0.001*(10<K&K<=50)+1e-04*(K>50)
      delta_outer <- delta_inner*10
      
      tt0 <- Sys.time()
      # gam
      gam(1, X, y, length(y), m, delta_inner, delta_outer, Max_iter, band_select, K_band, 
           C1=C1[,min(which(sub_streams==K),11)])
      tt1 <- Sys.time()
      
      # store
      band<-cbind(band,h)
      time <-c(time, as.numeric(tt1-tt0,unit='secs'))
      beta0_store <- c(beta0_store,beta0_est)
      beta_store <- cbind(beta_store,beta_est[,1:d])
      
      print(paste0('K=',K,', m=',m,', time=', round(difftime(tt1, tt0, units = 'secs'),3)))
      par(mfrow=c(2,2))
      plot(beta_est[,1],type='l',main=paste0('K=',K)); plot(beta_est[,2],type='l')
      plot(beta_est[,3], type='l'); plot(beta_est[,4],type='l')
      
      save(time, beta0_store, beta_store, band,
           file = paste0('res/credit/credit_batch_K',K,'.Rdata'))

    }
    
  }
  
  save(time, beta0_store, beta_store, band,
       file = 'res/credit/credit_full_batch.Rdata')
}
