setwd('/home/yangy/OGAM_code_and_data')
source('FNS/FNS_SmoBack.R')

#### parameter
{
  pd1 <- 1; pd2 <- 3
  d <- 2
  m <- 25 # No evalpoints
  eval_vec <- seq(0.05, 0.95, length.out = m)
  
  Max_iter <- 10
  
  K_band <- 600
  G <- rep(1,d)
  L_theta <- 5; L_sigma <- 5
  L <- 10
  link <- 'logit'
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
  load('res/flight/online_constants_for_bandwidths.Rdata')
  start <- 0
  
  for(year in 1996:2004){
    
    { 
      load(paste0('datasets/flight/', year, '_gam.Rdata'))
      tab <- table(df$Date)
      dates <- as.numeric(names(tab))[tab>=100]
      if(year>2000){dates <- as.numeric(names(tab))[tab>=4000]}
      Kmax <- length(dates) + start
    }
    
    for (K in (start+1):Kmax) {
      
      # generate data
      date <- dates[(K-start)]
      data <- df[(df$Date == date), ]
      X <- cbind(data$CRSDepTime,data$HistDelayRate)
      y <- data$Delayed
      n <- length(y)

      # delta
      delta_inner <- 0.01*(K<=10)+0.001*(10<K&K<=30)+1e-04*(K>30)
      delta_outer <- delta_inner
      
      tt0 <- Sys.time()
      # online gam
      ogam(K, X, y, n, m, delta_inner, delta_outer, Max_iter, band_select, K_band, C1=C1[,min(K,K_band)], L=L)
      tt1 <- Sys.time()
      
      # store
      band<-cbind(band,h)
      time <-c(time, as.numeric(tt1-tt0,unit='secs'))
      beta0_store <- c(beta0_store,beta0_est)
      beta_store <- cbind(beta_store,beta_est[,1:d])

      print(paste0('K=',K,', time=', round(difftime(tt1, tt0, units = 'secs'),3)))
      par(mfrow=c(1,2))
      plot(beta_store[,2*K-1],type='l')
      plot(beta_store[,2*K], type='l')
    }
    
    start <- Kmax

  }
  
  save(time, beta0_store, beta_store, band,
       file = 'res/flight/flight_full_online1.Rdata')
}

# batch
sub_streams <- c(1,seq(10,50,10),seq(100,1000,100), seq(1500,2500,500),3283)
{
  # stored information
  time <- c()
  beta0_store <-c()
  beta_store <-c()
  band <- c()
  X <- c(); y <- c()
  band_select <- FALSE
  load('res/flight/online_constants_for_bandwidths.Rdata')
  start <- 0
  
  for(year in 1996:2004){
    
    { 
      load(paste0('datasets/flight/', year, '_gam.Rdata'))
      tab <- table(df$Date)
      dates <- as.numeric(names(tab))[tab>=100]
      if(year>2000){dates <- as.numeric(names(tab))[tab>=4000]}
      Kmax <- length(dates) + start
    }
    
    for (K in (start+1):Kmax) {
      
      # generate data
      date <- dates[(K-start)]
      data <- df[(df$Date == date), ]
      X <- rbind(X, cbind(data$CRSDepTime,data$HistDelayRate))
      y <- c(y, data$Delayed)
      n <- length(y)
      
      if(K %in% sub_streams){

        # delta
        delta_inner <- 0.01*(K<=10)+0.001*(10<K&K<=30)+1e-04*(K>30)
        delta_outer <- delta_inner
        
        # online gam
        tt0 <- Sys.time()
        gam(1, X, y, n, delta_inner, delta_outer, band_select, K_band, C1=C1[,which(sub_streams==min(K,K_band))])
        tt1 <- Sys.time()
        
        # store
        band<-cbind(band,h)
        time <-c(time, as.numeric(tt1-tt0,unit='secs'))
        beta0_store <- c(beta0_store,beta0_est)
        beta_store <- cbind(beta_store,beta_est[,1:d])
        
      }

    }
    
    start <- Kmax

  }
  
  save(time, beta0_store, beta_store,band,
       file = 'res/flight_full_batch.Rdata')
}

