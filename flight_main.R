setwd('.../OGAM')
source('FNS/FNS_SmoBack.R')

#### parameter
{
  pd1 <- 1; pd2 <- 3
  d <- 2
  m <- 25 # No evalpoints
  eval_vec <- seq(0.05, 0.95, length.out = m)
  Max_iter <- 10
  N <- 0
  pd1 <- 1  
  K_band <- 300
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
      N <- N+n

      # delta
      delta_inner <- 0.01*(K<=10)+0.001*(10<K&K<=30)+1e-04*(K>30)
      delta_outer <- delta_inner
      
      t0<-Sys.time()
      
      if(K==1){
        
        ## some initial values
        h <- rep(1,d)
        ## stored statistics for the main regression
        U_ <- array(0, dim = c(pd1*d+1,m^d,L))
        V_ <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L)) 
        Q_ <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L))
        if(link!='identity'){
          R_ <- array(0, dim = c(pd1*d+1,pd1*d+1,pd1*d+1,m^d,L))
        }else{
          R_ <- 0
        }
      }  
      
      # compute bandwidth and candidates
      {
        Ch <- C1[,min(K,K_band)]
        h <- sapply(1:d, function(i){min(Ch[i]*N^(-1/5),h[i])})
        eta <- sapply(1:L, function(l){((L-l+1)/L)^(1/5) * h}) #dim = d*L
      }
      
      # initial parametric estimate and combination rule
      if(K==1){
        
        initial_res <- initial(X,y,h,eval_vec,pd1)
        beta0_est <- initial_res$beta0; beta_est <- initial_res$beta
        
        idx <- 1:L
        centrds <- eta
        
      }else{
        
        idx<-sapply(1:L,function(l){
          which.min(abs(eta[1,l] - centrds[1,]))
        })
        for(i in 1:d){
          centrds[i,] <- (centrds[i,idx] * (N-n[K]) + eta[i,] * n[K]) / N
        }
        
      }
      
      # backfitting
      {
        res_beta <- backfitting(beta0_est, beta_est,U_[,,idx[1]],V_[,,,idx[1]],Q_[,,,idx[1]],R_[,,,,idx[1]], n, N, h, pd1)
        if(res_beta$delta < delta_outer){
          beta0_est <- res_beta$beta0; beta_est <- res_beta$beta
          print(paste0('backfitting modified : delta=', res_beta$delta))
        }
        rm(res_beta)
      }
      # print('backfitting OK')
      
      # update statistics
      {
        res_update <- update_stats(U_,V_,Q_,R_, beta0_est, beta_est, eta, idx, n, N, pd1, L)
        U_ <- res_update$U_
        V_ <- res_update$V_
        Q_ <- res_update$Q_
        R_ <- res_update$R_
        rm(res_update)
      }
      t1<-Sys.time()      
      # store
      band<-cbind(band,h)
      time <-c(time, as.numeric(t1-t0,unit='secs'))
      beta0_store <- c(beta0_store,beta0_est)
      beta_store <- cbind(beta_store,beta_est[,1:d])

      print(paste0('K=',K,', time=', round(difftime(t1, t0, units = 'secs'),3)))
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
  load('res/flight/batch_constants_for_bandwidths.Rdata')
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
      N <- n+N
      
      if(K %in% sub_streams){

        # delta
        delta_inner <- 0.01*(K<=10)+0.001*(10<K&K<=30)+1e-04*(K>30)
        delta_outer <- delta_inner
        
        # online gam
        t0 <- Sys.time()
        {
          ## some initial values
          h <- rep(1,d)
          ## stored statistics for the main regression
          U_ <- array(0, dim = c(pd1*d+1,m^d,L))
          V_ <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L)) 
          Q_ <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L))
          if(link!='identity'){
            R_ <- array(0, dim = c(pd1*d+1,pd1*d+1,pd1*d+1,m^d,L))
          }else{
            R_ <- 0
          }
        }  
        
        # compute bandwidth and candidates
        Ch <- C1[,min(K,K_band)]
        h <- Ch*N^(-1/5)
        
        # initial parametric estimate and combination rule
        initial_res <- initial(X,y,h,eval_vec,pd1)
        beta0_est <- initial_res$beta0; beta_est <- initial_res$beta
        
        # backfitting
        {
          res_beta <- backfitting(beta0_est, beta_est,U_[,,1],V_[,,,1],Q_[,,,1],R_[,,,,1], N, N, h, pd1)
          if(res_beta$delta < delta_outer){
            beta0_est <- res_beta$beta0; beta_est <- res_beta$beta
            print(paste0('backfitting modified : delta=', res_beta$delta))
          }
          rm(res_beta)
        }
        t1 <- Sys.time()
        
        # store
        band<-cbind(band,h)
        time <-c(time, as.numeric(t1-t0,unit='secs'))
        beta0_store <- c(beta0_store,beta0_est)
        beta_store <- cbind(beta_store,beta_est[,1:d])
        
      }

    }
    
    start <- Kmax

  }
  
  save(time, beta0_store, beta_store,band,
       file = 'res/flight_full_batch.Rdata')
}

