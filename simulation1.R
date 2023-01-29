setwd('.../OGAM')

library(foreach)
library(doParallel)

source('FNS/FNS_SmoBack.R')
source('FNS/FNS_DataGene_Simu1.R')

# define parameters
{
  R <- 100
  Kmax <- 1000
  sub_streams <- c(1,seq(20,Kmax,20))
  set.seed(2020)
  sds <- ceiling(runif(R)*1e6)
  n <- ceiling(rnorm(Kmax,500,10))
  d <- 2
  m <- 40 # No evalpoints
  eval_vec <- seq(0.05, 0.95, length.out = m)
  N <- 0
  pd1 <- 1
  Max_iter <- 50
  beta_true <- cbind(beta1_fun(eval_vec), beta2_fun(eval_vec),
                     beta1_fun_deri(eval_vec), beta2_fun_deri(eval_vec))
  link <- 'log'
  K_band <- 200 # time to stop update the constant for bandwidth 
  G <- rep(0.5,d)
  L_theta <- 5; L_sigma <- 5
}

############################ main regression ###############
######### online method
#### Input: the candidate sequence lengths
#### Output: the computing times for each update, the integrated mean squared errors and the selected bandwidths.
load('res/sim1/online_constants_for_bandwidths.Rdata')
for(L in c(3,5,10))
{
  Mcl<-100
  cl<-makeCluster(Mcl)
  registerDoParallel(cl)
  oln<- foreach(sd=sds) %dopar%{
    
    library(MASS)
    ## seed
    set.seed(sd) 
    sub_sds<-ceiling(runif(Kmax)*10^5)
    ll <- which(sds==sd)
    
    ## stored information
    time <- c()
    rss <- c()
    band <- c()
    
    for (K in 1:Kmax) {
      
      set.seed(sub_sds[K])
      # generate data
      data <- gene_data(n[K])
      X <- data$x
      y <- data$y
      N <- N + n[K]
      rm(data)
      print(paste('K=',K))

      # delta
      delta_inner <- 0.01*(K<=10)+0.001*(10<K&K<=30)+1e-04*(K>30)
      delta_outer <- delta_inner

      # online gam
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
        Ch <- C1[,min(K,K_band),ll]
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
        res_beta <- backfitting(beta0_est, beta_est,U_[,,idx[1]],V_[,,,idx[1]],Q_[,,,idx[1]],R_[,,,,idx[1]], n[K], N, h, pd1)
        if(res_beta$delta < delta_outer){
          beta0_est <- res_beta$beta0; beta_est <- res_beta$beta
          print(paste0('backfitting modified : delta=', res_beta$delta))
        }
        rm(res_beta)
      }
      # print('backfitting OK')
      
      # update statistics
      {
        res_update <- update_stats(U_,V_,Q_,R_, beta0_est, beta_est, eta, idx, n[K], N, pd1, L)
        U_ <- res_update$U_
        V_ <- res_update$V_
        Q_ <- res_update$Q_
        R_ <- res_update$R_
        rm(res_update)
      }
      t1<-Sys.time()
      
      # store
      rss1 <- sapply(1:d,function(i){var(beta_true[,i]-beta_est[,i])})
      rss <- cbind(rss,rss1)
      band<-cbind(band,h)
      time <-c(time, as.numeric(t1-t0,unit='secs'))
    }
    
    save(time, rss, band, file = paste('res/sim1/online_gam_L',L,'_',ll,'.Rdata',sep=''))
    
  }
  stopCluster(cl)
}


######### batch method
#### Input: L=1 correponds to the batch method with no candidate bandwidth
#### Output: the computing times for each update, the integrated mean squared errors and the selected bandwidths.
load('res/sim1/batch_constants_for_bandwidths.Rdata')
L <- 1
{
  Mcl<-100
  cl<-makeCluster(Mcl)
  registerDoParallel(cl)
  btch<- foreach(sd=sds) %dopar%{
    
    library(MASS)
    ## seed
    set.seed(sd) 
    sub_sds<-ceiling(runif(Kmax)*10^5)
    ll <- which(sds==sd)
    
    ## stored information
    time <- c()
    rss <- c()
    band <- c()
    X <- c(); y <- c()

    for (K in 1:Kmax) {
      
      set.seed(sub_sds[K])
      # generate data
      data <- gene_data(n[K])
      X <- rbind(X, data$x)
      y <- c(y, data$y)
      N <- N+n[K]
      rm(data)
      print(paste('K=',K))

      # delta
      delta_inner <- 0.01*(K<=10)+0.001*(10<K&K<=30)+1e-04*(K>30)
      delta_outer <- delta_inner

      if(K%in%sub_streams){
        # batch gam
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
        Ch <- C1[,min(K,K_band),ll]
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
        rss1 <- sapply(1:d,function(i){var(beta_true[,i]-beta_est[,i])})
        rss <- cbind(rss,rss1)
        band<-cbind(band,h)
        time <-c(time, as.numeric(t1-t0,unit='secs'))
      }
    }
    
    save(time, rss, band, file = paste('res/sim1/batch_gam_',ll,'.Rdata',sep=''))
    
  }
  stopCluster(cl)
}

