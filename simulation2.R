setwd('/home/yangy/gam')

library(foreach)
library(doParallel)

# load related functions
source('FNS/FNS_SmoBack.R')
source('FNS/FNS_DataGene_Simu2.R')

### define input parameters
#' @param R number of simulated replicates
#' @param Kmax the total number of blocks
#' @param sub_stream time to conduct batch estimate
#' @param sds the seeds indices which are generated randomly
#' @param n the sample size of each data block
#' @param d the model dimension
#' @param m number of evaluation points
#' @param pd1=1 corresponds to local linear smoothing
#' @param Max_iter the maximal iteration steps for the algorithm
#' @param beta_true values of underlying true function at the evaluation points
#' @param link the link function for the generalized additive model
#' @param K_band time to stop update the constant for bandwidth
{
  R <- 100
  Kmax <- 600
  sub_streams <- c(1,seq(20,Kmax,20))
  set.seed(2020)
  n <- ceiling(rnorm(Kmax,500,10)); n[1] <- 2000
  d <- 4
  m <- 10 
  eval_vec <- seq(0.05, 0.95, length.out = m)
  pd1 <- 1
  Max_iter <- 50
  N <- 0
  beta_true <- cbind(beta1_fun(eval_vec), beta2_fun(eval_vec), 
                     beta3_fun(eval_vec), beta4_fun(eval_vec))
  link <- 'identity'
  K_band=200
}

############################### main regression #########################
######### online method
#### Input: 
## L: the candidate sequence lengths
## N: initialize the accumulated sample size
#### Output: 
## time: the computing times for each update
## rss: the integrated mean squared errors
## band: the selected bandwidthsload('res/sim2/online_constants_for_bandwidths.Rdata')
N <- 0
for(L in c(3,5, 10))
{
  Mcl<-100
  cl<-makeCluster(Mcl)
  registerDoParallel(cl)
  oln<- foreach(sd=1:100) %dopar%{
    
    library(MASS)
    set.seed(sd) 
    
    ## stored information
    time <- c()
    rss <- c()
    band <- c()
    N <- 0
    
    for (K in 1:Kmax) {
      
      # delta_stop
      {
        delta_inner <- 0.01*(K<=10)+0.001*(10<K&K<=30)+1e-04*(K>30)
        delta_outer <- delta_inner
      }
      
      # generate data
      {
        data <- gene_data(n[K])
        X <- data$x
        y <- data$y
        N <- N + n[K]
        rm(data)
      }  
      
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
        Ch <- C1[,min(K,K_band),sd]
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
      band <- cbind(band,h)
      time <- c(time, as.numeric(t1-t0,unit = 'secs'))
    }
    
    save(time, rss, file = paste0('res/sim2/online_gam_L',L,'_',sd,'.Rdata'))    
  }
  stopCluster(cl)
}

######### batch method
#### Input: 
## L: the candidate sequence lengths
## N: initialize the accumulated sample size
#### Output: 
## time: the computing times for each update
## rss: the integrated mean squared errors
## band: the selected bandwidthsload('res/sim2/batch_constants_for_bandwidths.Rdata')
L <- 1
N <- 0
{
  
  Mcl<-100
  cl<-makeCluster(Mcl)
  registerDoParallel(cl)
  btch <- foreach(sd=1:100) %dopar%{
    
    library(MASS)
    set.seed(sd) 
    
    ## stored information
    time <- c()
    rss <- c()
    band <- c()
    X <- c(); y <- c()
    
    for (K in 1:Kmax) {
      
      # delta_stop
      {
        delta_inner <- 0.01*(K<=10)+0.001*(10<K&K<=30)+1e-04*(K>30)
        delta_outer <- delta_inner
      }
      
      # generate data
      {
        data <- gene_data(n[K])
        X <- rbind(X, data$x)
        y <- c(y,data$y)
        N <- N + n[K]
        rm(data)
      }  
      
      
      if(K%in%sub_streams){
        
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
        Ch <- C1[,min(K,K_band),sd]
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
        
        print(paste0('K=',K,'time=',as.numeric(t1-t0,unit = 'secs')))
      }
      
    }
    
    save(time, rss, file = paste0('res/sim2/batch_gam_',sd,'_0.Rdata'))
    
  }
  stopCluster(cl)
}
