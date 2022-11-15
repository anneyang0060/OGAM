setwd('/home/yangy/OGAM_code_and_data')

library(foreach)
library(doParallel)

source('FNS/FNS_SmoBack.R')
source('FNS/FNS_DataGene_Simu2.R')

# parameters
{
  R <- 100
  Kmax <- 600
  sub_streams <- c(1,seq(20,Kmax,20))
  set.seed(2020)
  sds <- ceiling(runif(R)*1e6)
  n <- ceiling(rnorm(Kmax,500,10))
  n[1] <- 2000
  d <- 4
  m <- 10 # No evalpoints

  Max_iter <- 50
  N <- 0
  beta_true <- cbind(beta1_fun(eval_vec), beta2_fun(eval_vec), beta3_fun(eval_vec), 
                     beta4_fun(eval_vec), beta5_fun(eval_vec))
  link <- 'identity'
  
  K_band=200
  G <- rep(0.5,d)
  L_theta <- 3; L_sigma <- 3
}

############################### main regression #########################
######### online 
load('res/sim2/online_constants_for_bandwidths.Rdata')
band_select <- FALSE
for(L in c(3,5,10))
{
  
  Mcl<-100
  cl<-makeCluster(Mcl)
  registerDoParallel(cl)
  oln<- foreach(sd=sds) %dopar%{

    library(MASS)
    set.seed(sd) 
    sub_sds<-ceiling(runif(Kmax)*10^5)
    ll <- which(sds==sd)
    
    ## stored information
    time <- c()
    rss <- c()
    band <- c()
    
    for (K in 1:Kmax) {
      
      set.seed(sub_sds[K])
      
      # delta_stop
      {
        delta_stop_inner <- 0.01*(K<=10)+0.001*(10<K&K<=30)+1e-04*(K>30)
        delta_stop_outer <- delta_stop_inner
      }
      
      # generate data
      {
        data <- gene_data(n[K])
        X <- data$x
        y <- data$y
        N <- N + n[K]
        rm(data)
      }  
      
      # online gam
      t0<-Sys.time()
      ogam(K, X, y, n, m, delta_inner, delta_outer, Max_iter, band_select, K_band, C1=C1[,min(K,K_band)], L=L)
      t1<-Sys.time()
      
      # store
      rss1 <- sapply(1:d,function(i){var(beta_true[,i]-beta_est[,i])})
      rss <- cbind(rss,rss1)
      band <- cbind(band,h)
      time <- c(time, as.numeric(t1-t0,unit = 'secs'))
    }
    
    save(time, rss, band, file = paste0('res/sim2/online_gam_L',L,'_',ll,'.Rdata'))
    
  }
  stopCluster(cl)
}

######### batch 
load('res/sim2/batch_constants_for_bandwidths.Rdata')
band_select <- FALSE
L <- 1
{
  
  Mcl<-100
  cl<-makeCluster(Mcl)
  registerDoParallel(cl)
  btch <- foreach(sd=sds) %dopar%{
    
    library(MASS)
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
      
      # delta_stop
      {
        delta_stop_inner <- 0.01*(K<=10)+0.001*(10<K&K<=30)+1e-04*(K>30)
        delta_stop_outer <- delta_stop_inner
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
        # batch gam
        t0<-Sys.time()
        gam(1, X, y, n, delta_inner, delta_outer, band_select, K_band, C1=C1[,which(sub_streams==min(K,K_band))])
        t1<-Sys.time()
        
        # store
        rss1 <- sapply(1:d,function(i){var(beta_true[,i]-beta_est[,i])})
        rss <- cbind(rss,rss1)
        band<-cbind(band,h)
        time <-c(time, as.numeric(t1-t0,unit='secs'))

        print(paste0('K=',K,'time=',as.numeric(t1-t0,unit = 'secs')))
      }
      
    }
    
    save(time, rss, band, file = paste0('res/sim2/batch_gam_',ll,'.Rdata'))
    
  }
  stopCluster(cl)
}
