setwd('/home/yangy/OGAM_code_and_data')

library(foreach)
library(doParallel)

source('FNS/FNS_SmoBack.R')
source('FNS/FNS_DataGene_Simu1.R')

# parameters
{
  R <- 100
  Kmax <- 1000
  sub_streams <- c(1,seq(20,Kmax,20))
  set.seed(2020)
  sds <- ceiling(runif(R)*1e6)
  n <- ceiling(rnorm(Kmax,500,10))
  d <- 2
  m <- 40 # No evalpoints

  Max_iter <- 50
  beta_true <- cbind(beta1_fun(eval_vec), beta2_fun(eval_vec),
                     beta1_fun_deri(eval_vec), beta2_fun_deri(eval_vec))
  link <- 'log'
  
  K_band <- 200 # time to stop update the constant for bandwidth 
  G <- rep(0.5,d)
  L_theta <- 5; L_sigma <- 5
}

############################ main regression ###############
######### online
load('res/sim1/online_constants_for_bandwidths.Rdata')
band_select <- FALSE
for(L in c(3,5,10))
{
  Mcl<-50
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
      rm(data)
      print(paste('K=',K))

      # delta
      delta_inner <- 0.01*(K<=10)+0.001*(10<K&K<=30)+1e-04*(K>30)
      delta_outer <- delta_inner

      # online gam
      t0<-Sys.time()
      ogam(K, X, y, n, m, delta_inner, delta_outer, Max_iter, band_select, K_band, C1=C1[,min(K,K_band)], L=L)
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


######### batch
load('res/sim1/batch_constants_for_bandwidths.Rdata')
band_select <- FALSE
L <- 1
{
  Mcl<-50
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
      rm(data)
      print(paste('K=',K))

      # delta
      delta_inner <- 0.01*(K<=10)+0.001*(10<K&K<=30)+1e-04*(K>30)
      delta_outer <- delta_inner

      if(K%in%sub_streams){
        # batch gam
        t0<-Sys.time()
        ogam(1, X, y, n, delta_inner, delta_outer, band_select, K_band, C1=C1[,which(sub_streams==min(K,K_band))])
        t1<-Sys.time()
        
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

