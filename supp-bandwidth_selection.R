setwd('/home/yangy/gam')

library(foreach)
library(doParallel)

############################### bandwidth selection for Simulation1 #########################
rm(list = ls())
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
  pd1 <- 1; pd2 <- 3
  d <- 2
  m <- 40 # No evalpoints
  eval_vec <- seq(0.05, 0.95, length.out = m)
  
  Max_iter <- 50
  beta_true <- cbind(beta1_fun(eval_vec), beta2_fun(eval_vec),
                     beta1_fun_deri(eval_vec), beta2_fun_deri(eval_vec))
  link <- 'log'
  
  # beta0 <- 0
  # beta <- matrix(0, m, (pd1+1)*d) # beta1,beta2,beta11,beta12
  # h <- rep(1,d)
  # load('constants for bandwidths.Rdata')
  K_band <- 200
  G <- rep(0.5,d)
  L_theta <- 5; L_sigma <- 5
}

######### online band select
{
  {
    U_theta <- array(0, dim = c(pd2*d+1,m^d,L_theta))
    V_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,m^d,L_theta)) 
    Q_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,m^d,L_theta))
    R_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,pd2*d+1,m^d,L_theta))
    U_sigma <- array(0, dim = c(pd1*d+1,m^d,L_sigma))
    V_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L_sigma)) 
    Q_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L_sigma))
    R_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,pd1*d+1,m^d,L_sigma))
    C <- c(); theta_store <- c(); sigma_store <- c()
  }  
  
  Mcl<-100
  cl<-makeCluster(Mcl)
  registerDoParallel(cl)
  oln<- foreach(sd=sds) %dopar%{
    
    library(MASS)
    set.seed(sd) 
    sub_sds<-ceiling(runif(Kmax)*10^5)
    ll <- which(sds==sd)
    
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
        print(paste('K=',K))
      }  
      
      # bandwidth selection
      {
        
        delta_stop_inner_sigma <- delta_stop_inner*10
        delta_stop_outer_sigma <- delta_stop_inner_sigma
        delta_stop_inner_theta <- delta_stop_inner_sigma*10
        delta_stop_outer_theta <- delta_stop_inner_theta
        
        h_theta <- G * N^(-1/7); h_sigma <- G * N^(-1/5)
        eta_theta <- sapply(1:L_theta, function(l){((L_theta-l+1)/L_theta)^(1/7) * h_theta})#dim: d*L
        eta_sigma <- sapply(1:L_sigma, function(l){((L_sigma-l+1)/L_sigma)^(1/5) * h_sigma})
        
        if(K==1){
          
          # generate initial estimate and bandwidth
          initial_res <- initial(X,y,h_theta,eval_vec,pd2)
          beta0_theta_old <- initial_res$beta0; beta_theta_old <- initial_res$beta
          beta0_theta <- initial_res$beta0; beta_theta <- initial_res$beta
          idx_theta <- 1:L_theta
          centrds_theta <- eta_theta
          rm(initial_res)
          
          initial_res <- initial(X,y,h_sigma,eval_vec,pd1)
          beta0_sigma_old <- initial_res$beta0; beta_sigma_old <- initial_res$beta
          beta0_sigma <- initial_res$beta0; beta_sigma <- initial_res$beta
          idx_sigma <- 1:L_sigma
          centrds_sigma <- eta_sigma
          rm(initial_res)
          
        }else{
          
          # update index for pseudo-bandwidth and the centroids
          idx_theta<-sapply(1:L_theta,function(l){
            which.min(abs(eta_theta[1,l] - centrds_theta[1,]))
          })# idx of all d dimensions are the same, it is sufficient to compuute the first dim
          for(i in 1:d){
            centrds_theta[i,] <- (centrds_theta[i,idx_theta] * (N-n[K]) + eta_theta[i,] * n[K]) / N
          }
          
          idx_sigma<-sapply(1:L_sigma,function(l){
            which.min(abs(eta_sigma[1,l] - centrds_sigma[1,]))
          })
          for(i in 1:d){
            centrds_sigma[i,] <- (centrds_sigma[i,idx_sigma] * (N-n[K]) + eta_sigma[i,] * n[K]) / N
          }
          
        }
        
        # backfitting
        {
          res_beta <- backfitting(beta0_theta, beta_theta,U_theta[,,idx_theta[1]],V_theta[,,,idx_theta[1]],
                                  Q_theta[,,,idx_theta[1]],R_theta[,,,,idx_theta[1]], n[K], N, h_theta, pd2)
          if(res_beta$delta < delta_stop_outer_theta){
            beta0_theta <- res_beta$beta0
            beta_theta <- res_beta$beta
          }
          rm(res_beta)
          
          res_beta <- backfitting(beta0_sigma, beta_sigma,U_sigma[,,idx_sigma[1]],V_sigma[,,,idx_sigma[1]],
                                  Q_sigma[,,,idx_sigma[1]],R_sigma[,,,,idx_sigma[1]], n[K], N, h_sigma, pd1)
          if(res_beta$delta < delta_stop_outer_sigma){
            beta0_sigma <- res_beta$beta0
            beta_sigma <- res_beta$beta
          }
          rm(res_beta)
        }
        
        # compute the constants
        {
          deri <- beta_theta[,(2*d+1):(3*d)]^2
          beta_hat <- matrix(0,nrow(X),d); sigma2 <- rep(0,d)
          for(i in 1:d){
            beta_hat[,i] <- predict(smooth.spline(eval_vec,beta_sigma[,i]),X[,i])$y
          }
          # sigma2 <- mean(exp(rowSums(beta_hat[,(1:d)])+beta0_sigma))
          for(i in 1:d){
            sigma2[i] <- mean(exp(beta_hat[,(1:d)[-i]]))*mean(exp(-beta_sigma[,i]))*exp(-beta0_sigma)
          }
          
          Ch <- (15 * sigma2/colMeans(deri))^(1/5)
          C <- cbind(C, Ch)
          theta_store <- cbind(theta_store, colMeans(deri))
          sigma_store <- cbind(sigma_store, sigma2) 
        }
        
        # update statistics
        {
          res_update <- update_stats(U_theta,V_theta,Q_theta,R_theta, beta0_theta,beta_theta,eta_theta,
                                     idx_theta,n[K],N,pd2,L_theta)
          U_theta <- res_update$U_
          V_theta <- res_update$V_
          Q_theta <- res_update$Q_
          R_theta <- res_update$R_
          rm(res_update)
          
          res_update <- update_stats(U_sigma,V_sigma,Q_sigma,R_sigma, beta0_sigma,beta_sigma,eta_sigma,
                                     idx_sigma,n[K],N,pd1,L_sigma)
          U_sigma <- res_update$U_
          V_sigma <- res_update$V_
          Q_sigma <- res_update$Q_
          R_sigma <- res_update$R_
          rm(res_update)
        }
        
      }
      
      print(paste('K=',K))
    }
    
    save(C, theta_store, sigma_store, 
         file = paste('res/sim1/online_constant_for_bandwidth_sd',ll,'.Rdata',sep=''))
    
  }
  stopCluster(cl)
  
  C1 <- array(0, dim=c(d,Kmax,R))
  for(ll in 1:R){
    load(paste('res/sim1/online_constant_for_bandwidth_sd',ll,'.Rdata',sep=''))
    C1[,,ll] <- C
    file.remove(paste('res/sim1/online_constant_for_bandwidth_sd',ll,'.Rdata',sep=''))
  }
  save(C1, file='res/sim1/online_constants_for_bandwidths.Rdata')
  
}

######### batch band select
{
  {
    U_theta <- array(0, dim = c(pd2*d+1,m^d))
    V_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,m^d)) 
    Q_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,m^d))
    R_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,pd2*d+1,m^d))
    U_sigma <- array(0, dim = c(pd1*d+1,m^d))
    V_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d)) 
    Q_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d))
    R_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,pd1*d+1,m^d))
    X <- c(); y <- c(); N <- 0
    C<- c(); theta_store <- c(); sigma_store <- c()
  }  
  
  Mcl<-100
  cl<-makeCluster(Mcl)
  registerDoParallel(cl)
  oln<- foreach(sd=sds) %dopar%{
    
    library(MASS)
    set.seed(sd) 
    sub_sds<-ceiling(runif(Kmax)*10^5)
    ll <- which(sds==sd)
    
    
    for (K in 1:500) {
      
      set.seed(sub_sds[K])
      
      # delta_stop
      {
        delta_stop_inner <- 0.01*(K<=10)+0.001*(10<K&K<=30)+1e-04*(K>30)
        delta_stop_outer <- delta_stop_inner
        delta_stop_inner_sigma <- delta_stop_inner*10
        delta_stop_outer_sigma <- delta_stop_inner_sigma
        delta_stop_inner_theta <- delta_stop_inner_sigma*10
        delta_stop_outer_theta <- delta_stop_inner_theta
      }
      
      # generate data
      {
        data <- gene_data(n[K])
        X <- rbind(X, data$x)
        y <- c(y, data$y)
        N <- N + n[K]
        rm(data)
        # print(paste('K=',K))
      }  
      
      # bandwidth selection
      if(K %in% sub_streams){
        
        h_theta <- G * N^(-1/7); h_sigma <- G * N^(-1/5)
        
        # generate initial estimate and bandwidth
        {
          initial_res <- initial(X,y,h_theta,eval_vec,pd2)
          beta0_theta_old <- initial_res$beta0; beta_theta_old <- initial_res$beta
          beta0_theta <- initial_res$beta0; beta_theta <- initial_res$beta
          rm(initial_res)
          
          initial_res <- initial(X,y,h_sigma,eval_vec,pd1)
          beta0_sigma_old <- initial_res$beta0; beta_sigma_old <- initial_res$beta
          beta0_sigma <- initial_res$beta0; beta_sigma <- initial_res$beta
          rm(initial_res)
        }
        
        
        # backfitting
        {
          res_beta <- backfitting(beta0_theta, beta_theta,U_theta,V_theta,
                                  Q_theta,R_theta, N, N, h_theta, pd2)
          if(res_beta$delta < delta_stop_outer_theta){
            beta0_theta <- res_beta$beta0
            beta_theta <- res_beta$beta
          }
          rm(res_beta)
          
          res_beta <- backfitting(beta0_sigma, beta_sigma,U_sigma,V_sigma,
                                  Q_sigma,R_sigma, N, N, h_sigma, pd1)
          if(res_beta$delta < delta_stop_outer_sigma){
            beta0_sigma <- res_beta$beta0
            beta_sigma <- res_beta$beta
          }
          rm(res_beta)
        }
        
        # compute the constants
        {
          deri <- beta_theta[,(2*d+1):(3*d)]^2
          beta_hat <- matrix(0,nrow(X),d); sigma2 <- rep(0,d)
          for(i in 1:d){
            beta_hat[,i] <- predict(smooth.spline(eval_vec,beta_sigma[,i]),X[,i])$y
          }
          # sigma2 <- mean(exp(rowSums(beta_hat[,(1:d)])+beta0_sigma))
          for(i in 1:d){
            sigma2[i] <- mean(exp(beta_hat[,(1:d)[-i]]))*mean(exp(-beta_sigma[,i]))*exp(-beta0_sigma)
          }
          
          Ch <- (15 * sigma2/colMeans(deri))^(1/5)
          C <- cbind(C, Ch)
          theta_store <- cbind(theta_store, colMeans(deri))
          sigma_store <- cbind(sigma_store, sigma2) 
        }
        
      }
      
    }
    
    save(C, theta_store, sigma_store, 
         file = paste('res/sim1/batch_constant_for_bandwidth_sd',ll,'.Rdata',sep=''))
    
  }
  stopCluster(cl)
  
  C1 <- array(0, dim=c(d,sum(sub_streams<=500),R))
  for(ll in 1:R){
    load(paste('res/sim1/batch_constant_for_bandwidth_sd',ll,'.Rdata',sep=''))
    C1[,,ll] <- C
    file.remove(paste('res/sim1/batch_constant_for_bandwidth_sd',ll,'.Rdata',sep=''))
  }
  save(C1, file='res/sim1/batch_constants_for_bandwidths.Rdata')
}

#######################################################################


############################### bandwidth selection for Simulation2 #########################
rm(list = ls())
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
  pd1 <- 1; pd2 <- 3
  d <- 4
  m <- 10 # No evalpoints
  eval_vec <- seq(0.05, 0.95, length.out = m)
  
  Max_iter <- 50
  N <- 0
  beta_true <- cbind(beta1_fun(eval_vec), beta2_fun(eval_vec), beta3_fun(eval_vec), 
                     beta4_fun(eval_vec), beta5_fun(eval_vec))
  link <- 'identity'
  
  K_band=200
  G <- rep(0.5,d)
  L_theta <- 3; L_sigma <- 3
  time<-c()
  beta0_store<-c()
  beta_store<-c()
  rss <- c()
  band <- c()
}

######### online band select
{
  
  Mcl<-50
  cl<-makeCluster(Mcl)
  registerDoParallel(cl)
  oln<- foreach(sd=sds) %dopar%{
    
    {
      U_theta <- array(0, dim = c(pd2*d+1,m^d,L_theta))
      V_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,m^d,L_theta)) 
      Q_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,m^d,L_theta))
      R_theta <- 0
      U_sigma <- array(0, dim = c(pd1*d+1,m^d,L_sigma))
      V_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L_sigma)) 
      Q_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L_sigma))
      R_sigma <- 0
      C <- c(); theta_store <- c(); sigma_store <- c()
      
    }  
    
    library(MASS)
    set.seed(sd) 
    sub_sds<-ceiling(runif(Kmax)*10^5)
    ll <- which(sds==sd)
    
    for (K in 1:K_band) {
      
      set.seed(sub_sds[K])
      
      # delta_stop
      {
        delta_stop_inner <- 0.01*(K<=10)+0.001*(10<K&K<=30)+1e-04*(K>30)
        delta_stop_outer <- delta_stop_inner
        delta_stop_inner_sigma <- delta_stop_inner*10
        delta_stop_outer_sigma <- delta_stop_inner_sigma
        delta_stop_inner_theta <- delta_stop_inner_sigma*10
        delta_stop_outer_theta <- delta_stop_inner_theta
      }
      
      # generate data
      {
        data <- gene_data(n[K])
        X <- data$x
        y <- data$y
        N <- N + n[K]
        rm(data)
        print(paste('K=',K))
      }  
      
      t0 <- Sys.time()
      # bandwidth selection
      {
        h_theta <- G * N^(-1/7); h_sigma <- G * N^(-1/5)
        eta_theta <- sapply(1:L_theta, function(l){((L_theta-l+1)/L_theta)^(1/7) * h_theta})#dim: d*L
        eta_sigma <- sapply(1:L_sigma, function(l){((L_sigma-l+1)/L_sigma)^(1/5) * h_sigma})
        
        if(K==1){
          
          C<- c()
          
          # generate initial estimate and bandwidth
          initial_res <- initial(X,y,h_theta,eval_vec,pd2)
          beta0_theta_old <- initial_res$beta0; beta_theta_old <- initial_res$beta
          beta0_theta <- initial_res$beta0; beta_theta <- initial_res$beta
          idx_theta <- 1:L_theta
          centrds_theta <- eta_theta
          rm(initial_res)
          
          initial_res <- initial(X,y,h_sigma,eval_vec,pd1)
          beta0_sigma_old <- initial_res$beta0; beta_sigma_old <- initial_res$beta
          beta0_sigma <- initial_res$beta0; beta_sigma <- initial_res$beta
          idx_sigma <- 1:L_sigma
          centrds_sigma <- eta_sigma
          rm(initial_res)
          
        }else{
          
          # update index for pseudo-bandwidth and the centroids
          idx_theta<-sapply(1:L_theta,function(l){
            which.min(abs(eta_theta[1,l] - centrds_theta[1,]))
          })# idx of all d dimensions are the same, it is sufficient to compuute the first dim
          for(i in 1:d){
            centrds_theta[i,] <- (centrds_theta[i,idx_theta] * (N-n[K]) + eta_theta[i,] * n[K]) / N
          }
          
          idx_sigma<-sapply(1:L_sigma,function(l){
            which.min(abs(eta_sigma[1,l] - centrds_sigma[1,]))
          })
          for(i in 1:d){
            centrds_sigma[i,] <- (centrds_sigma[i,idx_sigma] * (N-n[K]) + eta_sigma[i,] * n[K]) / N
          }
          
        }
        
        # backfitting
        {
          res_beta <- backfitting(beta0_theta, beta_theta,U_theta[,,idx_theta[1]],V_theta[,,,idx_theta[1]],
                                  Q_theta[,,,idx_theta[1]],R_theta, n[K], N, h_theta, pd2)
          if(res_beta$delta < delta_stop_outer_theta){
            beta0_theta <- res_beta$beta0
            beta_theta <- res_beta$beta
          }
          rm(res_beta)
          
          res_beta <- backfitting(beta0_sigma, beta_sigma,U_sigma[,,idx_sigma[1]],V_sigma[,,,idx_sigma[1]],
                                  Q_sigma[,,,idx_sigma[1]],R_sigma[,,,,idx_sigma[1]], n[K], N, h_sigma, pd1)
          if(res_beta$delta < delta_stop_outer_sigma){
            beta0_sigma <- res_beta$beta0
            beta_sigma <- res_beta$beta
          }
          rm(res_beta)
        }
        
        # compute the constants
        {
          deri <- beta_theta[,(2*d+1):(3*d)]^2
          beta_hat <- matrix(0,nrow(X),d); sigma2 <- rep(0,d)
          for(i in 1:d){
            beta_hat[,i] <- predict(smooth.spline(eval_vec,beta_sigma[,i]),X[,i])$y
          }
          sigma2 <- mean((y-rowSums(beta_hat[,(1:d)])-beta0_sigma)^2)
          # for(i in 1:d){
          #   sigma2[i] <- mean(exp(beta_hat[,(1:d)[-i]]))*mean(exp(beta_sigma[,i]))*exp(beta0_sigma)
          # }
          
          Ch <- (15 * sigma2/colMeans(deri))^(1/5)
          C <- cbind(C, Ch)
          theta_store <- cbind(theta_store, colMeans(deri))
          sigma_store <- c(sigma_store, sigma2) 
        }
        
        # update statistics
        {
          res_update <- update_stats(U_theta,V_theta,Q_theta,R_theta, beta0_theta,beta_theta,eta_theta,
                                     idx_theta,n[K],N,pd2,L_theta)
          U_theta <- res_update$U_
          V_theta <- res_update$V_
          Q_theta <- res_update$Q_
          R_theta <- res_update$R_
          rm(res_update)
          
          res_update <- update_stats(U_sigma,V_sigma,Q_sigma,R_sigma, beta0_sigma,beta_sigma,eta_sigma,
                                     idx_sigma,n[K],N,pd1,L_sigma)
          U_sigma <- res_update$U_
          V_sigma <- res_update$V_
          Q_sigma <- res_update$Q_
          R_sigma <- res_update$R_
          rm(res_update)
        }
        
      }
      t1 <- Sys.time()
      
      # print(paste0('time=',as.numeric(t1-t0)))
      # print(paste('K=',K))
    }
    
    save(C, theta_store, sigma_store, 
         file = paste('res/sim2/online_constant_for_bandwidth_sd',ll,'.Rdata',sep=''))
    
  }
  stopCluster(cl)
  
  C_o <- array(0, dim=c(d,K_band,R))
  for(ll in 1:R){
    load(paste('res/sim2/online_constant_for_bandwidth_sd',ll,'.Rdata',sep=''))
    C_o[,,ll] <- C
    file.remove(paste('res/sim2/online_constant_for_bandwidth_sd',ll,'.Rdata',sep=''))
  }
  C1 <- C_o
  save(C1, file='res/sim2/online_constants_for_bandwidths.Rdata')
  
}

######### batch band select
{
  Mcl<-50
  cl<-makeCluster(Mcl)
  registerDoParallel(cl)
  bth<- foreach(sd=sds) %dopar%{
    
    {
      U_theta <- array(0, dim = c(pd2*d+1,m^d))
      V_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,m^d)) 
      Q_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,m^d))
      R_theta <- 0
      U_sigma <- array(0, dim = c(pd1*d+1,m^d))
      V_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d)) 
      Q_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d))
      R_sigma <- 0
    }  
    
    C <- c(); theta_store <- c(); sigma_store <- c()
    X <- c(); y <- c()
    library(MASS)
    set.seed(sd) 
    sub_sds<-ceiling(runif(Kmax)*10^5)
    ll <- which(sds==sd)
    
    for (K in 1:K_band) {
      
      set.seed(sub_sds[K])
      
      # delta_stop
      {
        delta_stop_inner <- 0.01*(K<=10)+0.001*(10<K&K<=30)+1e-04*(K>30)
        delta_stop_outer <- delta_stop_inner
      }
      
      # generate data
      {
        data <- gene_data(n[K])
        X <- rbind(X,data$x)
        y <- c(y,data$y)
        N <- N + n[K]
        rm(data)
        print(paste('K=',K))
      }  
      
      t0 <- Sys.time()
      # bandwidth selection
      if(K %in% sub_streams){
        h_theta <- G * N^(-1/7); h_sigma <- G * N^(-1/5)
        eta_theta <- sapply(1:L_theta, function(l){((L_theta-l+1)/L_theta)^(1/7) * h_theta})#dim: d*L
        eta_sigma <- sapply(1:L_sigma, function(l){((L_sigma-l+1)/L_sigma)^(1/5) * h_sigma})
        
        {

          # generate initial estimate and bandwidth
          initial_res <- initial(X,y,h_theta,eval_vec,pd2)
          beta0_theta_old <- initial_res$beta0; beta_theta_old <- initial_res$beta
          beta0_theta <- initial_res$beta0; beta_theta <- initial_res$beta
          rm(initial_res)
          
          initial_res <- initial(X,y,h_sigma,eval_vec,pd1)
          beta0_sigma_old <- initial_res$beta0; beta_sigma_old <- initial_res$beta
          beta0_sigma <- initial_res$beta0; beta_sigma <- initial_res$beta
          rm(initial_res)
          
        }
        
        
        # backfitting
        {
          res_beta <- backfitting(beta0_theta, beta_theta,U_theta,V_theta,
                                  Q_theta,R_theta, N, N, h_theta, pd2)
          if(res_beta$delta < delta_stop_outer){
            beta0_theta <- res_beta$beta0
            beta_theta <- res_beta$beta
          }
          rm(res_beta)
          
          res_beta <- backfitting(beta0_sigma, beta_sigma,U_sigma,V_sigma,
                                  Q_sigma,R_sigma, N, N, h_sigma, pd1)
          if(res_beta$delta < delta_stop_outer){
            beta0_sigma <- res_beta$beta0
            beta_sigma <- res_beta$beta
          }
          rm(res_beta)
        }
        
        # compute the constants
        {
          deri <- beta_theta[,(2*d+1):(3*d)]^2
          beta_hat <- matrix(0,nrow(X),d); sigma2 <- rep(0,d)
          for(i in 1:d){
            beta_hat[,i] <- predict(smooth.spline(eval_vec,beta_sigma[,i]),X[,i])$y
          }
          sigma2 <- mean((y-rowSums(beta_hat[,(1:d)])-beta0_sigma)^2)
          
          Ch <- (15 * sigma2/colMeans(deri))^(1/5)
          C <- cbind(C, Ch)
          theta_store <- cbind(theta_store, colMeans(deri))
          sigma_store <- c(sigma_store, sigma2) 
        }
        
        
      }
      t1 <- Sys.time()
      
      # print(paste0('time=',as.numeric(t1-t0)))
      # print(paste('K=',K))
    }
    
    save(C, theta_store, sigma_store, 
         file = paste('res/sim2/batch_constant_for_bandwidth_sd',ll,'.Rdata',sep=''))
    
  }
  stopCluster(cl)
  
  C_b <- array(0, dim=c(d,sum(sub_streams<=K_band),R))
  for(ll in 1:R){
    load(paste('res/sim2/batch_constant_for_bandwidth_sd',ll,'.Rdata',sep=''))
    C_b[,,ll] <- C
    file.remove(paste('res/sim2/batch_constant_for_bandwidth_sd',ll,'.Rdata',sep=''))
  }
  C1 <- C_b
  save(C1, file='res/sim2/batch_constants_for_bandwidths.Rdata')
  
}

#######################################################################


############################### bandwidth selection for flight #########################
rm(list = ls())
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

######### online band select
{
  {
    U_theta <- array(0, dim = c(pd2*d+1,m^d,L_theta))
    V_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,m^d,L_theta)) 
    Q_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,m^d,L_theta))
    R_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,pd2*d+1,m^d,L_theta))
    U_sigma <- array(0, dim = c(pd1*d+1,m^d,L_sigma))
    V_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L_sigma)) 
    Q_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L_sigma))
    R_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,pd1*d+1,m^d,L_sigma))
    C <- c()
  }  

  library(MASS)
  for(year in 1996:2004){
    
    { 
      load(paste0('datasets/flight/', year, '_gam.Rdata'))
      tab <- table(df$Date)
      dates <- as.numeric(names(tab))[tab>=100]
      if(year>2000){dates <- as.numeric(names(tab))[tab>=4000]}
      Kmax <- length(dates) + start
    }
    
    for (K in (start+1):Kmax) {
      
      tt0 <- Sys.time()
      # generate data
      date <- dates[(K-start)]
      data <- df[(df$Date == date), ]
      X <- cbind(data$CRSDepTime,data$HistDelayRate)
      y <- data$Delayed
      n <- length(y) 
      
      # delta_stop
      {
        delta_stop_inner <- 0.1*(K<=10)+0.01*(10<K&K<=30)+1e-03*(K>30)
        delta_stop_outer <- delta_stop_inner
      }

      # bandwidth selection
      {
        h_theta <- G * N^(-1/7); h_sigma <- G * N^(-1/5)
        eta_theta <- sapply(1:L_theta, function(l){((L_theta-l+1)/L_theta)^(1/7) * h_theta})#dim: d*L
        eta_sigma <- sapply(1:L_sigma, function(l){((L_sigma-l+1)/L_sigma)^(1/5) * h_sigma})
        
        if(K==1){
          
          # generate initial estimate and bandwidth
          initial_res <- initial(X,y,h_theta,eval_vec,pd2)
          beta0_theta_old <- initial_res$beta0; beta_theta_old <- initial_res$beta
          beta0_theta <- initial_res$beta0; beta_theta <- initial_res$beta
          idx_theta <- 1:L_theta
          centrds_theta <- eta_theta
          rm(initial_res)
          
          initial_res <- initial(X,y,h_sigma,eval_vec,pd1)
          beta0_sigma_old <- initial_res$beta0; beta_sigma_old <- initial_res$beta
          beta0_sigma <- initial_res$beta0; beta_sigma <- initial_res$beta
          idx_sigma <- 1:L_sigma
          centrds_sigma <- eta_sigma
          rm(initial_res)
          
        }else{
          
          # update index for pseudo-bandwidth and the centroids
          idx_theta<-sapply(1:L_theta,function(l){
            which.min(abs(eta_theta[1,l] - centrds_theta[1,]))
          })# idx of all d dimensions are the same, it is sufficient to compuute the first dim
          for(i in 1:d){
            centrds_theta[i,] <- (centrds_theta[i,idx_theta] * (N-n[K]) + eta_theta[i,] * n[K]) / N
          }
          
          idx_sigma<-sapply(1:L_sigma,function(l){
            which.min(abs(eta_sigma[1,l] - centrds_sigma[1,]))
          })
          for(i in 1:d){
            centrds_sigma[i,] <- (centrds_sigma[i,idx_sigma] * (N-n[K]) + eta_sigma[i,] * n[K]) / N
          }
          
        }
        
        # backfitting
        {
          res_beta <- backfitting(beta0_theta, beta_theta,U_theta[,,idx_theta[1]],V_theta[,,,idx_theta[1]],
                                  Q_theta[,,,idx_theta[1]],R_theta[,,,,idx_theta[1]], n[K], N, h_theta, pd2)
          if(res_beta$delta < delta_stop_outer_theta){
            beta0_theta <- res_beta$beta0
            beta_theta <- res_beta$beta
          }
          rm(res_beta)
          
          res_beta <- backfitting(beta0_sigma, beta_sigma,U_sigma[,,idx_sigma[1]],V_sigma[,,,idx_sigma[1]],
                                  Q_sigma[,,,idx_sigma[1]],R_sigma[,,,,idx_sigma[1]], n[K], N, h_sigma, pd1)
          if(res_beta$delta < delta_stop_outer_sigma){
            beta0_sigma <- res_beta$beta0
            beta_sigma <- res_beta$beta
          }
          rm(res_beta)
        }
        
        # compute the constants
        {
          deri <- beta_theta[,(2*d+1):(3*d)]^2
          beta_hat <- matrix(0,nrow(X),d); sigma2 <- rep(0,d)
          for(i in 1:d){
            beta_hat[,i] <- predict(smooth.spline(eval_vec,beta_sigma[,i]),X[,i])$y
          }
          # sigma2 <- mean(exp(rowSums(beta_hat[,(1:d)])+beta0_sigma))
          for(i in 1:d){
            sigma2[i] <- mean(exp(beta_hat[,(1:d)[-i]]))*mean(exp(-beta_sigma[,i]))*exp(-beta0_sigma)
          }
          
          Ch <- (15 * sigma2/colMeans(deri))^(1/5)
          C <- cbind(C, Ch)
        }
        
        # update statistics
        {
          res_update <- update_stats(U_theta,V_theta,Q_theta,R_theta, beta0_theta,beta_theta,eta_theta,
                                     idx_theta,n[K],N,pd2,L_theta)
          U_theta <- res_update$U_
          V_theta <- res_update$V_
          Q_theta <- res_update$Q_
          R_theta <- res_update$R_
          rm(res_update)
          
          res_update <- update_stats(U_sigma,V_sigma,Q_sigma,R_sigma, beta0_sigma,beta_sigma,eta_sigma,
                                     idx_sigma,n[K],N,pd1,L_sigma)
          U_sigma <- res_update$U_
          V_sigma <- res_update$V_
          Q_sigma <- res_update$Q_
          R_sigma <- res_update$R_
          rm(res_update)
        }
        
      }
      
      print(paste('K=',K))
    }
    
    start <- Kmax
  }
  C1 <- C
  save(C1,  file = 'res/flight/online_constants_for_bandwidths.Rdata')
}

######### batch band select
{
  {
    U_theta <- array(0, dim = c(pd2*d+1,m^d))
    V_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,m^d)) 
    Q_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,m^d))
    R_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,pd2*d+1,m^d))
    U_sigma <- array(0, dim = c(pd1*d+1,m^d))
    V_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d)) 
    Q_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d))
    R_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,pd1*d+1,m^d))
    X <- c(); y <- c(); N <- 0; C<- c()
  }  
  
  for(year in 1996:2004){
    
    { 
      load(paste0('datasets/flight/', year, '_gam.Rdata'))
      tab <- table(df$Date)
      dates <- as.numeric(names(tab))[tab>=100]
      if(year>2000){dates <- as.numeric(names(tab))[tab>=4000]}
      Kmax <- length(dates) + start
    }
    
    for (K in (start+1):Kmax) {
      
      tt0 <- Sys.time()
      # generate data
      date <- dates[(K-start)]
      data <- df[(df$Date == date), ]
      X <- cbind(data$CRSDepTime,data$HistDelayRate)
      y <- data$Delayed
      n <- length(y)
      
      # delta_stop
      {
        delta_stop_inner <- 0.01*(K<=10)+0.001*(10<K&K<=30)+1e-04*(K>30)
        delta_stop_outer <- delta_stop_inner
      }

      # bandwidth selection
      if(K %in% sub_streams){
        
        h_theta <- G * N^(-1/7); h_sigma <- G * N^(-1/5)
        
        # generate initial estimate and bandwidth
        {
          initial_res <- initial(X,y,h_theta,eval_vec,pd2)
          beta0_theta_old <- initial_res$beta0; beta_theta_old <- initial_res$beta
          beta0_theta <- initial_res$beta0; beta_theta <- initial_res$beta
          rm(initial_res)
          
          initial_res <- initial(X,y,h_sigma,eval_vec,pd1)
          beta0_sigma_old <- initial_res$beta0; beta_sigma_old <- initial_res$beta
          beta0_sigma <- initial_res$beta0; beta_sigma <- initial_res$beta
          rm(initial_res)
        }

        # backfitting
        {
          res_beta <- backfitting(beta0_theta, beta_theta,U_theta,V_theta,
                                  Q_theta,R_theta, N, N, h_theta, pd2)
          if(res_beta$delta < delta_stop_outer_theta){
            beta0_theta <- res_beta$beta0
            beta_theta <- res_beta$beta
          }
          rm(res_beta)
          
          res_beta <- backfitting(beta0_sigma, beta_sigma,U_sigma,V_sigma,
                                  Q_sigma,R_sigma, N, N, h_sigma, pd1)
          if(res_beta$delta < delta_stop_outer_sigma){
            beta0_sigma <- res_beta$beta0
            beta_sigma <- res_beta$beta
          }
          rm(res_beta)
        }
        
        # compute the constants
        {
          deri <- beta_theta[,(2*d+1):(3*d)]^2
          beta_hat <- matrix(0,nrow(X),d); sigma2 <- rep(0,d)
          for(i in 1:d){
            beta_hat[,i] <- predict(smooth.spline(eval_vec,beta_sigma[,i]),X[,i])$y
          }
          # sigma2 <- mean(exp(rowSums(beta_hat[,(1:d)])+beta0_sigma))
          for(i in 1:d){
            sigma2[i] <- mean(exp(beta_hat[,(1:d)[-i]]))*mean(exp(-beta_sigma[,i]))*exp(-beta0_sigma)
          }
          
          Ch <- (15 * sigma2/colMeans(deri))^(1/5)
          C <- cbind(C, Ch)
        }
        
      }
      
    }
    
    start <- Kmax
  }
  C1 <- C
  save(C1,  file = 'res/flight/batch_constants_for_bandwidths.Rdata')
}


#######################################################################


############################### bandwidth selection for credit #########################
source('FNS/FNS_SmoBack_credit.R')

# read in data
{
  if(!file.exists('datasets/credit/P2P_Macro_Data_processed.Rdata')){
    
    library(haven)
    
    data <- read_dta('P2P_Macro_Data.dta')
    data <- data[,c('badloan','grade_','loan_amnt','int_rate','total_pymnt','loan_status')]
    data <- data[data$loan_status!="Current",]
    data <- data[,c('badloan','grade_','loan_amnt','int_rate','total_pymnt')]
    
    gc()
    
    data$int_rate <- log(data$int_rate)
    data$total_pymnt <- log(data$total_pymnt)
    data <- data[which(data$total_pymnt>6.6&data$total_pymnt<10.8),]
    
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
  Xtrain <- as.matrix(data[idx,c('loan_amnt','int_rate','dti','total_pymnt')])
  ytest <- data$badloan[(1:Nfull)[-idx]]
  Xtest <- as.matrix(data[(1:Nfull)[-idx],c('loan_amnt','int_rate','dti','total_pymnt')])
  rm(idx)
}

#### parameter
{
  pd1 <- 1; pd2 <- 3
  d <- 4
  
  m <- 16 # No evalpoints
  eval_vec <- seq(0.05, 0.95, length.out = m)
  
  Max_iter <- 50
  K_band <- 200
  L_theta <- 5; L_sigma <- 5
  link <- 'logit'
  
  nn <- rep(1000,867)
  nn[1] <- 3000
  Kmax <- length(nn)
  nn[Kmax] <- Ntrain-sum(nn[1:(Kmax-1)])
  
}

######### online band select
{

  {
    U_theta <- array(0, dim = c(pd2*d+1,m^d,L_theta))
    V_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,m^d,L_theta)) 
    Q_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,m^d,L_theta))
    R_theta <- 0
    U_sigma <- array(0, dim = c(pd1*d+1,m^d,L_sigma))
    V_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L_sigma)) 
    Q_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L_sigma))
    R_sigma <- 0
    C <- c()
  }  
  
  for (K in 1:K_band) {
    
    X <- as.matrix(Xtrain[(NN+1):(NN+nn[K]),])
    y <- ytrain[(NN+1):(NN+nn[K])]
    
    # delta
    delta_inner <- 0.01*(K<=10)+0.001*(10<K&K<=50)+1e-04*(K>50)
    # delta_outer <- 0.2*(K<=10)+0.1*(10<K&K<=30)+0.01*(30<K&K<=50)+1e-03*(K>50)
    delta_outer <- 0.1*(K<=20)+0.01*(20<K&K<=50)+1e-03*(K>50)
    
    # bandwidth selection
    h_theta <- G * N^(-1/7); h_sigma <- G * N^(-1/5)
    eta_theta <- sapply(1:L_theta, function(l){((L_theta-l+1)/L_theta)^(1/7) * h_theta})#dim: d*L
    eta_sigma <- sapply(1:L_sigma, function(l){((L_sigma-l+1)/L_sigma)^(1/5) * h_sigma})
    
    if(K==1){
      
      C<- c()
      
      # generate initial estimate and bandwidth
      initial_res <- initial(X,y,h_theta,eval_vec,pd2)
      beta0_theta_old <- initial_res$beta0; beta_theta_old <- initial_res$beta
      beta0_theta <- initial_res$beta0; beta_theta <- initial_res$beta
      idx_theta <- 1:L_theta
      centrds_theta <- eta_theta
      rm(initial_res)
      
      initial_res <- initial(X,y,h_sigma,eval_vec,pd1)
      beta0_sigma_old <- initial_res$beta0; beta_sigma_old <- initial_res$beta
      beta0_sigma <- initial_res$beta0; beta_sigma <- initial_res$beta
      idx_sigma <- 1:L_sigma
      centrds_sigma <- eta_sigma
      rm(initial_res)
      
    }else{
      
      # update index for pseudo-bandwidth and the centroids
      idx_theta<-sapply(1:L_theta,function(l){
        which.min(abs(eta_theta[1,l] - centrds_theta[1,]))
      })# idx of all d dimensions are the same, it is sufficient to compuute the first dim
      for(i in 1:d){
        centrds_theta[i,] <- (centrds_theta[i,idx_theta] * (N-n[K]) + eta_theta[i,] * n[K]) / N
      }
      
      idx_sigma<-sapply(1:L_sigma,function(l){
        which.min(abs(eta_sigma[1,l] - centrds_sigma[1,]))
      })
      for(i in 1:d){
        centrds_sigma[i,] <- (centrds_sigma[i,idx_sigma] * (N-n[K]) + eta_sigma[i,] * n[K]) / N
      }
      
    }
    
    # backfitting
    {
      res_beta <- backfitting(beta0_theta, beta_theta,U_theta[,,idx_theta[1]],V_theta[,,,idx_theta[1]],
                              Q_theta[,,,idx_theta[1]],R_theta, n[K], N, h_theta, pd2)
      if(res_beta$delta < delta_stop_outer_theta){
        beta0_theta <- res_beta$beta0
        beta_theta <- res_beta$beta
      }
      rm(res_beta)
      
      res_beta <- backfitting(beta0_sigma, beta_sigma,U_sigma[,,idx_sigma[1]],V_sigma[,,,idx_sigma[1]],
                              Q_sigma[,,,idx_sigma[1]],R_sigma[,,,,idx_sigma[1]], n[K], N, h_sigma, pd1)
      if(res_beta$delta < delta_stop_outer_sigma){
        beta0_sigma <- res_beta$beta0
        beta_sigma <- res_beta$beta
      }
      rm(res_beta)
    }
    
    # compute the constants
    {
      deri <- beta_theta[,(2*d+1):(3*d)]^2
      beta_hat <- matrix(0,nrow(X),d); sigma2 <- rep(0,d)
      for(i in 1:d){
        beta_hat[,i] <- predict(smooth.spline(eval_vec,beta_sigma[,i]),X[,i])$y
      }
      sigma2 <- mean((y-rowSums(beta_hat[,(1:d)])-beta0_sigma)^2)
      # for(i in 1:d){
      #   sigma2[i] <- mean(exp(beta_hat[,(1:d)[-i]]))*mean(exp(beta_sigma[,i]))*exp(beta0_sigma)
      # }
      
      Ch <- (15 * sigma2/colMeans(deri))^(1/5)
      C <- cbind(C, Ch)
    }
    
    # update statistics
    {
      res_update <- update_stats(U_theta,V_theta,Q_theta,R_theta, beta0_theta,beta_theta,eta_theta,
                                 idx_theta,n[K],N,pd2,L_theta)
      U_theta <- res_update$U_
      V_theta <- res_update$V_
      Q_theta <- res_update$Q_
      R_theta <- res_update$R_
      rm(res_update)
      
      res_update <- update_stats(U_sigma,V_sigma,Q_sigma,R_sigma, beta0_sigma,beta_sigma,eta_sigma,
                                 idx_sigma,n[K],N,pd1,L_sigma)
      U_sigma <- res_update$U_
      V_sigma <- res_update$V_
      Q_sigma <- res_update$Q_
      R_sigma <- res_update$R_
      rm(res_update)
    }
    
  }
  
  C1 <- C
  save(C1, file='res/credit/online_constants_for_bandwidths.Rdata')
  
}

######### batch band select
{
  
  {
    U_theta <- array(0, dim = c(pd2*d+1,m^d))
    V_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,m^d)) 
    Q_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,m^d))
    R_theta <- 0
    U_sigma <- array(0, dim = c(pd1*d+1,m^d))
    V_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d)) 
    Q_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d))
    R_sigma <- 0
    C <- c(); X <- c(); y <- c()
  }  
  
  for (K in 1:K_band) {
    
    X <- rbind(X, as.matrix(Xtrain[(NN+1):(NN+nn[K]),]))
    y <- c(y, ytrain[(NN+1):(NN+nn[K])])
    
    delta_inner <- 0.01*(K<=10)+0.001*(10<K&K<=50)+1e-04*(K>50)
    # delta_outer <- 0.2*(K<=10)+0.1*(10<K&K<=30)+0.01*(30<K&K<=50)+1e-03*(K>50)
    delta_outer <- 0.1*(K<=20)+0.01*(20<K&K<=50)+1e-03*(K>50)
    
    # bandwidth selection
    if(K %in% sub_streams){
      
      h_theta <- G * N^(-1/7); h_sigma <- G * N^(-1/5)
      eta_theta <- sapply(1:L_theta, function(l){((L_theta-l+1)/L_theta)^(1/7) * h_theta})#dim: d*L
      eta_sigma <- sapply(1:L_sigma, function(l){((L_sigma-l+1)/L_sigma)^(1/5) * h_sigma})
      
      # generate initial estimate
      {
        initial_res <- initial(X,y,h_theta,eval_vec,pd2)
        beta0_theta_old <- initial_res$beta0; beta_theta_old <- initial_res$beta
        beta0_theta <- initial_res$beta0; beta_theta <- initial_res$beta
        rm(initial_res)
        
        initial_res <- initial(X,y,h_sigma,eval_vec,pd1)
        beta0_sigma_old <- initial_res$beta0; beta_sigma_old <- initial_res$beta
        beta0_sigma <- initial_res$beta0; beta_sigma <- initial_res$beta
        rm(initial_res)
      }
      
      # backfitting
      {
        res_beta <- backfitting(beta0_theta, beta_theta,U_theta,V_theta,
                                Q_theta,R_theta, N, N, h_theta, pd2)
        if(res_beta$delta < delta_stop_outer){
          beta0_theta <- res_beta$beta0
          beta_theta <- res_beta$beta
        }
        rm(res_beta)
        
        res_beta <- backfitting(beta0_sigma, beta_sigma,U_sigma,V_sigma,
                                Q_sigma,R_sigma, N, N, h_sigma, pd1)
        if(res_beta$delta < delta_stop_outer){
          beta0_sigma <- res_beta$beta0
          beta_sigma <- res_beta$beta
        }
        rm(res_beta)
      }
      
      # compute the constants
      {
        deri <- beta_theta[,(2*d+1):(3*d)]^2
        beta_hat <- matrix(0,nrow(X),d); sigma2 <- rep(0,d)
        for(i in 1:d){
          beta_hat[,i] <- predict(smooth.spline(eval_vec,beta_sigma[,i]),X[,i])$y
        }
        sigma2 <- mean((y-rowSums(beta_hat[,(1:d)])-beta0_sigma)^2)
        
        Ch <- (15 * sigma2/colMeans(deri))^(1/5)
        C <- cbind(C, Ch)
      }
      
    }
    
  }
  
  
  C1 <- C
  save(C1, file='res/credit/batch_constants_for_bandwidths.Rdata')
  
}
#######################################################################
