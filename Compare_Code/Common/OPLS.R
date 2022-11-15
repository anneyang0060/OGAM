#
# The following code is adapted/copied from the supplemental material of the paper:
#     Yang and Yao (2022+): Online Estimation for Functional Data. 
#                           https://doi.org/10.1080/01621459.2021.2002158
#' helper functions
{  
  Epan <- function(z)
  { 
    return( 3/4 * (1-z^2) * (abs(z)<1) ) 
  } 
  
  online_LCub <- function(x, y, eval, h, L, res_list, N, n, d,K)
  {
    eta <- sapply(1:L, function(l){ ((L-l+1) / L) ^ (1/(6+d)) * h}) 
    if(K>1){
      idx <- sapply(1:L,function(l){which.min(abs(eta[l] - res_list$centroids))})}
    else{idx <- 1:L}
    res_list$centroids <- (res_list$centroids[idx] * (N-n) + eta * n) / N
    { 
      EV <- length(eval)
      for(l in 1:L){
        Pnew <- array(0, dim = c(4,4,EV)); qnew <- matrix(0,4,EV)
        for(i in 1:EV){
          side <- cbind(1, x - eval[i], (x - eval[i])^2, (x - eval[i])^3)
          K_vec <- Epan((x - eval[i])/eta[l])/eta[l]
          for(nr in 1:4){
            for(nc in 1:4){
              Pnew[nr,nc,i] <- sum(K_vec*side[,nr]*side[,nc])/n }
            qnew[nr,i] <-  sum(K_vec*side[,nr]*y)/n
          }
        }
        res_list$P[,,,l] <- (res_list$P[,,,idx[l]] * (N-n) + Pnew * n) / N
        res_list$q[,,l] <- (res_list$q[,,idx[l]] * (N-n) + qnew * n) / N
      }
    }
    return(res_list)
  }
  
  online_LL <- function(x, y, eval, h, L, res_list, N, n, d,K)
  {
    eta <- sapply(1:L, function(l){ ((L-l+1) / L) ^ (1/(d+4)) * h})
    # print(h)
    # print(eta)
    if(K>1){
      idx <- sapply(1:L,function(l){which.min(abs(eta[l] - res_list$centroids))})}
    else{
      idx <- 1:L
    }
    # print(idx)
    # print(res_list$centroids)
    res_list$centroids <- (res_list$centroids[idx] * (N-n) + eta * n) / N
    {
      EV <- length(eval) 
      for(l in 1:L){
        Pnew <- array(0, dim = c(2,2,EV)); qnew <- matrix(0,2,EV)
        for(i in 1:EV){
          side <- cbind(1, x - eval[i])
          K_vec <- Epan((x - eval[i])/eta[l])/eta[l]
          Pnew[,,i] <- matrix(c(
            sum(K_vec*side[,1]^2), sum(K_vec*side[,1]*side[,2]),
            sum(K_vec*side[,2]*side[,1]),  sum(K_vec*side[,2]^2)
          ),2,2) / n
          qnew[,i] <- matrix(c(
            sum(K_vec*side[,1]*y), sum(K_vec*side[,2]*y)
          ),2,1) / n
        }
        res_list$P[,,,l] <- (res_list$P[,,,idx[l]] * (N-n) + Pnew * n) / N
        res_list$q[,,l] <- (res_list$q[,,idx[l]] * (N-n) + qnew * n) / N
      }
    }
    # print(res_list$centroids)
    return(res_list)
  }

  batch_LCub <- function(x, y, eval, h, N, d)
  {
    if(d==1){ 
      EV <- length(eval)
      P <- array(0, dim = c(4,4,EV)); q <- matrix(0,4,EV)
      for(i in 1:EV){
        side <- cbind(1, x - eval[i], (x - eval[i])^2, (x - eval[i])^3)
        K_vec <- Epan((x - eval[i])/h)/h
        for(nr in 1:4){
          for(nc in 1:4){
            P[nr,nc,i] <- sum(K_vec*side[,nr]*side[,nc])/N
          }
          q[nr,i] <-  sum(K_vec*side[,nr]*y)/N
        }
      }
      
      sec_deri <- sapply(1:EV, function(i){
        2*(solve(P[,,i]+diag(1e-12,4)) %*% matrix(q[,i],4,1))[3]
      })
      index1 <- 1  
      index2 <- 50
      theta <- sum(P[1,1,index1:index2]
                    * sec_deri[index1:index2]^2)/sum(P[1,1,index1:index2])}
    return(theta)
  }
  
  batch_LL <- function(x, y, eval, h, N, d)
  {
    if(d==1){
      EV <- length(eval)
      P <- array(0, dim = c(2,2,EV)); q <- matrix(0,2,EV)
      for(i in 1:EV){
        side <- cbind(1, x - eval[i])
        K_vec <- Epan((x - eval[i])/h)/h
        P[,,i] <- matrix(c(
          sum(K_vec*side[,1]^2), sum(K_vec*side[,1]*side[,2]),
          sum(K_vec*side[,2]*side[,1]),  sum(K_vec*side[,2]^2)
        ),2,2) / N
        q[,i] <-  matrix(c(
          sum(K_vec*side[,1]*y), sum(K_vec*side[,2]*y)
        ),2,1) / N
      }
      den <- sapply(1:EV, function(i){P[1,1,i]})
      est <- sapply(1:EV, function(i){
        (solve(P[,,i]+diag(1e-12,2)) %*% matrix(q[,i],2,1))[1]
      })}
    res<- list(den,est)
    names(res) <- c('den','est')
    return(res)
  }
  
  RMISEfunc<-function(t,f1,f2)
  {
    x=t
    y=(f1-f2)^2
    return(mean(y))
  }
}

#' function to estimate regression function with one-pass local smoothing method OPLS method) 
#' @param data_batch the data in the batch form
#' @param domain the t domain
#' @param m true regression function
#' @return a list containing RMISE, time of maintaining and querying
batchrunlocal=function(data_batch,domain,m){
  # common parameters
  {
    a <- domain[1]+0.05; b <- domain[2]-0.05
    EV1 <- 50; 
    eval_mu <- seq(a,b,length.out = EV1)
    G <- 0.7
    K <- length(data_batch)
    Kmax <- K
    K1 <- 300
  }
  # initialize
  {
    mus <- c(); h1 <- c(); sigma <- c(); maintaintime=c(); querytime=c()
    L1 <- 20
    N <- 0
    h0 <- 0.7
    theta_mu <- 0; sigma_mu_new <- 0
    
    res_theta_mu <- list()
    res_theta_mu$centroids <- rep(0, L1)
    res_theta_mu$P <- array(0, dim = c(4,4,EV1,L1))
    res_theta_mu$q <- array(0, dim = c(4,EV1,L1))
    
    res_sigma_mu1 <- list()
    res_sigma_mu1$centroids <- rep(0, L1)
    res_sigma_mu1$P <- array(0, dim = c(2,2,EV1,L1))
    res_sigma_mu1$q <- array(0, dim = c(2,EV1,L1))
    res_mu <- res_sigma_mu1
  }
  # update
  for(K in 1:Kmax)
  {
    print(paste('K=',K))
    t0 <- Sys.time()
    # generate data
    { 
      x <- data_batch[[K]][,1]
      y <-data_batch[[K]][,2]
      NK=length(x)
      N <- N + NK
    }
    # sigma_mu 
    if(K<=K1)
    {
      h_theta_mu <- G * N^(-1/7)
      res_theta_mu <- online_LCub(x=x, y=y, eval=eval_mu, h=h_theta_mu, L=L1,res_list=res_theta_mu, N=N, n=NK, d=1,K=K)
      mu_sec_deri <- sapply(1:EV1, function(i){2*(solve(res_theta_mu$P[,,i,1]+diag(1e-12,4)) %*% matrix(res_theta_mu$q[,i,1],4,1))[3]})
      theta_mu <-  sum(res_theta_mu$P[1,1,1:50 ,1] * mu_sec_deri[1:50]^2)/sum(res_theta_mu$P[1,1,1:50 ,1])
      
      h_sigma_mu <- G * N^(-1/5)
      res_sigma_mu1 <- online_LL(x, y,eval_mu, h_sigma_mu, L1, res_sigma_mu1, N, NK, 1,K)
      mu <- sapply(1:EV1, function(i){(solve(res_sigma_mu1$P[,,i,1]+diag(1e-12,2)) %*% matrix(res_sigma_mu1$q[,i,1],2,1))[1]})
      sp <- smooth.spline(eval_mu,mu)
      mu_est<-predict(sp,x)$y
      sigma_mu_new <- (N-NK)/N * sigma_mu_new + NK/N * mean((y-mu_est)^2)
      # sigma_mu_new <- mean((y-mu_est)^2)
    }
    
    # estimate
    {
      h_mu <- min((15 * sigma_mu_new / theta_mu)^(1/5) * N^(-1/5),h0)
      h0 <- h_mu
      res_mu <- online_LL(x, y,eval_mu, h_mu, L1, res_mu, N, NK, 1,K)
      t1=Sys.time()
      mu <- sapply(1:EV1, function(i){(solve(res_mu$P[,,i,1]+diag(1e-12,2)) %*% matrix(res_mu$q[,i,1],2,1))[1]})
      t2=Sys.time()
    }
    # results
    {
      mus <- cbind(mus, mu)
      h1 <- c(h1,h_mu)
      sigma <- c(sigma,sigma_mu_new)
      maintaintime <- c(maintaintime,difftime(t1,t0,units = 'secs'))
      querytime <- c(querytime,difftime(t2,t1,units = 'secs'))
    }
  }
  RMISE=rep(0,K)
  mut=m(eval_mu)
  for(i in 1:K){RMISE[i]=RMISEfunc(eval_mu,mus[,i],mut)}
  R=list(RMISE,maintaintime,querytime,h1)
  return(R)
}

fullrun=function(data_batch,domain,m){
  # common parameters
  {
    a <- domain[1]+0.05; b <- domain[2]-0.05
    EV1 <- 50; 
    eval_mu <- seq(a,b,length.out = EV1)
    G <- 0.7
    K <- length(data_batch)
    Kmax <- K
    K1 <- 300
  }
  # initialize
  {
    mus <- c(); h1 <- c(); sigma <- c(); maintaintime=c(); querytime=c()
    N <- 0; x <- c(); y <- c()
  }
  # update
  for(K in 1:Kmax)
  {
    print(paste('K=',K))
    t0 <- Sys.time()
    # generate data
    { 
      x <- c(x, data_batch[[K]][,1])
      y <- c(y, data_batch[[K]][,2])
      N <- length(y)
    }
    # sigma_mu 
    if(K<=K1)
    {
      {
        h_theta_mu <- G * N^(-1/7)
        theta_mu <- batch_LCub(x=x, y=y, eval=eval_mu, h=h_theta_mu, N=N, d=1)
        h_sigma_mu <- G * N^(-1/5)
        res_sigma_mu1 <- batch_LL(x=x, y=y, eval=eval_mu, h=h_sigma_mu, N=N, d=1)
        mu <- res_sigma_mu1$est
        sp <- smooth.spline(eval_mu,mu)
        mu_est<-predict(sp,x)$y
        sigma_mu_new <- mean((y-mu_est)^2)
      }
    }
    # estimate
    {
      h_mu <- (15 * sigma_mu_new / theta_mu)^(1/5) * N^(-1/5)  
      res_mu <- batch_LL(x=x, y=y, eval=eval_mu, h=h_mu, N=N, d=1)
      t1=Sys.time()
      mu <- res_mu$est
      t2=Sys.time()
    }
    # results
    {
      mus <- cbind(mus, mu)
      h1 <- c(h1,h_mu)
      maintaintime <- c(maintaintime,difftime(t1,t0,units = 'secs'))
      querytime <- c(querytime,difftime(t2,t1,units = 'secs'))
    }
  }
  RMISE=rep(0,K)
  mut=m(eval_mu)
  for(i in 1:K){RMISE[i]=RMISEfunc(eval_mu,mus[,i],mut)}
  R=list(RMISE,maintaintime,querytime,h1)
  return(R)
}



#' function to estimate regression function with shrinked and fixed constant for bandwidth selection
#' @param data_batch the data in the batch form
#' @param domain the t domain
#' @param m true regression function
#' @return a list containing RMISE, time of maintaining and querying
batchrunshrink=function(data_batch,domain,m){
  # common parameters
  {
    a <- domain[1]+0.05; b <- domain[2]-0.05
    EV1 <- 50; 
    eval_mu <- seq(a,b,length.out = EV1)
    G <- 0.7
    K <- length(data_batch)
    Kmax <- K
    K1 <- 300
  }
  # initialize
  {
    mus <- c(); h1 <- c(); sigma <- c(); maintaintime=c(); querytime=c()
    L1 <- 1
    N <- 0
    h0 <- 0.7
    theta_mu <- 0; sigma_mu_new <- 0
    
    res_theta_mu <- list()
    res_theta_mu$centroids <- rep(0, L1)
    res_theta_mu$P <- array(0, dim = c(4,4,EV1,L1))
    res_theta_mu$q <- array(0, dim = c(4,EV1,L1))
    
    res_sigma_mu1 <- list()
    res_sigma_mu1$centroids <- rep(0, L1)
    res_sigma_mu1$P <- array(0, dim = c(2,2,EV1,L1))
    res_sigma_mu1$q <- array(0, dim = c(2,EV1,L1))
    res_sigma_mu2 <- res_sigma_mu1
    res_mu <- res_sigma_mu1
  }
  # update
  for(K in 1:Kmax)
  {
    print(paste('K=',K))
    t0 <- Sys.time()
    # generate data
    { 
      x <- data_batch[[K]][,1]
      y <-data_batch[[K]][,2]
      NK=length(x)
      N <- N + NK
    }
    # sigma_mu 
    if(K<=K1)
    {
      {
        h_theta_mu <- G * N^(-1/7)
        res_theta_mu <- online_LCub(x=x, y=y, eval=eval_mu, h=h_theta_mu, L=L1,res_list=res_theta_mu, N=N, n=NK, d=1,K=K)
        mu_sec_deri <- sapply(1:EV1, function(i){2*(solve(res_theta_mu$P[,,i,1]+diag(1e-12,4)) %*% matrix(res_theta_mu$q[,i,1],4,1))[3]})
        theta_mu <-  sum(res_theta_mu$P[1,1,1:50 ,1] * mu_sec_deri[1:50]^2)/sum(res_theta_mu$P[1,1,1:50 ,1])
        h_sigma_mu <- G * N^(-1/5)
        res_sigma_mu1 <- online_LL(x, y,eval_mu, h_sigma_mu, L1, res_sigma_mu1, N, NK, 1,K)
        mu <- sapply(1:EV1, function(i){(solve(res_sigma_mu1$P[,,i,1]+diag(1e-12,2)) %*% matrix(res_sigma_mu1$q[,i,1],2,1))[1]})
        sp <- smooth.spline(eval_mu,mu)
        mu_est<-predict(sp,x)$y
        sigma_mu_new <- (N-NK)/N * sigma_mu_new + NK/N * mean((y-mu_est)^2)
        # sigma_mu_new <- mean((y-mu_est)^2)
      }
    }
    # estimate
    {
      h_mu <- sqrt(0.3)*min((15 * sigma_mu_new / theta_mu)^(1/5) * N^(-1/5),h0)
      h0 <- h_mu/sqrt(0.3)
      res_mu <- online_LL(x, y,eval_mu, h_mu, L1, res_mu, N, NK, 1,K)
      t1=Sys.time()
      mu <- sapply(1:EV1, function(i){(solve(res_mu$P[,,i,1]+diag(1e-12,2)) %*% matrix(res_mu$q[,i,1],2,1))[1]})
      t2=Sys.time()
    }
    # results
    {
      mus <- cbind(mus, mu)
      h1 <- c(h1,h_mu)
      sigma <- c(sigma,sigma_mu_new)
      maintaintime <- c(maintaintime,difftime(t1,t0,units = 'secs'))
      querytime <- c(querytime,difftime(t2,t1,units = 'secs'))
    }
  }
  RMISE=rep(0,K)
  mut=m(eval_mu)
  for(i in 1:K){RMISE[i]=RMISEfunc(eval_mu,mus[,i],mut)}
  R=list(RMISE,maintaintime,querytime)
  return(R)
}
