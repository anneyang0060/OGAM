### initial: use glm to get an intial estimate, 
###          return beta (a (pd+1)*d-dim matrx) whose columns is composed of:
###          beta_1,...,beta_d,beta'_1,...,beta'_d, beta"_1,...,beta"_d,...
initial <- function(X, y, h, eval_vec,pd){
  
  d <- ncol(X); M <- length(eval_vec)
  beta <- matrix(0,length(eval_vec), (pd+1)*d)
  X <- cbind(X,X^2,X^3)
  if(link=='logit'){model_glm<-glm(y~X,binomial(link = 'logit'))}
  if(link=='log'){model_glm<-glm(y~X,poisson(link = 'log'))}
  if(link=='identity'){model_glm<-glm(y~X,gaussian(link = "identity"))}
  c_para <- unname(model_glm$coefficients)
  
  beta0 <- c_para[1]
  
  for(i in 1:d){
    beta[,i] <- (cbind(eval_vec,eval_vec^2,eval_vec^3)
                 %*%matrix(c_para[c(i+1,i+1+d,i+1+2*d)],3,1))
    beta[,i+d] <- (h[i]*cbind(1,2*eval_vec,3*eval_vec^2)
                   %*%matrix(c_para[c(i+1,i+1+d,i+1+2*d)],3,1))
    if(pd > 1){
      beta[,i+2*d] <- h[i]^2*cbind(2,6*eval_vec)%*%matrix(c_para[c(i+1+d,i+1+2*d)],2,1)/2
      beta[,i+3*d] <- h[i]^3*6*c_para[i+1+2*d]/6
    }
  }
  
  v <- beta0
  for(i in 1:d){
    v <- v+ rep(rep(1:M,each=M^(i-1)),M^(d-i))
  }
  v <- exp(v)
  for(i in 1:d){
    c_tmpt <- mean(rep(rep(1:M,each=M^(i-1)),M^(d-i)) * v) / mean(v)
    beta0 <- beta0 + c_tmpt
    beta[,i] <- beta[,i] - c_tmpt
  }
  
  beta0 <- beta0 + sum(colMeans(beta[,1:d]))
  beta[,1:d] <- beta[,1:d] - matrix(rep(colMeans(beta[,1:d]),each=M),ncol=d)
  
  res <- list(beta0,beta)
  names(res) <- c('beta0', 'beta')
  return(res)
}

### functions to compute derivatives of Q functions w.r.t. beta
m_fun <- function(beta){
  if(link=='logit'){y <- 1/(1+exp(-beta))}
  if(link=='log'){y <- exp(beta)}
  if(link=='identity'){y <- beta}
  return(y)
}
q1 <- function(y1, beta){
  y <- y1 - m_fun(beta)
  return(y)
}
q2 <- function(beta){
  if(link=='logit'){y <- -m_fun(beta) * (1-m_fun(beta))}
  if(link=='log'){y <- -m_fun(beta)}
  if(link=='identity'){y <- -1}
  return(y)
}
q3 <- function(beta){
  if(link=='logit'){y <- -m_fun(beta) * (1-m_fun(beta)) * (1-2*m_fun(beta))}
  if(link=='log'){y <- -m_fun(beta)}
  return(y)
}

Epan<-function(z){3/4*(1-z^2)*(abs(z)<1)} 

## compute_PQR: Given beta(x), compute statistics P,Q,R
compute_PQR<-function(beta0, beta, n, h, pd){
  
  P <- array(0, dim = c(pd*d+1, m^d))
  Q <- array(0, dim = c(pd*d+1, pd*d+1, m^d))
  if(link!='identity'){R <- array(0, dim = c(pd*d+1, pd*d+1, pd*d+1, m^d))}
  
  # Cartesian product of evaluation points
  idx_eval <- c()
  for(i in 1:d){
    idx_eval <- cbind(idx_eval, rep(rep(1:m,each=m^(i-1)),m^(d-i)))
  }
  # print(beta0)
  # print(beta)
  
  for(i in 1:nrow(idx_eval))
    {

    idx1 <- idx_eval[i,]
    beta1 <- beta0; x <- c()
    for(si in 1:d){
      beta1 <- beta1 + beta[idx1[si],si]
      x <- c(x,eval_vec[idx1[si]])
    }
    for(si1 in 1:pd)
      for(si2 in 1:d){
        beta1 <- c(beta1,beta[idx1[si2],si1*d+si2])
      }
    # print(x)
    # print(X[1:5,])
    
    side <- rep(1,n)
    for(si in 1:pd){side<-cbind(side,(X - (matrix( rep(x,each=n), n, d )))^si)}
    beta_expand <- as.vector(side %*% matrix(beta1,(pd*d+1),1)) # beta(X_{ki},x) defined as in eq(3) in the main context
    p <- 1/n 
    for(si in 1:d){p <- p* Epan(as.vector(side[,si+1])/h[si])/h[si]}
    valid_idx <- which(p>0) 
    # print(side[1:5,])
    
    if(length(valid_idx)>0){
      
      p <- p[valid_idx]; beta_expand <- beta_expand[valid_idx]; y1 <- y[valid_idx]
      side <- matrix(side[valid_idx,],length(valid_idx),pd*d+1) # matrix  of X(x)
      
      P[,i] <- -t(side)%*%matrix(q1(y1,beta_expand)*p,ncol=1)
      # if(length(valid_idx)>1){mm <- diag(q2(beta_expand) * p)}else{mm <- q2(beta_expand) * p}
      
      if(link!='identity'){

        for(ui in 1:(pd*d+1))
          for(uj in 1:(pd*d+1)){
            Q[ui,uj,i] <- -sum(q2(beta_expand) * p *side[,ui]*side[,uj])
            for(uk in 1:(pd*d+1)){
              R[ui,uj,uk,i] <- -0.5*sum(q3(beta_expand) * p *side[,ui]*side[,uj]*side[,uk])
            }
          }
        # Q[,,i] <- -t(side)%*%mm%*%side
        # for(uk in 1:(pd*d+1)){
        #   R[,,uk,i] <- -0.5*t(side)%*%mm%*%(side*matrix(rep(side[,uk],ncol(side)),ncol=ncol(side)))
        # }
        
      }else{
        
        for(ui in 1:(pd*d+1))
          for(uj in 1:(pd*d+1)){
            Q[ui,uj,i] <- -sum(q2(beta_expand) * p *side[,ui]*side[,uj])
          }
        # Q[,,i] <- -t(side)%*%mm%*%side
        # R <- 0
      }
      
    }
    
  }
  
  res <- list(P,Q,R)
  names(res) <- c('P', 'Q', 'R')
  return(res)
  
}

## compute_G_M_xi: Given beta(x), compute statistics G,M,xi
# should input UVQR evaluated at the closest candidate chain
compute_G_M_xi0 <- function(U_, V_, Q_, R_, P, Q, R, beta0, beta, n, N, h, pd){ 
  
  ### compute beta evaluated at all evaluation points
  idx_eval <- c()
  for(i in 1:d){
    idx_eval <- cbind(idx_eval, rep(rep(1:m,each=m^(i-1)),m^(d-i)))
  }
  
  ### compute G,M
  G <- matrix(0, pd*d+1, m^d)
  M <- array(0, dim = c(pd*d+1, pd*d+1, m^d))
  R1 <- rep(0, pd*d+1); R2 <- matrix(0,pd*d+1,pd*d+1)
  for(i in 1:m^d){
    
    idx1 <- idx_eval[i,]
    beta1 <- beta0; x <- c()
    for(si in 1:d){
      beta1 <- beta1 + beta[idx1[si],si]
      x <- c(x,eval_vec[idx1[si]])
    }
    for(si1 in 1:pd)
      for(si2 in 1:d){
        beta1 <- c(beta1,beta[idx1[si2],si1*d+si2])
      }
    
    if(link!='identity'){
      for(j1 in 1:(pd*d+1)){
        for(j2 in 1:(pd*d+1)){
          R2[j1,j2] <-  matrix(R_[j1,j2,,i], nrow=1) %*% matrix(beta1,ncol=1)
        }
        R1[j1] <- matrix(beta1,nrow=1) %*% matrix(R2[j1,], ncol=1)
      }
    }    
    
    G[,i] <- (N-n)/N * (U_[,i]+matrix(beta1,nrow=1) %*% V_[,,i]+R1) + n/N * P[,i]
    M[,,i] <- (N-n)/N * (V_[,,i]+2*R2) + n/N * Q[,,i]
    
  }
  
  ### compute xi0
  xi0 <- -mean(G[1,])/mean(M[1,1,])
  
  res <- list(G,M,xi0)
  names(res) <- c('G','M','xi0')
  return(res)
}

# should input UVQR evaluated at the closest candidate chain
backfitting <- 
  function(beta0,beta,U_,V_,Q_,R_, n, N, h, pd){
    
    delta_outer1 <- 1
    m_iter <- 0
    # Cartesian product of evaluation points
    idx_eval <- c()
    for(i in 1:d){
      idx_eval <- cbind(idx_eval, rep(rep(1:m,each=m^(i-1)),m^(d-i)))
    }
    
    # compute index used to update xi
    idx_xi <- matrix(0, pd+1, d)
    for(i in 1:d){idx_xi[,i] <- seq(i,(pd+1)*d,d)}
    
    # index extracted from G and M to update xi
    idx_M <- matrix(0, pd+1, d); idx_M[1,] <- 1
    idx_M[2:(pd+1),] <- idx_xi[1:pd,]+1
    
    ### outer loop
    while(delta_outer1 > delta_outer & m_iter<Max_iter){
      
      beta0_old <- beta0; beta_old <- beta
      # print(beta0)
      # print(beta)
      
      tt0 <- Sys.time()
      res_PQR <- compute_PQR(beta0, beta, n, h, pd)
      P <- res_PQR$P; Q <- res_PQR$Q; R <- res_PQR$R
      tt1 <- Sys.time()
      # print(paste0('PQR time=',difftime(tt1,tt0)))
      # print(P)
      
      # compute_G_M_xi0
      res <- compute_G_M_xi0(U_, V_, Q_, R_, P, Q, R, beta0, beta, n, N, h, pd)
      G <- res$G; M <- res$M; xi0 <- res$xi0
      # print('GM OK')
      
      ########### inner loop !!!!!!
      r <- 0
      delta_inner1 <- 1
      xi <- array(0, dim = c((pd+1)*d,m));xi_new <- array(0, dim = c((pd+1)*d,m))
      while(delta_inner1 > delta_inner & r < Max_iter){
        
        for(j in 1:d){
          
          for(i in 1:m){
            
            idx1 <- which(idx_eval[,j]==i)# index used to compute marginal integral
            
            eq_right <- -rowMeans(G[idx_M[,j],idx1])-xi0*rowMeans(M[idx_M[,j],1,idx1])
            for(ell in (1:d)[-j]){
              eq_right <- eq_right - rowMeans(sapply(idx1,function(ii){
                M[idx_M[,j],idx_M[,ell],ii]%*%matrix(xi_new[idx_xi[,ell],idx_eval[ii,ell]],ncol=1)
              }))
            }
            
            M_int <- rowMeans(M[idx_M[,j],idx_M[,j],idx1], dims = 2)
            M_inv <- try(solve(M_int), silent=TRUE)
            # dd <- c(dd, det(M_inv))
            if('try-error' %in% class(M_inv)){
              next
            }
            
            xi_new[idx_xi[,j],i] <- M_inv %*% eq_right
            
          }
          
        }
        
        r <- r + 1
        delta_inner1 <- sum(sapply(1:d, function(ii){sqrt(mean((xi[,ii]-xi_new[,ii])^2))}))
        
        print(paste('r=', r, 'delta_inner=', delta_inner1))
        if(delta_inner1>10/(1+(K>1))) break
        
        xi <- xi_new
      }
      
      if(delta_inner1>10/(1+(K>1))) break
      if(delta_outer1>10/(1+(K>1))) break
      
      
      # out update
      {
        ## compute the constants
        c <- rep(0,d)
        for(j in 1:d){
          q_xi_int <- rowMeans(sapply(idx_eval,function(ii){
            ((N-n)/N*Q_[1,idx_M[,j],ii] + n/N*Q[1,idx_M[,j],ii])*xi[idx_xi[,j],idx_eval[ii,j]]
          }))
          c[j] <- (mean(M[1,1,]))^{-1}*sum(q_xi_int)
        }
        
        ## update beta
        beta0 <- beta0_old + xi0 + sum(c)
        beta[,1:d] <- beta_old[,1:d] + t(xi[1:d,]) - matrix(rep(c,each=m),m,d)
        beta[,(d+1):((pd+1)*d)] <- beta_old[,(d+1):((pd+1)*d)] + t(xi[(d+1):((pd+1)*d),])
        
        # delta_out <- (max(abs(beta[,1:2]-beta_old[,1:2])))
        delta_outer1 <- sum(sapply(1:d, function(ii){sqrt(mean((beta[,ii]-beta_old[,ii])^2))}))
      }
      
      m_iter <- m_iter+1
      print(paste('m=_iter',m_iter,'delta_out=',delta_outer1))
    }
    print(paste('m=_iter',m_iter,'delta_out=',delta_outer1))
    
    res <- list(beta0, beta,  delta_outer1)
    names(res) <-c('beta0', 'beta', 'delta')
    return(res)
  }  

# should input UVQR evaluated at all candidates
update_stats <- function(U_,V_,Q_,R_, beta0,beta,eta,idx,n,N, pd,L){
  
  U <- rep(0, pd*d+1); P <- U
  Q <- matrix(0, pd*d+1, pd*d+1); V <- Q
  if(link!='identity'){
    R <- array(0, dim = c(pd*d+1, pd*d+1, pd*d+1))
  }else{
    R <- 0
  }
  
  idx_eval <- c()
  for(i in 1:d){
    idx_eval <- cbind(idx_eval, rep(rep(1:m,each=m^(i-1)),m^(d-i)))
  }
  
  for(l in 1:L){
    
    h <- eta[,l]
    res_PQR <- compute_PQR(beta0, beta, n, h, pd)
    P <- res_PQR$P; Q <- res_PQR$Q; R <- res_PQR$R
    
    for(i in 1:nrow(idx_eval)){
      
      
      idx1 <- idx_eval[i,]
      beta1 <- beta0; x <- c()
      for(si in 1:d){
        beta1 <- beta1 + beta[idx1[si],si]
        x <- c(x,eval_vec[idx1[si]])
      }
      for(si1 in 1:pd)
        for(si2 in 1:d){
          beta1 <- c(beta1,beta[idx1[si2],si1*d+si2])
        }
      
      if(link!='identity'){
        
        ### compute U, V
        for(ui in 1:(pd*d+1)){
          U[ui] <- P[ui,i] - sum(beta1*Q[ui,,i]) + matrix(beta1,1,pd*d+1) %*% R[ui,,,i] %*% matrix(beta1,pd*d+1,1)
          V[ui,] <- Q[ui,,i] - 2*matrix(beta1,1,pd*d+1)%*%R[ui,,,i]
        }
        
        ### update U,V,Q,R
        U_[,i,l] <- (N-n)/N * U_[,i,idx[l]] + n/N * U
        V_[,,i,l] <- (N-n)/N * V_[,,i,idx[l]] + n/N * V
        Q_[,,i,l] <- (N-n)/N * Q_[,,i,idx[l]] + n/N * Q[,,i]
        R_[,,,i,l] <- (N-n)/N * R_[,,,i,idx[l]] + n/N * R[,,,i]
        
      }else{
        
        ### compute U, V
        for(ui in 1:(pd*d+1)){
          U[ui] <- P[ui,i] - sum(beta1*Q[ui,,i])
        }
        
        ### update U,V,Q,R
        U_[,i,l] <- (N-n)/N * U_[,i,idx[l]] + n/N * U
        V_[,,i,l] <- (N-n)/N * V_[,,i,idx[l]] + n/N * Q[,,i]
        Q_[,,i,l] <- (N-n)/N * Q_[,,i,idx[l]] + n/N * Q[,,i]
        R_ <- 0
        
      }
      
      
    }
    
  }
  
  res <- list(U_,V_,Q_,R_)
  names(res) <- c('U_', 'V_', 'Q_', 'R_')
  return(res)
  
}


#' function to conduct the online smooth backfitting 
#' @param K the index of currently observed data block
#' @param X,y the observations
#' @param n the sample size of the block
#' @param m number of evaluation points
#' @param delta_inner,delta_outer permission error at training end for the inner and outer loops
#' @param Max_iter maximum number of iterations
#' @param band_select =TRUE if adopt data driven bandwidth selection; =FALSE if use pre-specified bandwidth
#' @param C1 pre-specified constant for bandwidth, is used when band_select=FALSE
#' @return a sequence of global statisstics and estimates
# online gam
ogam <- function(K, X, y, n, m, delta_inner, delta_outer, Max_iter, band_select, K_band, C1=NULL, L=1, L_theta=0, L_sigma=0, pd1=1, pd2=3){
  
  eval_vec <<- seq(0.05, 0.95, length.out = m)
  G <- rep(1,d)
  
  if(K==1){
    
    ## some initial values
    N <<- 0; h <<- rep(1,d)
    ## stored statistics for the main regression
    U_ <<- array(0, dim = c(pd1*d+1,m^d,L))
    V_ <<- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L)) 
    Q_ <<- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L))
    if(link!='identity'){
      R_ <<- array(0, dim = c(pd1*d+1,pd1*d+1,pd1*d+1,m^d,L))
    }else{
      R_ <<- 0
    }
    ## stored statistics for bandwidth selection
    if(band_select==TRUE){
      U_theta <<- array(0, dim = c(pd2*d+1,m^d,L_theta))
      V_theta <<- array(0, dim = c(pd2*d+1,pd2*d+1,m^d,L_theta))
      Q_theta <<- array(0, dim = c(pd2*d+1,pd2*d+1,m^d,L_theta))
      U_sigma <<- array(0, dim = c(pd1*d+1,m^d,L_sigma))
      V_sigma <<- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L_sigma))
      Q_sigma <<- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L_sigma))
      R_sigma <<- array(0, dim = c(pd1*d+1,pd1*d+1,pd1*d+1,m^d,L_sigma))
      R_theta <<- array(0, dim = c(pd2*d+1,pd2*d+1,pd2*d+1,m^d,L_theta))
      if(link!='identity'){
        R_sigma <<- array(0, dim = c(pd1*d+1,pd1*d+1,pd1*d+1,m^d,L_sigma))
        R_theta <<- array(0, dim = c(pd2*d+1,pd2*d+1,pd2*d+1,m^d,L_theta))
      }else{
        R_sigma <<- 0; R_theta <<- 0
      }
    }
    
  }  
  
  N <<- N + n
  if(K > K_band){band_select <- FALSE}
  if(K <= K_band & band_select==FALSE){Ch <<- C1}
  
  # bandwidth selection
  if(band_select==TRUE){
    
    h_theta <- G * N^(-1/7); h_sigma <- G * N^(-1/5)
    eta_theta <- sapply(1:L_theta, function(l){((L_theta-l+1)/L_theta)^(1/7) * h_theta})#dim: d*L
    eta_sigma <- sapply(1:L_sigma, function(l){((L_sigma-l+1)/L_sigma)^(1/5) * h_sigma})
    
    if(K==1){
      
      C <- c()
      
      # generate initial estimate and bandwidth
      initial_res <- initial(X,y,h_theta,eval_vec,pd2)
      beta0_theta_old <- initial_res$beta0; beta_theta_old <- initial_res$beta
      beta0_theta <- initial_res$beta0; beta_theta <- initial_res$beta
      idx_theta <- 1:L_theta
      centrds_theta <<- eta_theta
      rm(initial_res)
      
      initial_res <- initial(X,y,h_sigma,eval_vec,pd1)
      beta0_sigma_old <- initial_res$beta0; beta_sigma_old <- initial_res$beta
      beta0_sigma <- initial_res$beta0; beta_sigma <- initial_res$beta
      idx_sigma <- 1:L_sigma
      centrds_sigma <<- eta_sigma
      rm(initial_res)
      
    }else{
      
      # update index for pseudo-bandwidth and the centroids
      idx_theta <- sapply(1:L_theta,function(l){
        which.min(abs(eta_theta[1,l] - centrds_theta[1,]))
      })# idx of all d dimensions are the same, it is sufficient to compuute the first dim
      for(i in 1:d){
        centrds_theta[i,] <<- (centrds_theta[i,idx_theta] * (N-n) + eta_theta[i,] * n) / N
      }

      idx_sigma <- sapply(1:L_sigma,function(l){
        which.min(abs(eta_sigma[1,l] - centrds_sigma[1,]))
      })
      for(i in 1:d){
        centrds_sigma[i,] <<- (centrds_sigma[i,idx_sigma] * (N-n) + eta_sigma[i,] * n) / N
      }
    }
    
    # backfitting
    {
      # tt0 <- Sys.time()
      res_beta <- backfitting(beta0_theta, beta_theta,U_theta[,,idx_theta[1]],V_theta[,,,idx_theta[1]],
                              Q_theta[,,,idx_theta[1]],R_theta[,,,,idx_theta[1]], n, N, h_theta, pd2)
      if(K==1){rm(beta0_theta, beta_theta)}
      if(res_beta$delta < delta_outer){
        beta0_theta <<- res_beta$beta0
        beta_theta <<- res_beta$beta
      }
      rm(res_beta)
      # tt1 <- Sys.time()
      # print(paste0('theta backfitting:',round(difftime(tt1,tt0,units = 'mins'),3),'mins'))
      
      # tt0 <- Sys.time()
      res_beta <- backfitting(beta0_sigma, beta_sigma,U_sigma[,,idx_sigma[1]],V_sigma[,,,idx_sigma[1]],
                              Q_sigma[,,,idx_sigma[1]],R_sigma[,,,,idx_sigma[1]], n, N, h_sigma, pd1)
      if(K==1){rm(beta0_sigma, beta_sigma)}
      if(res_beta$delta < delta_outer){
        beta0_sigma <<- res_beta$beta0
        beta_sigma <<- res_beta$beta
      }
      rm(res_beta)
      # tt1 <- Sys.time()
      # print(paste0('sigma backfitting:',round(difftime(tt1,tt0,units = 'mins'),3),'mins'))
    }
    
    # compute the constants
    {
      beta_hat <- matrix(0,nrow(X),d); sigma2 <- rep(0,d)
      for(i in 1:d){
        beta_hat[,i] <- predict(smooth.spline(eval_vec,beta_sigma[,i]),X[,i])$y
      }
      if(link=='log'){
        for(i in 1:d){
          sigma2[i] <- mean(exp(beta_hat[,(1:d)[-i]]))*mean(exp(-beta_sigma[,i]))*exp(-beta0_sigma)
        }
      }
      if(link=='identity'){
        sigma2 <- mean((y-rowSums(beta_hat[,(1:d)])-beta0_sigma)^2)
      }
      if(link=='logit'){
        sigma2 <- 4*mean(1/(1-exp(-beta0_sigma-rowSums(beta_hat)))*(1-1/(1+exp(-beta0_sigma-rowSums(beta_hat)))))
      }
      Ch <<- (15 * sigma2/colMeans(beta_theta[,(2*d+1):(3*d)]^2))^(1/5)
    }
    
    # update statistics
    {
      res_update <- update_stats(U_theta,V_theta,Q_theta,R_theta, beta0_theta,beta_theta,eta_theta,
                                 idx_theta,n,N,pd2,L_theta)
      U_theta <<- res_update$U_
      V_theta <<- res_update$V_
      Q_theta <<- res_update$Q_
      R_theta <<- res_update$R_
      rm(res_update)
      
      res_update <- update_stats(U_sigma,V_sigma,Q_sigma,R_sigma, beta0_sigma,beta_sigma,eta_sigma,
                                 idx_sigma,n,N,pd1,L_sigma)
      U_sigma <<- res_update$U_
      V_sigma <<- res_update$V_
      Q_sigma <<- res_update$Q_
      R_sigma <<- res_update$R_
      rm(res_update)
    }
    
  }
  
  # compute bandwidth and candidates
  {
    h <<- sapply(1:d, function(i){min(Ch[i]*N^(-1/5),h[i])})
    eta <- sapply(1:L, function(l){((L-l+1)/L)^(1/5) * h}) #dim = d*L
  }
  

  # initial parametric estimate and combination rule
  if(K==1){
    
    initial_res <- initial(X,y,h,eval_vec,pd1)
    beta0_est <<- initial_res$beta0; beta_est <<- initial_res$beta
    
    idx <- 1:L
    centrds <<- eta
    
  }else{
    
    idx<-sapply(1:L,function(l){
      which.min(abs(eta[1,l] - centrds[1,]))
    })
    for(i in 1:d){
      centrds[i,] <- (centrds[i,idx] * (N-n) + eta[i,] * n) / N
    }
    centrds <<- centrds
    
  }

  t0 <<- Sys.time()
  # backfitting
  {
    res_beta <- backfitting(beta0_est, beta_est,U_[,,idx[1]],V_[,,,idx[1]],Q_[,,,idx[1]],R_[,,,,idx[1]], n, N, h, pd1)
    if(res_beta$delta < delta_outer){
      beta0_est <<- res_beta$beta0; beta_est <<- res_beta$beta
      print(paste0('backfitting modified : delta=', res_beta$delta))
    }
    rm(res_beta)
  }
  # print('backfitting OK')
  
  # update statistics
  {
    res_update <- update_stats(U_,V_,Q_,R_, beta0_est, beta_est, eta, idx, n, N, pd1, L)
    U_ <<- res_update$U_
    V_ <<- res_update$V_
    Q_ <<- res_update$Q_
    R_ <<- res_update$R_
    rm(res_update)
  }
  t1 <<- Sys.time()
  # print('update statistics OK')
  
}

# batch gam
gam <- function(K, X, y, n, m, delta_inner, delta_outer, Max_iter, band_select, K_band, C1=NULL, L=1, L_theta=1, L_sigma=1, pd1=1, pd2=3){
  
  eval_vec <<- seq(0.05, 0.95, length.out = m)
  G <- rep(1,d)
  
  {
    ## some initial values
    N <<- 0; h <<- rep(1,d)
    ## stored statistics for the main regression
    U_ <- array(0, dim = c(pd1*d+1,m^d,L))
    V_ <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L)) 
    Q_ <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L))
    if(link!='identity'){
      R_ <- array(0, dim = c(pd1*d+1,pd1*d+1,pd1*d+1,m^d,L))
    }else{
      R_ <- 0
    }
    ## stored statistics for bandwidth selection
    if(band_select==TRUE){
      U_theta <- array(0, dim = c(pd2*d+1,m^d,L_theta))
      V_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,m^d,L_theta))
      Q_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,m^d,L_theta))
      U_sigma <- array(0, dim = c(pd1*d+1,m^d,L_sigma))
      V_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L_sigma))
      Q_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,m^d,L_sigma))
      R_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,pd1*d+1,m^d,L_sigma))
      R_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,pd2*d+1,m^d,L_theta))
      if(link!='identity'){
        R_sigma <- array(0, dim = c(pd1*d+1,pd1*d+1,pd1*d+1,m^d,L_sigma))
        R_theta <- array(0, dim = c(pd2*d+1,pd2*d+1,pd2*d+1,m^d,L_theta))
      }else{
        R_sigma <- 0; R_theta <- 0
      }
    }
    
  }  
  
  N <<- N + n
  if(band_select==FALSE){Ch <<- C1}
  
  # bandwidth selection
  if(band_select==TRUE){
    
    h_theta <- G * N^(-1/7); h_sigma <- G * N^(-1/5)
    eta_theta <- sapply(1:L_theta, function(l){((L_theta-l+1)/L_theta)^(1/7) * h_theta})#dim: d*L
    eta_sigma <- sapply(1:L_sigma, function(l){((L_sigma-l+1)/L_sigma)^(1/5) * h_sigma})
    
    {
      C <- c()
      
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
      res_beta <- backfitting(beta0_theta, beta_theta,U_theta[,,idx_theta[1]],V_theta[,,,idx_theta[1]],
                              Q_theta[,,,idx_theta[1]],R_theta[,,,,idx_theta[1]], n, N, h_theta, pd2)
      if(K==1){rm(beta0_theta, beta_theta)}
      if(res_beta$delta < delta_outer){
        beta0_theta <<- res_beta$beta0
        beta_theta <<- res_beta$beta
      }
      rm(res_beta)

      res_beta <- backfitting(beta0_sigma, beta_sigma,U_sigma[,,idx_sigma[1]],V_sigma[,,,idx_sigma[1]],
                              Q_sigma[,,,idx_sigma[1]],R_sigma[,,,,idx_sigma[1]], n, N, h_sigma, pd1)
      if(K==1){rm(beta0_sigma, beta_sigma)}
      if(res_beta$delta < delta_outer){
        beta0_sigma <<- res_beta$beta0
        beta_sigma <<- res_beta$beta
      }
      rm(res_beta)
    }
    
    # compute the constants
    {
      beta_hat <- matrix(0,nrow(X),d); sigma2 <- rep(0,d)
      for(i in 1:d){
        beta_hat[,i] <- predict(smooth.spline(eval_vec,beta_sigma[,i]),X[,i])$y
      }
      if(link=='log'){
        for(i in 1:d){
          sigma2[i] <- mean(exp(beta_hat[,(1:d)[-i]]))*mean(exp(-beta_sigma[,i]))*exp(-beta0_sigma)
        }
      }
      if(link=='identity'){
        sigma2 <- mean((y-rowSums(beta_hat[,(1:d)])-beta0_sigma)^2)
      }
      if(link=='logit'){
        sigma2 <- 4*mean(1/(1-exp(-beta0_sigma-rowSums(beta_hat)))*(1-1/(1+exp(-beta0_sigma-rowSums(beta_hat)))))
      }
      Ch <<- (15 * sigma2/colMeans(beta_theta[,(2*d+1):(3*d)]^2))^(1/5)
    }
    
  }
  
  # compute bandwidth and candidates
  h <<- Ch*N^(-1/5)

  # initial parametric estimate and combination rule
  initial_res <- initial(X,y,h,eval_vec,pd1)
  beta0_est <<- initial_res$beta0; beta_est <<- initial_res$beta

  # backfitting
  {
    res_beta <- backfitting(beta0_est, beta_est,U_[,,idx[1]],V_[,,,idx[1]],Q_[,,,idx[1]],R_[,,,,idx[1]], n, N, h, pd1)
    if(res_beta$delta < delta_outer){
      beta0_est <<- res_beta$beta0; beta_est <<- res_beta$beta
      print(paste0('backfitting modified : delta=', res_beta$delta))
    }
    rm(res_beta)
  }

}

