#' helper functions

###################  Part 1 functions for tuning parameters  ######################
#' @param t the time points 
#' @param y observed values in time points t
#' @param Cseq the candidate sequence of coefficient C
#' @param C the coefficient to calculate the basis number
#' @param hseq the candidate sequence of h
#' @param ext the value of extension margin
#' @param domain the t domain 
#' @return a optimal value of h and coefficient of rho_n
tune <- function(t,y,Cseq,C,hseq,ext,domain)
{
  Kfold <- ifelse(length(t)<20,yes=length(t),no=5)
  n=length(t)
  weig=1/n
  cvo <- cv.partition(length(t),Kfold)
  err <- lapply(1:cvo$num.test.sets,
                function(k){
                  trIdx <- cvo$training[[k]]
                  teIdx <- cvo$test[[k]]
                  sapply(hseq,function(h){
                    aux.mat <- compute.aux.matrices(q=floor(n^h/C),t=t[trIdx],y=y[trIdx],domain=extend(domain,ext))
                    sapply(Cseq,function(c){
                      ahat <- solve((c/n^((7*h+1)/2))*aux.mat$W+aux.mat$H*weig,aux.mat$G*weig)
                      sse(ahat=ahat,ext=ext,domain=domain,t=t[teIdx],y=y[teIdx])})},simplify=TRUE)
                })
  err <- as.matrix(Reduce('+',err))
  I <- which(err==min(err),arr.ind=T)
  hopt=hseq[I[1,2]]
  optC=Cseq[I[nrow(I),1]]
  R=list(hopt,optC)
  return(R)
}

#' @param n sample size
#' @param K number of partition
#' @return a list of data
cv.partition <- function(n,K)
{
  if(K <= 1) stop('K must be larger than 1')
  I <- sample.int(n)
  a <- n %% K
  b <- (n-a)/K
  trunk.size <- c(rep(b+1,a),rep(b,K-a))
  start.idx <- c(1,1+cumsum(trunk.size[1:(K-1)]))
  end.idx <- start.idx + trunk.size - 1
  test <- lapply(1:K,function(k){
    I[start.idx[k]:end.idx[k]]
  })
  training <- lapply(1:K,function(k){I[-(start.idx[k]:end.idx[k])]})
  return(list(training=training,test=test,num.test.sets=length(test)))
}

#' function to extend the domain 
#' @param ext the value of extension margin
#' @param domain the t domain
#' @return a extended domain
extend <- function(domain,ext)
{
  return(c(domain[1]-(domain[2]-domain[1])*ext,
           domain[2]+(domain[2]-domain[1])*ext))
}

#' function to predict the regression function in the extended domain
#' @param tobs the time points 
#' @param ahat the estimated coefficients of basis expansion
#' @param ext the value of extension margin
#' @param domain the t domain
#' @return a vector
pred <- function(tobs,ahat,ext,domain)
{
  if(ext > 0) D <- c(domain[1]-(domain[2]-domain[1])*ext,
                     domain[2]+(domain[2]-domain[1])*ext)
  else D <- domain
  B <- evaluate.basis(K=length(ahat),grid=tobs,domain=D)
  return(c(B %*% ahat))
}

#' function to calculate error sum of squares (SSE)
#' @param ahat the estimated coefficients of basis expansion
#' @param ext the value of extension margin
#' @param domain the t domain
#' @param t time points of observations 
#' @param y the observations in time points t
#' @return the value of SSE
sse <- function(ahat,ext,domain,t,y)
{
  E <- sapply(1:length(t),function(i){
    tmp <- pred(t[[i]],ahat,ext,domain)
    sum((tmp-y[[i]])^2)
  })
  return(sum(E))
}

################### Part 2 functions for generating basis functions (Fourier basis)  ######################

#' function to generate Fourier basis 
#' @param K the number of basis functions
#' @param m the number of equispaced points on domain
#' @param domain the domain on which basis functions are defined
#' @param grid a vector specifying the time points to evaluate the basis functions. If grid is supplied, then m is ignored
#' @return a m by K matrix, where rows index basis functions while columns index points in the grid.
evaluate.basis <- function(K,
                           m = 51,
                           domain = c(0,1),
                           grid = seq(domain[1], domain[2], length.out = m))
{
  stopifnot(is.numeric(K) && length(K) == 1 && K > 0)
  m <- length(grid)
  a <- domain[1]
  b <- domain[2]
  x <- 0
  y <- 1
  alpha <- (b*x-a*y)/(x-y)
  beta <- (a-b)/(x-y)
  pts <- (grid - alpha) / beta
  res <- sapply(1:K, function(k)
    if (k == 1) rep(1, m)
    else if (k %% 2 == 0)  sqrt(2) * cos(k * pi * pts)
    else sqrt(2) * sin((k - 1) * pi * pts))
  res <- res / sqrt(beta) 
  res <- matrix(res, ncol = K) 
  return(res)
}

#' function to generate the derivative functions of basis 
#' @param K the number of basis functions
#' @param m the number of equispaced points on domain
#' @param domain the domain on which basis functions are defined
#' @param grid a vector specifying the time points to evaluate the basis functions. If grid is supplied, then m is ignored
#' @param r the order of derivative
#' @return a m by K matrix
deriv.fourier <- function(K,
                          m = 51,
                          domain = c(0,1),
                          grid = seq(domain[1], domain[2], length.out = m),
                          r=2)
{
  if (min(grid) < domain[1] || max(grid) > domain[2])
    stop('some points in grid are outside of domain.')
  m <- length(grid)
  V <- matrix(0,m,K)
  L <- domain[2] - domain[1]
  a <- domain[1]
  for(k in 1:K)
  {
    if(k==1)
    {
      if(r==0) V[,k] <- 1
    }
    else if(k %% 2 == 0) 
    {
      s <- 1
      if(r %% 4 == 1 || r %% 4 == 2) s <- -1
      if(r %% 2 == 1)
        V[,k] <- s * sqrt(2/L) * (k*pi/L)^r * sin(k*pi*(grid-a)/L)
      else
        V[,k] <- s * sqrt(2/L) * (k*pi/L)^r * cos(k*pi*(grid-a)/L)
    }
    else 
    {
      s <- 1
      if(r %% 4 == 2 || r %% 4 == 3) s <- -1
      if(r %% 2 == 1)
        V[,k] <- s * sqrt(2/L) * ((k-1)*pi/L)^r * cos((k-1)*pi*(grid-a)/L)
      else
        V[,k] <- s * sqrt(2/L) * ((k-1)*pi/L)^r * sin((k-1)*pi*(grid-a)/L)
    }
  }
  return(V)
}

#' function to generate regular points in the domain
#' @param m the number of points
#' @param domain the t domain
#' @param h the margin of the first and last points to the boundaries of the domain
#' @return a sequence of time points in the domain
regular.grid <- function(m=100,domain=c(0,1),h=1/(2*m))
{
  seq(domain[1]+h,domain[2]-h,length.out=m)
}

################### Part 3 functions for regression function estimation  ######################

#' function to calculate the summary statistics for estimation
#' @param t time points of observations 
#' @param y observed values in time points t
#' @param q the number of basis functions
#' @param domain the t domain
#' @return a list containing three summary statistics
compute.aux.matrices <- function(t,y,q,domain=c(0,1))
{
   n <- length(t)
   B=evaluate.basis(K=q,domain=domain,grid=t)
   H=t(B)%*%B
   G=t(B)%*%y
   M <- 1000
   pts <- regular.grid(m=1000,domain=domain)
   W <- deriv.fourier(K=q,grid=pts,r=2,domain=domain)
   W <- t(W) %*% W / M
   return(list(W=W,G=G,H=H))
}

#' function to calculate the summary statistic G
#' @param t time points  
#' @param y observed values in time points t
#' @param q the number of basis functions
#' @param domain the t domain
#' @return summary statistic of G
computeG <- function(t,y,q,domain=c(0,1))
{
  G=t(evaluate.basis(K=q,domain=domain,grid=t))%*%y
  return(G)
}

#' function to calculate the penalty term W
#' @param q the number of basis functions
#' @param domain the t domain
#' @return penalty matrix W
computeW <- function(q,domain)
{
  M <- 1000
  pts <- regular.grid(m=1000,domain=domain)
  W <- deriv.fourier(K=q,grid=pts,r=2,domain=domain)
  W <- t(W) %*% W / M
  return(W)
}

#' ################### Part 4 predict and evaluate  ##################################################

#' function to make predictions
#' @param regfunc.obj the object obtained by calling \code{reg.func}
#' @param newt a vector containing time points of observations to be evaluated.
#' @return a vector of estimation
predict.regfunc <- function(regfunc.obj,newt)
{
  domain <- regfunc.obj$domain
  ext <- regfunc.obj$ext
  if(ext > 0)
    D <- c(domain[1]-(domain[2]-domain[1])*ext,
           domain[2]+(domain[2]-domain[1])*ext)
  else D <- domain
  B <- evaluate.basis(K=length(regfunc.obj$ahat),grid=newt,domain=D)
  return(c(B %*% regfunc.obj$ahat))
}

#' other helper functions to verify phase transition and compare tuning methods between full data driven method and semi data driven method 

############################### full data driven method ################################
#' function to select h
#' @param t time points of observations 
#' @param y observed values in time points t 
#' @param hseq the candidate sequence of h
#' @param rho the candidate sequence of rho
#' @param ext the value of extension margin
#' @param domain the t domain
#' @param C the coefficient to calculate the basis number
#' @return optimal h
tuneh <- function(t,y,hseq,rhoseq,ext,domain,C)
{
  Kfold <- ifelse(length(t)<20,yes=length(t),no=5)
  n=length(t)
  weig=1/n
  cvo <- mcfda::cv.partition(length(t),Kfold)
  err <- lapply(1:cvo$num.test.sets,
                function(k){
                  trIdx <- cvo$training[[k]]
                  teIdx <- cvo$test[[k]]
                  sapply(hseq,function(h){
                    aux.mat <- compute.aux.matrices(q=floor(n^h/C),t=t[trIdx],y=y[trIdx],domain=extend(domain,ext))
                    sapply(rhoseq,function(ri){
                      ahat <- solve(ri*aux.mat$W+aux.mat$H*weig,aux.mat$G*weig)
                      sse(ahat=ahat,ext=ext,domain=domain,t=t[teIdx],y=y[teIdx])})},simplify=TRUE)})
  err <- as.matrix(Reduce('+',err))
  I <- which(err==min(err),arr.ind=T)
  hopt=hseq[I[1,2]]
  optrho=rhoseq[I[nrow(I),1]]
  R=list(hopt,optrho)
  return(R)
}

#' function to select coefficient of rho
#' @param t time points of observations (a list form)
#' @param y observed values in time points t (a list form)
#' @param rho the candidate sequence of rho
#' @param ext the value of extension margin
#' @param q the number of basis functions
#' @param domain the t domain
#' @return a optimal rho
tune1p <- function(t,y,rhoseq,ext,q,domain)
{
  if(is.null(rhoseq)){rhoseq=10^seq(-8,2,length.out=21)}
  Kfold <- ifelse(length(t)<20,yes=length(t),no=5)
  n=length(t)
  weig=1/n
  cvo <- mcfda::cv.partition(length(t),Kfold)
  err <- lapply(1:cvo$num.test.sets,
                function(k){
                  trIdx <- cvo$training[[k]]
                  teIdx <- cvo$test[[k]]
                  aux.mat <- compute.aux.matrices(q=q,t=t[trIdx],y=y[trIdx],domain=extend(domain,ext))
                  sapply(rhoseq,function(ri){
                    ahat <- solve(ri*aux.mat$W+aux.mat$H*weig,aux.mat$G*weig)
                    sse(ahat=ahat,ext=ext,domain=domain,t=t[teIdx],y=y[teIdx])},simplify=TRUE)})
  err <- as.matrix(Reduce('+',err))
  I <- which(err==min(err),arr.ind=T)
  optrho=rhoseq[I[nrow(I),1]]
  return(optrho)
}

#' tune function for the batches except the first batch
tunep <- function(H,G,nq,t,y,q,rhoseq,ext,domain)
{
  q1=length(nq)
  nq1=1/(nq-length(t))
  G2=diag(nq1)%*%G[1:q1]
  W=computeW(q=q1,domain=extend(domain,ext))
  err= sapply(rhoseq,function(ri){
    ahat <- solve(ri*W+H[1:q1,1:q1],G2)
    sse(ahat=ahat,ext=ext,domain=domain,t=t,y=y)})
  I <- which(err==min(err))
  if(length(I)>1){I=I[1]}
  rhoopt=rhoseq[I]
  return(rhoopt)
}

#' batch run code with full data driven method to tune parameters
batchrunfull=function(data_batch,ext=0,c0=0.5,domain=c(0,1),h=1/5,q0=5,C=1/2,m){
  ln=1000 #the number of discrete time points to predict density function
  EV=100  #the number of discrete time points to predict regression function
  tnew=seq(domain[1],domain[2],length=ln)
  tnew1=seq(domain[1],domain[2],length=EV)
  K=length(data_batch)  
  b=length(data_batch[[1]][[1]])
  n=K*b
  qstart=q0
  qend=floor(n^h/C)
  inc=qend-qstart
  Q=seq(qstart,qend,by=1)
  
  S=rep(0,inc)
  tauseq=rep(0,inc)
  B=rep(0,inc)
  tau=rep(0,inc)
  for( j in 1:inc)
  {
    S[j]= floor((C*Q[j+1])^(1/h))
    B[j]=ceiling(S[j]/b)
    tauseq[j]=floor(c0*S[j])
    tau[j]=ceiling(tauseq[j]/b)
  }
  
  lostn=rep(0,Q[length(Q)]) 
  nq=rep(0,Q[length(Q)]) 
  rhoseq=rep(0,K) 
  rhoopt=rep(0,K) 
  coe=seq(-2,0,length=20) 
  regest=matrix(0,K,EV) 
  
  lpreinc=length(which(tau==1))
  lpreB=length(which(B==1))
  preB=B
  if(lpreB>0){preB=B[-which(B==1)]}
  n=0
  tauseq=c(0,tauseq)
  maintime=c()
  querytime=c()
  #' Update
  {
    print(paste('K=',1))
    t=data_batch[[1]][,1]
    y=data_batch[[1]][,2]
    # maintain
    t0=Sys.time()
    q=qstart+lpreB
    if(lpreinc>0){q=qstart+lpreinc}
    n=n+length(t)
    G=computeG(t=t,y=y,q=q,domain=extend(domain,ext))
    G1=G[1:(qstart+lpreB)]
    nq[1:q]=n
    nq1=1/nq[1:(qstart+lpreB)]
    G2=diag(nq1)%*%G1 
    p=q
    theta=apply(evaluate.basis(K=p,domain=domain,grid=unlist(t)),2,sum)
    theta2=theta/n
    t1=Sys.time()
    # query
    # calculate H
    R1=list(ahat=theta2,ext=0,domain=domain)
    fhat=predict.regfunc(R1,tnew)*(domain[2]-domain[1])
    Phi=evaluate.basis(K=q,domain=extend(domain,ext),grid=tnew)
    H=t(Phi)%*%diag(fhat)%*%Phi/ln
    penmat=computeW(q=q,domain=extend(domain,ext))
    rhoopt[1]=tune1p(t=t,y=y,rho=NULL,ext=ext,q=q,domain=domain)
    rhoseq[1]=rhoopt[1]
    ahat= solve(H[1:(qstart+lpreB),1:(qstart+lpreB)]+rhoopt[1]*penmat[1:(qstart+lpreB),1:(qstart+lpreB)],G2)
    R=list(ahat=ahat,penmat=penmat,ext=ext,rho=rhoopt[1],G=G,domain=domain)
    class(R) <- 'regfunc'
    regest[1,]=predict.regfunc(R,tnew1)
    t2=Sys.time()
    q=qstart+lpreB
    maintime=c(maintime,difftime(t1,t0,units = 'secs'))
    querytime=c(querytime,difftime(t2,t1,units = 'secs'))
  }
  if(2%in%preB)
  {
    q=q+length(which(preB==2))
    preB=preB[-which(preB==2)]
  }
  for (i in 2:(preB[1]-1)){
    print(paste('K=',i))
    t=data_batch[[i]][,1] 
    y=data_batch[[i]][,2]
    # maintain
    t0 <- Sys.time()
    if(i %in% tau)
    {
      G=c(G,rep(0,length(which(tau==i))))
      theta=c(theta,rep(0,length(which(tau==i))))
      lostn[length(G)]=n
    }
    n=n+length(t)
    nq[1:q]=n-lostn[1:q]
    nq1=1/nq[1:q]
    p=q
    theta1=apply(evaluate.basis(K=length(theta),domain=domain,grid=unlist(t)),2,sum)
    theta=theta+theta1
    theta2=diag(nq1)%*%(theta[1:p])
    G=G+computeG(t=t,y=y,q=length(G),domain=extend(domain,ext))
    G2=diag(nq1)%*%G[1:q]
    t1 <- Sys.time()
    # query
    R1 <- list(ahat=theta2,ext=0,domain=domain)
    fhat=predict.regfunc(R1,tnew)*(domain[2]-domain[1])
    Phi=evaluate.basis(K=length(G),domain=extend(domain,ext),grid=tnew)
    H=t(Phi)%*%diag(fhat)%*% Phi/ln
    penmat=computeW(q=q,domain=extend(domain,ext))
    rhoseq[i]=rhoopt[i-1]
    rhoseq1=sapply(coe,function(s){res=10^s*rhoseq[i]})
    rhoopt[i]=tunep(H=H,G=G,nq=nq[1:q],t=t,y=y,q=length(G),rhoseq=rhoseq1,ext=ext,domain=domain)
    ahat= solve(H[1:q,1:q]+rhoopt[i]*penmat,G2)
    R=list(ahat=ahat,penmat=penmat,rho=rhoopt[i],ext=ext,G=G,domain=domain)
    class(R) <- 'regfunc'
    regest[i,]=predict.regfunc(R,tnew1)
    t2 <- Sys.time()
    maintime=c(maintime,difftime(t1,t0,units = 'secs'))
    querytime=c(querytime,difftime(t2,t1,units = 'secs'))
  }
  q=q+1
  if(length(which(preB==3))>1)
  {
    preB=preB[-1]
    q=q+length(which(preB==3))-1
  }
  for(j in 1:(length(preB)-1)){
    for(i in (preB[j]):(preB[j+1]-1)) {
      print(paste('K=',i))
      t=data_batch[[i]][,1]  
      y=data_batch[[i]][,2] 
      # maintain
      t0 <- Sys.time()
      if(i %in% tau)
      {
        G=c(G,rep(0,length(which(tau==i))))
        theta=c(theta,rep(0,length(which(tau==i))))
        lostn[length(G)]=n
      }
      n=n+length(t)
      nq[1:q]=n-lostn[1:q]
      nq1=1/nq[1:q]
      p=q
      theta1=apply(evaluate.basis(K=length(theta), domain=domain,grid=unlist(t)),2,sum)
      theta=theta+theta1
      theta2=diag(nq1)%*%(theta[1:p])
      G=G+computeG(t=t,y=y,q=length(G),domain=extend(domain,ext))
      G2=diag(nq1)%*%G[1:q]
      t1 <- Sys.time()
      # query
      R1 <- list(ahat=theta2,ext=0,domain=domain)
      fhat=predict.regfunc(R1,tnew)*(domain[2]-domain[1])
      Phi=evaluate.basis(K=length(G), domain=extend(domain,ext),grid=tnew)
      H=t( Phi)%*%diag(fhat)%*% Phi/ln
      penmat=computeW(q=q,domain=extend(domain,ext))
      rhoseq[i]=rhoopt[i-1]
      rhoseq1=sapply(coe,function(s){res=10^s*rhoseq[i]})
      rhoopt[i]=tunep(H=H,G=G,nq=nq[1:q],t=t,y=y,q=length(G),rhoseq=rhoseq1,ext=ext,domain=domain)
      ahat= solve(H[1:q,1:q]+rhoopt[i]*penmat,G2)
      R=list(ahat=ahat,penmat=penmat,rho=rhoopt[i],ext=ext,G=G,domain=domain)
      class(R) <- 'regfunc'
      regest[i,]=predict.regfunc(R,tnew1)
      t2 <- Sys.time()
      maintime=c(maintime,difftime(t1,t0,units = 'secs'))
      querytime=c(querytime,difftime(t2,t1,units = 'secs'))
    }
    q=q+1
  }
  if(i<K){
    for(i in (B[length(B)]):K) {
      print(paste('K=',i))
      t=data_batch[[i]][,1]  
      y=data_batch[[i]][,2] 
      t0 <- Sys.time()
      n=n+length(t)
      nq[1:q]=n-lostn[1:q]
      nq1=1/nq[1:q]
      p=q
      theta1=apply(evaluate.basis(K=length(theta), domain=domain,grid=unlist(t)),2,sum)
      theta=theta+theta1
      theta2=diag(nq1)%*%(theta[1:p])
      G=G+computeG(t=t,y=y,q=length(G),domain=extend(domain,ext))
      G2=diag(nq1)%*%G[1:q]
      t1 <- Sys.time()
      # query
      R1 <- list(ahat=theta2,ext=0,domain=domain)
      fhat=predict.regfunc(R1,tnew)*(domain[2]-domain[1])
      Phi=evaluate.basis(K=length(G), domain=extend(domain,ext),grid=tnew)
      H=t(Phi)%*%diag(fhat)%*%Phi/ln
      penmat=computeW(q=q,domain=extend(domain,ext))
      rhoseq[i]=rhoopt[i-1]
      rhoseq1=sapply(coe,function(s){res=10^s*rhoseq[i]})
      rhoopt[i]=tunep(H=H,G=G,nq=nq[1:q],t=t,y=y,q=length(G),rhoseq=rhoseq1,ext=ext,domain=domain)
      ahat= solve(H[1:q,1:q]+rhoopt[i]*penmat,G2)
      R <- list(ahat=ahat,penmat=penmat,rho=rhoopt[i],ext=ext,G=G,domain=domain)
      class(R) <- 'regfunc'
      regest[i,]=predict.regfunc(R,tnew1)
      t2 <- Sys.time()
      maintime=c(maintime,difftime(t1,t0,units = 'secs'))
      querytime=c(querytime,difftime(t2,t1,units = 'secs'))
    }
  }
  RMISE=rep(0,K)
  mut=m(tnew1)
  for(i in 1:K){
    RMISE[i]=RMISEfunc(tnew1,regest[i,],mut)}
  R=list(RMISE,maintime,querytime)
  return(R)
}

############################### phase transition of L2 convergence rate  ################################
#'function to calculate the estimation error with L2-norm
#'@param a1 true coefficient vector
#'@param a2 estimated coefficient vector
#'@param K the basis number
#'@param domain the t domain
#'@return estimation error
estl2=function(a1,a2,K,domain)
{
  pts=seq(domain[1],domain[2],length=1000)
  a <- domain[1]
  b <- domain[2]
  x <- 0
  y <- 1
  alpha <- (b*x-a*y)/(x-y)
  beta <- (a-b)/(x-y)
  pts <- (pts - alpha) / beta
  res <- sapply(1:K, function(k)
    if (k == 1) rep(1, 1000)
    else if (k %% 2 == 0)  sqrt(2) * cos(k * pi * pts)
    else sqrt(2) * sin((k - 1) * pi * pts))
  res <- res / sqrt(beta) # normalization
  res <- matrix(res, ncol = K) 
  errf=(res%*%(a1-a2))^2
  result=sqrt(pracma::trapz(pts,errf))
  return(result)
}

#'function to calculate the approximation error with L2-norm
#'@param a1 true coefficient vector
#'@param m true regression function
#'@param K the basis number
#'@param domain the t domain
#'@return approximation error
appl2=function(a1,m,K,domain)
{
  pts=seq(domain[1],domain[2],length=1000)
  a <- domain[1]
  b <- domain[2]
  x <- 0
  y <- 1
  alpha <- (b*x-a*y)/(x-y)
  beta <- (a-b)/(x-y)
  pts <- (pts - alpha) / beta
  res <- sapply(1:K, function(k)
    if (k == 1) rep(1, 1000)
    else if (k %% 2 == 0)  sqrt(2) * cos(k * pi * pts)
    else sqrt(2) * sin((k - 1) * pi * pts))
  res <- res / sqrt(beta) # normalization
  res <- matrix(res, ncol = K) 
  est=res%*%a1
  errf=(m(pts)-est)^2
  result=sqrt(pracma::trapz(pts,errf))
  return(result)
}

batchrunL2=function(data_batch,acoe,ext=0,c0=0.5,domain=c(0,1),h=1/5,q0=5,C=1/2,m,Crho=0.01){
  ln=1000 #the number of discrete time points to predict density function
  EV=100  #the number of discrete time points to predict regression function
  tnew=seq(domain[1],domain[2],length=ln)
  tnew1=seq(domain[1],domain[2],length=EV)
  K=length(data_batch)  
  b=length(data_batch[[1]][[1]])
  n=K*b
  qstart=q0
  qend=floor(n^h/C)
  inc=qend-qstart
  Q=seq(qstart,qend,by=1)
  S=rep(0,inc)
  tauseq=rep(0,inc)
  B=rep(0,inc)
  tau=rep(0,inc)
  for( j in 1:inc)
  {
    S[j]= floor((C*Q[j+1])^(1/h))
    B[j]=ceiling(S[j]/b)
    tauseq[j]=floor(c0*S[j])
    tau[j]=ceiling(tauseq[j]/b)
  }
  
  lostn=rep(0,Q[length(Q)]) 
  nq=rep(0,Q[length(Q)]) 
  rhoopt=rep(0,K) 
  regest=matrix(0,K,EV) 
  
  lpreinc=length(which(tau==1))
  lpreB=length(which(B==1))
  preB=B
  if(lpreB>0){preB=B[-which(B==1)]}
  n=0
  tauseq=c(0,tauseq)
  
  l2esterr=rep(0,K)
  l2approerr=rep(0,K)
  maintime=c()
  querytime=c()
  #' Update
  {
    print(paste('K=',1))
    t=data_batch[[1]][,1]
    y=data_batch[[1]][,2]
    # maintain
    t0=Sys.time()
    q=qstart+lpreB
    if(lpreinc>0){q=qstart+lpreinc}
    n=n+length(t)
    G=computeG(t=t,y=y,q=q,domain=extend(domain,ext))
    G1=G[1:(qstart+lpreB)]
    nq[1:q]=n
    nq1=1/nq[1:(qstart+lpreB)]
    G2=diag(nq1)%*%G1 
    p=q
    theta=apply(evaluate.basis(K=p,domain=domain,grid=unlist(t)),2,sum)
    theta2=theta/n
    t1=Sys.time()
    # query
    # calculate H
    R1 <- list(ahat=theta2,ext=0,domain=domain)
    fhat=predict.regfunc(R1,tnew)*(domain[2]-domain[1])
    Phi=evaluate.basis(K=q,domain=extend(domain,ext),grid=tnew)
    H=t(Phi)%*%diag(fhat)%*%Phi/ln
    penmat=computeW(q=(qstart+lpreB),domain=extend(domain,ext))
    rhoopt[1]=Crho/(n^((7*h+1)/2))
    ahat= solve(H[1:(qstart+lpreB),1:(qstart+lpreB)]+rhoopt[1]*penmat,G2)
    aq=acoe[1:(qstart+lpreB)]
    l2esterr[1]=estl2(a1=ahat,a2=aq,K=(qstart+lpreB),domain=extend(domain,ext))
    l2approerr[1]=appl2(a1=aq,m=m,K=(qstart+lpreB),domain=extend(domain,ext))
    R <- list(ahat=ahat,penmat=penmat,ext=ext,rho=rhoopt[1],G=G,domain=domain)
    class(R) <- 'regfunc'
    regest[1,]=predict.regfunc(R,tnew1)
    t2=Sys.time()
    q=qstart+lpreB
    maintime=c(maintime,difftime(t1,t0,units = 'secs'))
    querytime=c(querytime,difftime(t2,t1,units = 'secs'))
  }
  if(2%in%preB)
  {
    q=q+length(which(preB==2))
    preB=preB[-which(preB==2)]
  }
  for (i in 2:(preB[1]-1)){
    print(paste('K=',i))
    t=data_batch[[i]][,1] 
    y=data_batch[[i]][,2]
    # maintain
    t0 <- Sys.time()
    if(i %in% tau)
    {
      G=c(G,rep(0,length(which(tau==i))))
      theta=c(theta,rep(0,length(which(tau==i))))
      lostn[length(G)]=n
    }
    n=n+length(t)
    nq[1:q]=n-lostn[1:q]
    nq1=1/nq[1:q]
    p=q
    theta1=apply(evaluate.basis(K=length(theta),domain=domain,grid=unlist(t)),2,sum)
    theta=theta+theta1
    theta2=diag(nq1)%*%(theta[1:p])
    G=G+computeG(t=t,y=y,q=length(G),domain=extend(domain,ext))
    G2=diag(nq1)%*%G[1:q]
    t1 <- Sys.time()
    # query
    R1 <- list(ahat=theta2,ext=0,domain=domain)
    fhat=predict.regfunc(R1,tnew)*(domain[2]-domain[1])
    Phi=evaluate.basis(K=length(G),domain=extend(domain,ext),grid=tnew)
    H=t(Phi)%*%diag(fhat)%*% Phi/ln
    penmat=computeW(q=q,domain=extend(domain,ext))
    rhoopt[i]=Crho/(n^((7*h+1)/2))
    ahat= solve(H[1:q,1:q]+rhoopt[i]*penmat,G2)
    aq=acoe[1:q]
    l2esterr[i]=estl2(a1=ahat,a2=aq,K=q,domain=extend(domain,ext))
    l2approerr[i]=appl2(a1=aq,m=m,K=q,domain=extend(domain,ext))
    R <- list(ahat=ahat,penmat=penmat,rho=rhoopt[i],ext=ext,G=G,domain=domain)
    class(R) <- 'regfunc'
    regest[i,]=predict.regfunc(R,tnew1)
    t2 <- Sys.time()
    maintime=c(maintime,difftime(t1,t0,units = 'secs'))
    querytime=c(querytime,difftime(t2,t1,units = 'secs'))
  }
  q=q+1
  if(length(which(preB==3))>1)
  {
    preB=preB[-1]
    q=q+length(which(preB==3))-1
  }
  for(j in 1:(length(preB)-1)){
    for(i in (preB[j]):(preB[j+1]-1)) {
      print(paste('K=',i))
      t=data_batch[[i]][,1]  
      y=data_batch[[i]][,2] 
      # maintain
      t0 <- Sys.time()
      if(i %in% tau)
      {
        G=c(G,rep(0,length(which(tau==i))))
        theta=c(theta,rep(0,length(which(tau==i))))
        lostn[length(G)]=n
      }
      n=n+length(t)
      nq[1:q]=n-lostn[1:q]
      nq1=1/nq[1:q]
      p=q
      theta1=apply(evaluate.basis(K=length(theta), domain=domain,grid=unlist(t)),2,sum)
      theta=theta+theta1
      theta2=diag(nq1)%*%(theta[1:p])
      G=G+computeG(t=t,y=y,q=length(G),domain=extend(domain,ext))
      G2=diag(nq1)%*%G[1:q]
      t1 <- Sys.time()
      # query
      R1 <- list(ahat=theta2,ext=0,domain=domain)
      fhat=predict.regfunc(R1,tnew)*(domain[2]-domain[1])
      Phi=evaluate.basis(K=length(G), domain=extend(domain,ext),grid=tnew)
      H=t( Phi)%*%diag(fhat)%*% Phi/ln
      penmat=computeW(q=q,domain=extend(domain,ext))
      rhoopt[i]=Crho/(n^((7*h+1)/2))
      ahat= solve(H[1:q,1:q]+rhoopt[i]*penmat,G2)
      aq=acoe[1:q]
      l2esterr[i]=estl2(a1=ahat,a2=aq,K=q,domain=extend(domain,ext))
      l2approerr[i]=appl2(a1=aq,m=m,K=q,domain=extend(domain,ext))
      R <- list(ahat=ahat,penmat=penmat,rho=rhoopt[i],ext=ext,G=G,domain=domain)
      class(R) <- 'regfunc'
      regest[i,]=predict.regfunc(R,tnew1)
      t2 <- Sys.time()
      maintime=c(maintime,difftime(t1,t0,units = 'secs'))
      querytime=c(querytime,difftime(t2,t1,units = 'secs'))
    }
    q=q+1
  }
  if(i<K){
    for(i in (B[length(B)]):K) {
      print(paste('K=',i))
      t=data_batch[[i]][,1]  
      y=data_batch[[i]][,2] 
      t0 <- Sys.time()
      n=n+length(t)
      nq[1:q]=n-lostn[1:q]
      nq1=1/nq[1:q]
      p=q
      theta1=apply(evaluate.basis(K=length(theta), domain=domain,grid=unlist(t)),2,sum)
      theta=theta+theta1
      theta2=diag(nq1)%*%(theta[1:p])
      G=G+computeG(t=t,y=y,q=length(G),domain=extend(domain,ext))
      G2=diag(nq1)%*%G[1:q]
      t1 <- Sys.time()
      # query
      R1 <- list(ahat=theta2,ext=0,domain=domain)
      fhat=predict.regfunc(R1,tnew)*(domain[2]-domain[1])
      Phi=evaluate.basis(K=length(G), domain=extend(domain,ext),grid=tnew)
      H=t(Phi)%*%diag(fhat)%*%Phi/ln
      penmat=computeW(q=q,domain=extend(domain,ext))
      rhoopt[i]=Crho/(n^((7*h+1)/2))
      ahat= solve(H[1:q,1:q]+rhoopt[i]*penmat,G2)
      aq=acoe[1:q]
      l2esterr[i]=estl2(a1=ahat,a2=aq,K=q,domain=extend(domain,ext))
      l2approerr[i]=appl2(a1=aq,m=m,K=q,domain=extend(domain,ext))
      R <- list(ahat=ahat,penmat=penmat,rho=rhoopt[i],ext=ext,G=G,domain=domain)
      class(R) <- 'regfunc'
      regest[i,]=predict.regfunc(R,tnew1)
      t2 <- Sys.time()
      maintime=c(maintime,difftime(t1,t0,units = 'secs'))
      querytime=c(querytime,difftime(t2,t1,units = 'secs'))
    }
  }
  RMISE=rep(0,K)
  mut=m(tnew1)
  for(i in 1:K){
    RMISE[i]=RMISEfunc(tnew1,regest[i,],mut)}
  R=list(RMISE,l2esterr,l2approerr)
  return(R)
}

############################### phase transition of uniform convergence rate  ################################
#' function to calculate the RMISE between true function and estimated function with uniform rate.
#' @param t a vector containing time points of observations to be evaluated.
#' @param f1 the values of true function at time points (a vector)
#' @param f2 the values of estimated function at time points (a vector)
#' @return the value of RMISE
Unifunc<-function(t,f1,f2)
{
  err=abs(f1-f2)
  result=max(err)
  return(result)
}

#'function to calculate the estimation error with sup-norm
#'@param a1 true coefficient vector
#'@param a2 estimated coefficient vector
#'@param K the basis number
#'@param domain the t domain
#'@return estimation error
estuni=function(a1,a2,K,domain)
{
  pts=seq(domain[1],domain[2],length=10000)
  a <- domain[1]
  b <- domain[2]
  x <- 0
  y <- 1
  alpha <- (b*x-a*y)/(x-y)
  beta <- (a-b)/(x-y)
  pts <- (pts - alpha) / beta
  res <- sapply(1:K, function(k)
    if (k == 1) rep(1, 10000)
    else if (k %% 2 == 0)  sqrt(2) * cos(k * pi * pts)
    else sqrt(2) * sin((k - 1) * pi * pts))
  res <- res / sqrt(beta) # normalization
  res <- matrix(res, ncol = K) 
  result=max(abs(res%*%(a1-a2)))
  return(result)
}

#'function to calculate the approximation error with sup-norm
#'@param a1 true coefficient vector
#'@param m true regression function
#'@param K the basis number
#'@param domain the t domain
#'@return approximation error
appuni=function(a1,m,K,domain)
{
  pts=seq(domain[1],domain[2],length=10000)
  a <- domain[1]
  b <- domain[2]
  x <- 0
  y <- 1
  alpha <- (b*x-a*y)/(x-y)
  beta <- (a-b)/(x-y)
  pts <- (pts - alpha) / beta
  res <- sapply(1:K, function(k)
    if (k == 1) rep(1, 10000)
    else if (k %% 2 == 0)  sqrt(2) * cos(k * pi * pts)
    else sqrt(2) * sin((k - 1) * pi * pts))
  res <- res / sqrt(beta) # normalization
  res <- matrix(res, ncol = K) 
  est=res%*%a1
  result=max(abs(m(pts)-est))
  return(result)
}

batchrunUNI=function(data_batch,acoe,ext=0,c0=0.5,domain=c(0,1),h=1/5,q0=5,C=0.5,m,Crho=0.01){
  ln=1000 #the number of discrete time points to predict density function
  EV=100  #the number of discrete time points to predict regression function
  tnew=seq(domain[1],domain[2],length=ln)
  tnew1=seq(domain[1],domain[2],length=EV)
  K=length(data_batch)  
  b=length(data_batch[[1]][[1]])
  n=K*b
  qstart=q0
  qend=floor(n^h/C)
  inc=qend-qstart
  Q=seq(qstart,qend,by=1)
  S=rep(0,inc)
  tauseq=rep(0,inc)
  B=rep(0,inc)
  tau=rep(0,inc)
  for( j in 1:inc)
  {
    S[j]= floor((C*Q[j+1])^(1/h))
    B[j]=ceiling(S[j]/b)
    tauseq[j]=floor(c0*S[j])
    tau[j]=ceiling(tauseq[j]/b)
  }
  
  lostn=rep(0,Q[length(Q)]) 
  nq=rep(0,Q[length(Q)]) 
  rhoopt=rep(0,K) 
  regest=matrix(0,K,EV) 
  
  lpreinc=length(which(tau==1))
  lpreB=length(which(B==1))
  preB=B
  if(lpreB>0){preB=B[-which(B==1)]}
  n=0
  tauseq=c(0,tauseq)
  
  uniesterr=rep(0,K)
  uniapproerr=rep(0,K)
  maintime=c()
  querytime=c()
  #' Update
  {
    print(paste('K=',1))
    t=data_batch[[1]][,1]
    y=data_batch[[1]][,2]
    # maintain
    t0=Sys.time()
    q=qstart+lpreB
    if(lpreinc>0){q=qstart+lpreinc}
    n=n+length(t)
    G=computeG(t=t,y=y,q=q,domain=extend(domain,ext))
    G1=G[1:(qstart+lpreB)]
    nq[1:q]=n
    nq1=1/nq[1:(qstart+lpreB)]
    G2=diag(nq1)%*%G1 
    p=q
    theta=apply(evaluate.basis(K=p,domain=domain,grid=unlist(t)),2,sum)
    theta2=theta/n
    t1=Sys.time()
    # query
    # calculate H
    R1 <- list(ahat=theta2,ext=0,domain=domain)
    fhat=predict.regfunc(R1,tnew)*(domain[2]-domain[1])
    Phi=evaluate.basis(K=q,domain=extend(domain,ext),grid=tnew)
    H=t(Phi)%*%diag(fhat)%*%Phi/ln
    penmat=computeW(q=(qstart+lpreB),domain=extend(domain,ext))
    rhoopt[1]=Crho/(n^((7*h+1)/2))
    ahat=solve(H[1:(qstart+lpreB),1:(qstart+lpreB)]+rhoopt[1]*penmat,G2)
    aq=acoe[1:(qstart+lpreB)]
    uniesterr[1]=estuni(a1=ahat,a2=aq,K=(qstart+lpreB),domain=extend(domain,ext))
    uniapproerr[1]=appuni(a1=aq,m=m,K=(qstart+lpreB),domain=extend(domain,ext))
    R <- list(ahat=ahat,penmat=penmat,ext=ext,rho=rhoopt[1],G=G,domain=domain)
    class(R) <- 'regfunc'
    regest[1,]=predict.regfunc(R,tnew1)
    t2=Sys.time()
    q=qstart+lpreB
    maintime=c(maintime,difftime(t1,t0,units = 'secs'))
    querytime=c(querytime,difftime(t2,t1,units = 'secs'))
  }
  if(2%in%preB)
  {
    q=q+length(which(preB==2))
    preB=preB[-which(preB==2)]
  }
  for (i in 2:(preB[1]-1)){
    print(paste('K=',i))
    t=data_batch[[i]][,1] 
    y=data_batch[[i]][,2]
    # maintain
    t0 <- Sys.time()
    if(i %in% tau)
    {
      G=c(G,rep(0,length(which(tau==i))))
      theta=c(theta,rep(0,length(which(tau==i))))
      lostn[length(G)]=n
    }
    n=n+length(t)
    nq[1:q]=n-lostn[1:q]
    nq1=1/nq[1:q]
    p=q
    theta1=apply(evaluate.basis(K=length(theta),domain=domain,grid=unlist(t)),2,sum)
    theta=theta+theta1
    theta2=diag(nq1)%*%(theta[1:p])
    G=G+computeG(t=t,y=y,q=length(G),domain=extend(domain,ext))
    G2=diag(nq1)%*%G[1:q]
    t1 <- Sys.time()
    # query
    R1 <- list(ahat=theta2,ext=0,domain=domain)
    fhat=predict.regfunc(R1,tnew)*(domain[2]-domain[1])
    Phi=evaluate.basis(K=length(G),domain=extend(domain,ext),grid=tnew)
    H=t(Phi)%*%diag(fhat)%*% Phi/ln
    penmat=computeW(q=q,domain=extend(domain,ext))
    rhoopt[i]=Crho/(n^((7*h+1)/2))
    ahat= solve(H[1:q,1:q]+rhoopt[i]*penmat,G2)
    aq=acoe[1:q]
    uniesterr[i]=estuni(a1=ahat,a2=aq,K=q,domain=extend(domain,ext))
    uniapproerr[i]=appuni(a1=aq,m=m,K=q,domain=extend(domain,ext))
    R <- list(ahat=ahat,penmat=penmat,rho=rhoopt[i],ext=ext,G=G,domain=domain)
    class(R) <- 'regfunc'
    regest[i,]=predict.regfunc(R,tnew1)
    t2 <- Sys.time()
    maintime=c(maintime,difftime(t1,t0,units = 'secs'))
    querytime=c(querytime,difftime(t2,t1,units = 'secs'))
  }
  q=q+1
  if(length(which(preB==3))>1)
  {
    preB=preB[-1]
    q=q+length(which(preB==3))-1
  }
  for(j in 1:(length(preB)-1)){
    for(i in (preB[j]):(preB[j+1]-1)) {
      print(paste('K=',i))
      t=data_batch[[i]][,1]  
      y=data_batch[[i]][,2] 
      # maintain
      t0 <- Sys.time()
      if(i %in% tau)
      {
        G=c(G,rep(0,length(which(tau==i))))
        theta=c(theta,rep(0,length(which(tau==i))))
        lostn[length(G)]=n
      }
      n=n+length(t)
      nq[1:q]=n-lostn[1:q]
      nq1=1/nq[1:q]
      p=q
      theta1=apply(evaluate.basis(K=length(theta), domain=domain,grid=unlist(t)),2,sum)
      theta=theta+theta1
      theta2=diag(nq1)%*%(theta[1:p])
      G=G+computeG(t=t,y=y,q=length(G),domain=extend(domain,ext))
      G2=diag(nq1)%*%G[1:q]
      t1 <- Sys.time()
      # query
      R1 <- list(ahat=theta2,ext=0,domain=domain)
      fhat=predict.regfunc(R1,tnew)*(domain[2]-domain[1])
      Phi=evaluate.basis(K=length(G), domain=extend(domain,ext),grid=tnew)
      H=t( Phi)%*%diag(fhat)%*% Phi/ln
      penmat=computeW(q=q,domain=extend(domain,ext))
      rhoopt[i]=Crho/(n^((7*h+1)/2))
      ahat= solve(H[1:q,1:q]+rhoopt[i]*penmat,G2)
      aq=acoe[1:q]
      uniesterr[i]=estuni(a1=ahat,a2=aq,K=q,domain=extend(domain,ext))
      uniapproerr[i]=appuni(a1=aq,m=m,K=q,domain=extend(domain,ext))
      R <- list(ahat=ahat,penmat=penmat,rho=rhoopt[i],ext=ext,G=G,domain=domain)
      class(R) <- 'regfunc'
      regest[i,]=predict.regfunc(R,tnew1)
      t2 <- Sys.time()
      maintime=c(maintime,difftime(t1,t0,units = 'secs'))
      querytime=c(querytime,difftime(t2,t1,units = 'secs'))
    }
    q=q+1
  }
  if(i<K){
    for(i in (B[length(B)]):K) {
      print(paste('K=',i))
      t=data_batch[[i]][,1]  
      y=data_batch[[i]][,2] 
      t0 <- Sys.time()
      n=n+length(t)
      nq[1:q]=n-lostn[1:q]
      nq1=1/nq[1:q]
      p=q
      theta1=apply(evaluate.basis(K=length(theta), domain=domain,grid=unlist(t)),2,sum)
      theta=theta+theta1
      theta2=diag(nq1)%*%(theta[1:p])
      G=G+computeG(t=t,y=y,q=length(G),domain=extend(domain,ext))
      G2=diag(nq1)%*%G[1:q]
      t1 <- Sys.time()
      # query
      R1 <- list(ahat=theta2,ext=0,domain=domain)
      fhat=predict.regfunc(R1,tnew)*(domain[2]-domain[1])
      Phi=evaluate.basis(K=length(G), domain=extend(domain,ext),grid=tnew)
      H=t(Phi)%*%diag(fhat)%*%Phi/ln
      penmat=computeW(q=q,domain=extend(domain,ext))
      rhoopt[i]=Crho/(n^((7*h+1)/2))
      ahat= solve(H[1:q,1:q]+rhoopt[i]*penmat,G2)
      aq=acoe[1:q]
      uniesterr[i]=estuni(a1=ahat,a2=aq,K=q,domain=extend(domain,ext))
      uniapproerr[i]=appuni(a1=aq,m=m,K=q,domain=extend(domain,ext))
      R <- list(ahat=ahat,penmat=penmat,rho=rhoopt[i],ext=ext,G=G,domain=domain)
      class(R) <- 'regfunc'
      regest[i,]=predict.regfunc(R,tnew1)
      t2 <- Sys.time()
      maintime=c(maintime,difftime(t1,t0,units = 'secs'))
      querytime=c(querytime,difftime(t2,t1,units = 'secs'))
    }
  }
  UNISE=rep(0,K)
  mut=m(tnew1)
  for(i in 1:K){
    UNISE[i]=Unifunc(tnew1,regest[i,],mut)}
  R=list(UNISE,uniesterr,uniapproerr)
  return(R)
}