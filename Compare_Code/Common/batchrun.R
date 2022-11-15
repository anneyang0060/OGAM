#' function to estimate regression function with streaming data (in batch form)
#' @param data_batch the data in the batch form
#' @param ext the value of extension margin
#' @param c0 the value of c_{circle} defined in paper
#' @param domain the t domain
#' @param h parameter of function S
#' @param q0 the original basis number
#' @param C coefficient of function S
#' @param m true regression function; when analyzing the real data, m=NULL.
#' @param Crho coefficient of tuning parameter rho
#' @param rep  rep=5 for weekday data; rep=2 for weekend data; and rep=1 for other settings. 
#' @param seq  when results in particular batches are required (e.g. n/3,2n/3 and n), seq is a vector including the index of batches. Default value is NULL.
#' @return a list containing RMISE, time of maintaining and querying
batchrun=function(data_batch,ext=0,c0=0.5,domain=c(0,1),h=1/5,q0=5,C=1/2,m=NULL,Crho=1,rep=1,seq=NULL){
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
  rhoopt=rep(0,K) # storing optimal rho sequence
  regest=list() # storing the estimation results
  maintime=c()
  querytime=c()
  
  lpreinc=length(which(tau==1))
  lpreB=length(which(B==1))
  preB=B
  if(lpreB>0){preB=B[-which(B==1)]}
  n=0
  tauseq=c(0,tauseq)
  
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
    penmat=computeW(q=q,domain=extend(domain,ext))
    rhoopt[1]=Crho*rep/(n^((7*h+1)/2))
    ahat= solve(H[1:(qstart+lpreB),1:(qstart+lpreB)]+rhoopt[1]*penmat[1:(qstart+lpreB),1:(qstart+lpreB)],G2)
    R <- list(ahat=ahat,penmat=penmat,ext=ext,rho=rhoopt[1],G=G,domain=domain)
    class(R) <- 'regfunc'
    if(!is.null(m)){regest=list(predict.regfunc(R,tnew1))}
    else
    {
      if(1 %in% seq)
      {
        est=predict.regfunc(R,tnew1)
        regest=c(regest,list(est))
      }  
    }
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
  for (i in 2:(preB[1]-1))
    {
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
    rhoopt[i]=Crho*rep/(n^((7*h+1)/2))
    ahat= solve(H[1:q,1:q]+rhoopt[i]*penmat,G2)
    R <- list(ahat=ahat,penmat=penmat,rho=rhoopt[i],ext=ext,G=G,domain=domain)
    class(R) <- 'regfunc'
    if(!is.null(m)){regest=c(regest,list(predict.regfunc(R,tnew1)))}
    else
    {
      if(i %in% seq)
      {
        est=predict.regfunc(R,tnew1)
        regest=c(regest,list(est))
      }  
    }
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
      rhoopt[i]=Crho*rep/(n^((7*h+1)/2))
      ahat= solve(H[1:q,1:q]+rhoopt[i]*penmat,G2)
      R <- list(ahat=ahat,penmat=penmat,rho=rhoopt[i],ext=ext,G=G,domain=domain)
      class(R) <- 'regfunc'
      if(!is.null(m)){regest=c(regest,list(predict.regfunc(R,tnew1)))}
      else
      {
        if(i %in% seq)
        {
          est=predict.regfunc(R,tnew1)
          regest=c(regest,list(est))
        }  
      }
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
      rhoopt[i]=Crho*rep/(n^((7*h+1)/2))
      ahat= solve(H[1:q,1:q]+rhoopt[i]*penmat,G2)
      R <- list(ahat=ahat,penmat=penmat,rho=rhoopt[i],ext=ext,G=G,domain=domain)
      class(R) <- 'regfunc'
      if(!is.null(m)){regest=c(regest,list(predict.regfunc(R,tnew1)))}
      else
      {
        if(i %in% seq)
        {
          est=predict.regfunc(R,tnew1)
          regest=c(regest,list(est))
        }  
      }
      t2 <- Sys.time()
      maintime=c(maintime,difftime(t1,t0,units = 'secs'))
      querytime=c(querytime,difftime(t2,t1,units = 'secs'))
    }
  }
  if(!is.null(m))
  {
    RMISE=rep(0,K)
    mut=m(tnew1)
    for(i in 1:K){RMISE[i]=RMISEfunc(tnew1,regest[[i]],mut)}
    R=list(RMISE,maintime,querytime)  
  }
  else
  {
    est=predict.regfunc(R,tnew1)
    regest=c(regest,list(est))
    R=regest
  }
  return(R)
}