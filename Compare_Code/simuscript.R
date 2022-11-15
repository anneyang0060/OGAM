setwd('/home/yangy/gam/Compare_Code')
source('Common/helpers.R')
source('Common/OPLS.R')
source('Common/batchrun.R')
source('Common/datagene.R')


################################## define and load helper functions ################################
#' Routine to compare different one-pass algorithms
#' @param m true regression function
#' @param n sample size 
#' @param snr the signal-to-noise ratio 
#' @param domain the t domain
#' @param ext the value of extension margin
#' @param q0 the original basis number
#' @param n0 the number of observations to determine the coefficient of tuning parameter
#' @param C coefficient to determine the basis number q
simrun2=function(m,n,snr,domain,ext,q0,n0,C){
  data=datagene(m=m,snr=snr,n=n,domain=domain)
  # data to batch
  t0 <- Sys.time()
  t=data[[1]][1:n0]
  y=data[[2]][1:n0]
  Cseq=10^seq(-4,0,length.out=20)
  hseq=c(1/3,1/4,1/5,1/6,1/7)
  result=tune(t=t,y=y,Cseq=Cseq,C=C,hseq=hseq,ext=ext,domain=domain)
  h=result[[1]]
  Crho=result[[2]]
  t1 <- Sys.time()
  time=difftime(t1,t0,units = 'secs')
  batchsize=100
  databatch=datatobatch(data,batchsize)
  databatch1=datatobatch1(data,batchsize)
  
  basisresult=batchrun(data_batch=databatch,ext=ext,c0=0.5,domain=domain,h=h,q0=q0,C=C,m=m,Crho=Crho)
  asympresult=batchrunshrink(data_batch=databatch,domain=domain,m=m)
  localresult=batchrunlocal(data_batch=databatch,domain=domain,m=m)
  fullresult=fullrun(data_batch=databatch,domain=domain,m=m)
  
  # R=list(asympresult)
  R=list(basisresult,asympresult,localresult,fullresult,time)
  return(R)
}


################ define true regression function and parameters #########################
basisnum=8
acoe=rep(0,basisnum)
for(i in 1:basisnum)
{
  acoe[i]=i^{-1.5}
}
f=function(s,k){
  if(k==1) return(1)
  else if(k%%2==0) return(sqrt(2)*cos(k*pi*s))
  else return(sqrt(2)*sin((k-1)*pi*s))
}
m1=function(s){Reduce("+",lapply(1:basisnum,function(i){acoe[i]*f(s,i)}))}

n=100000
domain=c(0,1)
snr=2
n0=1000
C=1/2
q0=5
ext=c(0,0.1,0)
Nsim <- 100
sds <- 1:Nsim
  
#################  perform simulations for comparing different estimators ################

# Note: the simulations of DPSR was performed by the same synthetic data and by 
#       the Python version developed by and obtained through personal communications with
#       the authors of DPSR; the Python code thus cannot be made publicly here.
Mcl<-50
cl<-makeCluster(Mcl)
registerDoParallel(cl)
comp<- foreach(sd=sds) %dopar%{
  set.seed(sd)
  m1resultcom=simrun2(m=m1,n=n,snr=snr,domain=domain,ext=ext[1],C=C,q0=q0,n0=n0)
  save(m1resultcom, file = paste0('res/compare_',sd,'.Rdata'))
}
stopCluster(cl)

