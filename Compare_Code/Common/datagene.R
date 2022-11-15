#' function to generate standard deviation of noise via signal-to-noise ratio (snr)
#' @param data data from the regression function without noise
#' @param snr the value of snr
#' @return standard deviation of the noise
noisesd<-function(data,snr)
{
  datarms=mean(data)^2+var(data)
  sd=sqrt(datarms/snr)
  return(sd)
}

#' function to generate raw data
#' @param m true regression function
#' @param snr the value of snr
#' @param n sample size
#' @param domain the t domain
#' @return a list of data
datagene=function(m,snr,n,domain)
{
  t=runif(n,domain[1],domain[2])
  yclear=m(t)
  sdvalue=noisesd(yclear,snr)
  y=m(t)+rnorm(n,mean=0,sd=sdvalue)
  data=list(t,y)
  return(data)
}


#' function to transform raw data into batch form
#' @param data raw data in the form of data.frame
#' @param batchsize the number of observations in each batch
#' @return a list of data
datatobatch=function(data,batchsize)
{
  b=batchsize
  n=length(data[[1]])
  K=ceiling(n/b)
  data_batch=NULL
  for(i in 1:(K-1)) {
    batch_start<-(i-1)*b+1
    batch_end<-i*b
    data_batch[[i]] <- data.frame(data[[1]][batch_start: batch_end],data[[2]][batch_start: batch_end])
  }
  batch_start=(K-1)*b+1
  batch_end=n
  data_batch[[K]] <- data.frame(data[[1]][batch_start: batch_end],data[[2]][batch_start: batch_end])
  return(data_batch)
}

#' function to transform raw data into batch form (for non-streaming setting)
#' @param data raw data in the form of data.frame
#' @param batchsize the number of observations in each batch
#' @return a list of data
datatobatch1=function(data,batchsize)
{
  b=batchsize
  n=length(data[[1]])
  K=ceiling(n/b)
  data_batch=NULL
  for(i in 1:(K-1)) {
    batch_start<-1
    batch_end<-i*b
    data_batch[[i]] <- data.frame(data[[1]][batch_start: batch_end],data[[2]][batch_start: batch_end])
  }
  batch_start=1
  batch_end=n
  data_batch[[K]] <- data.frame(data[[1]][batch_start: batch_end],data[[2]][batch_start: batch_end])
  return(data_batch)
}