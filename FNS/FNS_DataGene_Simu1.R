norm1 <- 0.227; norm2 <- 0.317
# h1: 0.52; h2: 0.35
beta1_fun <- function(x){cos(sqrt(2)*pi*x)-norm1}
beta2_fun <- function(x){x/2+sin(2*pi*x)/2-norm2}
c_fun <- 2+norm1+norm2

beta1_fun_deri <- function(x){-sqrt(2)*pi*sin(sqrt(2)*pi*x)}
beta2_fun_deri <- function(x){1/2+pi*cos(2*pi*x)}
beta1_fun_deri2 <- function(x){-2*pi^2*cos(sqrt(2)*pi*x)}
beta2_fun_deri2 <- function(x){-2*pi^2*sin(2*pi*x)}

gene_data <- function(n){
  
  mean <- c(0,0)
  sigma <- matrix(c(10, 0.9, 0.9, 10),nrow=2,ncol=2)
  i<-0
  x<-c()
  while(i<n){
    x1 <- mvrnorm(10*n,mean,sigma) 
    x1 <- x1[abs(x1[,1]) < 1 & abs(x1[,2]) < 1,]
    i <- i + nrow(x1)
    x <- rbind(x,x1)
  }
  x <- abs(x[1:n,])
  
  mu <- exp(beta1_fun(x[,1]) + beta2_fun(x[,2]) + c_fun)
  # mu <- exp(sin(2*pi*x[,1]) + 3 * x[,2]^5 -2 * x[,2])
  y <- 0
  for(i in 1:n){
    y[i]<-rpois(1,mu[i])
  }
  # y <- mu
  data <- list(x,y)
  names(data) <- c('x','y')
  return(data)
  
}