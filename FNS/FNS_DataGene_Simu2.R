beta1_fun <-function(x){-4*(x-0.5)^2}
beta2_fun <- function(x){x^3}
beta3_fun <- function(x){sin(pi/2*x)}
beta4_fun <- function(x){cos(2*pi*x)}
beta5_fun <- function(x){sin(2*pi*x^2)}


gene_data <- function(n){
  
  x <- matrix(runif(n*4),n,4)
  e <- rnorm(n)
  y <- (beta1_fun(x[,1]) + beta2_fun(x[,2]) + beta3_fun(x[,3])
        + beta4_fun(x[,4]) + e)
    
  data <- list(x,y)
  names(data) <- c('x','y')
  return(data)
  
}