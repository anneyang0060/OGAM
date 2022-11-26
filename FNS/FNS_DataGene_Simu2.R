beta1_fun <-function(x){-4*(x-0.5)^2}
beta2_fun <- function(x){2*x^3}
beta3_fun <- function(x){sin(pi*x)}
beta4_fun <- function(x){cos(pi*x)}

beta1_fun_deri2 <- function(x){-8}
beta2_fun_deri2 <- function(x){12*x}
beta3_fun_deri2 <- function(x){-pi^2*cos(pi*x)}
beta4_fun_deri2 <- function(x){-pi^2*sin(pi*x)}

gene_data <- function(n){
  
  x <- matrix(runif(n*4),n,4)
  e <- rnorm(n)
  y <- (beta1_fun(x[,1]) + beta2_fun(x[,2]) + beta3_fun(x[,3])
        + beta4_fun(x[,4]) + e)
  
  data <- list(x,y)
  names(data) <- c('x','y')
  return(data)
  
}