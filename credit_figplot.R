setwd('/home/yangy/OGAM_code_and_data')
library(latex2exp)

d <- 4
m <- 10 # No evalpoints
eval_vec <- seq(0.1, 0.9, length.out = m)
Kmax <- 822
sub_streams <- c(1, seq(20,180,20), seq(200,700,100), 822)

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
ytest <- data$badloan[(1:Nfull)[-idx]]
Xtest <- as.matrix(data[(1:Nfull)[-idx],c('loan_amnt','int_rate','dti','annual_inc')])
rm(idx, data)


######################### generate table 2 ################################
### predict
load('res/credit/credit_online.Rdata')
pre <- c()

for(K in sub_streams){ 
  
  beta_est <- c()
  
  for(i in d:1){
    
    sp <- smooth.spline(eval_vec, beta_store[,d*K-(i-1)])
    beta_est <- cbind(beta_est, predict(sp, Xtest[,(d-i+1)])$y)
    
  }
  
  p1 <- exp(beta0_store[K]+rowSums(beta_est))/(1+exp(beta0_store[K]+rowSums(beta_est)))
  yhat <- p1>0.5
  pre <- c(pre, sum((yhat-ytest)==0)/length(ytest))
  
}
tab <- 1- pre 

load('res/credit/credit_batch.Rdata')
pre<-c()

for(K in 1:length(sub_streams)){
 
  beta_est <- c()
  
  for(i in d:1){
    
    sp <- smooth.spline(eval_vec, beta_store[,d*K-(i-1)])
    beta_est <- cbind(beta_est, predict(sp, Xtest[,(d-i+1)])$y)
    
  }
  
  p1 <- exp(beta0_store[K]+rowSums(beta_est))/(1+exp(beta0_store[K]+rowSums(beta_est)))
  yhat <- p1>0.5
  pre <- c(pre, sum((yhat-ytest)==0)/length(ytest))
  
}
tab <- rbind(tab,1-pre)
round(tab[,c(2,4,6,11,12,14,17)],4)


################


######################### generate Figure XXX--XXX !!! #############################
load('res/credit/credit_online.Rdata')
beta_o <- beta_store; to <- time
load('res/credit/credit_batch.Rdata')
beta_b <- beta_store; tb <- time

# component functions
{
  
  Ks<-c(6,17)
  Kos<-c(100,822)
  
  pdf('fig/credit_betas.pdf', height = 6 , width= 12)
  
  par(mai=c(0.5,0.5,0.4,0.4),omi=c(0.5,0.5,0.5,0),mfrow=c(2,4))
  x<-seq(0.05, 0.95,length.out = 25)
  eva<-seq(0.1, 0.9,length.out = 10)
  nc <- ncol(beta_o)
  beta_o1 <- c()
  for(i in 1:nc){
    sp <- smooth.spline(eva, beta_o[,i])
    beta_o1 <- cbind(beta_o1,predict(sp, x)$y)
  }
  nc <- ncol(beta_b)
  beta_b1 <- c()
  for(i in 1:nc){
    sp <- smooth.spline(eva, beta_b[,i])
    beta_b1 <- cbind(beta_b1,predict(sp, x)$y)
  }
  
  ylims <- c()
  for(i in d:1){
    ylims <- rbind(ylims, range(cbind(beta_b1[,d*Ks-(i-1)],beta_o1[,d*Kos-(i-1)])))
  }
  
  for(k in 2:1){
    for(i in d:1){
      plot(x,beta_b1[,d*Ks[k]-(i-1)],type='l',lty=2,xaxt='n',yaxt='n', 
           xlab='',  ylab = '',lwd=1.6, cex.main=2, ylim=ylims[d-i+1,],
           font.main = 1,family='serif')
      lines(x,beta_o1[,d*Kos[k]-(i-1)])
      axis(side=2,at=round(seq(min(ylims[d-i+1,]),max(ylims[d-i+1,]),length.out = 3),3),
           labels = round(seq(min(ylims[d-i+1,]),max(ylims[d-i+1,]),length.out = 3),3),
           font.axis = 1, cex.axis=1.6, family='serif')
      if(k==1){axis(side=1,at=round(seq(0,1,length.out = 5),2),font.axis = 1,
                      labels =round(seq(0,1,length.out = 5),2),cex.axis=1.6, family='serif')
        mtext(TeX("$\\x_{1}$"), side=1, line=3.5, adj=0.5, cex=1.6,  family='serif')}
      if(k==2){mtext(TeX("$\\beta_{1}$"), side=3, line=0.2, adj=0.5, cex=1.6, family='serif',las=1)}
      if(i==d){mtext(paste0('K=',Kos[3-k]), side=2, line=3.5, adj=0.5, cex=1.6,  family='serif')}
    }
  }
  dev.off()
}

# time
{
  pdf('fig/credit_time.pdf', height = 5, width= 7.5)
  par(mai=c(0.5,0.6,0.3,0.3),omi=c(0.5,0.5,0,0),mfrow=c(1,1))
  
  idx=1:5
  plot(sub_streams[idx],tb[idx],
       ylim=c(0,max(tb[idx])),xlim=range(sub_streams[idx]),
       type='l',lwd=1.5,lty=2,
       xlab='',ylab = '',xaxt='n',font.main=1,
       cex.lab=2.5,cex.axis=2,cex.main=2.5,family='serif')
  lines(1:max(sub_streams[idx]),to[1:max(sub_streams[idx])])
  mtext('time(secs)', side=2, line=3.5, adj=0.5, cex=2, family='serif')
  mtext('blocks', side=1, line=3.5, adj=0.5, cex=2, family='serif')
  axis(side=1,at=c(1,seq(20,80,20)), family='serif',
       labels = c(1,seq(20,80,20)),cex.axis=2)
  dev.off()
}









