setwd('/home/yangy/OGAM_code_and_data')
library(latex2exp)
######################### generate table 1 ################################
### predict
rm(list=ls())
load('datasets/flight/2000_gam.Rdata')
Xtest <- cbind(df$CRSDepTime,df$HistDelayRate)
ytest <- df$Delayed
tab <- c()
rm(df)

load('res/flight/flight_full_online.Rdata')
pre<-c()
M <- 25 # No evalpoints
eval_vec <- seq(0.05, 0.95, length.out = M)
sub_streams <- c(1,seq(10,50,10),seq(100,1000,100), seq(1500,2500,500),3283)

for(K in sub_streams){
  sp1 <- smooth.spline(eval_vec, beta_store[,2*K-1])
  sp2 <- smooth.spline(eval_vec, beta_store[,2*K])
  beta1 <- predict(sp1, Xtest[,1])$y
  beta2 <- predict(sp2, Xtest[,2])$y
  p1 <- exp(beta0_store[K]+beta1+beta2)/(1+exp(beta0_store[K]+beta1+beta2))
  yhat <- p1>0.5
  pre<-c(pre,sum((yhat-ytest)==0)/length(ytest))
}
tab <- 1-pre

load('res/flight/flight_full_batch.Rdata')
pre<-c()
M <- 25 # No evalpoints
eval_vec <- seq(0.05, 0.95, length.out = M)

for(K in 1:length(sub_streams)){
  sp1 <- smooth.spline(eval_vec, beta_store[,2*K-1])
  sp2 <- smooth.spline(eval_vec, beta_store[,2*K])
  beta1 <- predict(sp1, Xtest[,1])$y
  beta2 <- predict(sp2, Xtest[,2])$y
  p1 <- exp(beta0_store[K]+beta1+beta2)/(1+exp(beta0_store[K]+beta1+beta2))
  yhat <- p1>0.5
  pre<-c(pre,sum((yhat-ytest)==0)/length(ytest))
}
tab<-rbind(tab,1-pre)
round(tab[,c(1,2,4,6,7,11,16,18,10)],4)

######################### generate Figure XXX--XXX !!! #############################
rm(list=ls())
sub_streams <- c(1,seq(10,50,10),seq(100,1000,100), seq(1500,2500,500),3283)
load('res/flight/flight_full_online.Rdata')
beta_o <- beta_store; to <- time
load('res/flight/flight_full_batch.Rdata')
beta_b <- beta_store; tb <- time

# component functions
{
  
  Ks<-c(4,7,16,20)
  Kos<-c(30,100,1000,3283)
  ylim1 <- range(cbind(beta_b[,2*Ks-1],beta_o[,2*Kos-1]))
  ylim2 <- range(cbind(beta_b[,2*Ks],beta_o[,2*Kos]))
  x<-seq(0.05, 0.95,length.out = 25)
  
  pdf('fig/flights_betas.pdf', height = 6 , width= 12)
  par(mai=c(0.5,0.5,0.4,0.4),omi=c(0.5,0.5,0.5,0),mfrow=c(2,4))
  plot(x,beta_b[,2*Ks[1]-1],type='l',lty=2,xaxt='n',yaxt='n', 
       xlab='',  ylab = '',lwd=1.6, cex.main=2, ylim=ylim1,
       font.main = 1,family='serif')
  lines(x,beta_o[,2*Kos[1]-1])
  mtext(TeX("$\\beta_{1}$"), side=2, line=3.5, adj=0.5, cex=1.6, family='serif',las=1)
  mtext(paste('K=',sub_streams[Ks[1]],sep=''), side=3, line=0.2, adj=0.5, cex=1.6,  family='serif')
  axis(side=2,at=round(seq(min(ylim1),max(ylim1),length.out = 3),3),
       labels = round(seq(min(ylim1),max(ylim1),length.out = 3),3),
       font.axis = 1, cex.axis=1.6, family='serif')
  for(i in 2:4){
    plot(x,beta_b[,2*Ks[i]-1],type='l',lty=2,xaxt='n',yaxt='n', 
         xlab='',  ylab = '',lwd=1.6, cex.main=2, ylim=ylim1,
         font.main = 1,family='serif')
    lines(x,beta_o[,2*Kos[i]-1])
    mtext(paste('K=',sub_streams[Ks[i]],sep=''), side=3, line=0.2, adj=0.5, cex=1.6,  family='serif')
  }
  
  plot(x,beta_b[,2*Ks[1]],type='l',lty=2,xaxt='n',yaxt='n', 
       xlab='',  ylab = '',lwd=1.6, cex.main=2, ylim=ylim2,
       font.main = 1,family='serif')
  lines(x,beta_o[,2*Kos[1]])
  mtext(TeX("$\\beta_{2}$"), side=2, line=3.5, adj=0.5, cex=1.6, family='serif',las=1)
  mtext('x', side=1, line=3.5, adj=0.5, cex=1.6,  family='serif')
  axis(side=1,at=round(seq(0,1,length.out = 5),2),font.axis = 1,
       labels =round(seq(0,1,length.out = 5),2),cex.axis=1.6, family='serif')
  axis(side=2,at=round(seq(min(ylim2),max(ylim2),length.out = 3),3),
       labels = round(seq(min(ylim2),max(ylim2),length.out = 3),3),
       font.axis = 1, cex.axis=1.6, family='serif')
  for(i in 2:4){
    plot(x,beta_b[,2*Ks[i]],type='l',lty=2,xaxt='n',yaxt='n', 
         xlab='',  ylab = '',lwd=1.6, cex.main=2, ylim=ylim2,
         font.main = 1,family='serif')
    lines(x,beta_o[,2*Kos[i]])
    axis(side=1,at=round(seq(0,1,length.out = 5),2),font.axis = 1,
         labels =round(seq(0,1,length.out = 5),2),cex.axis=1.6, family='serif')
    mtext('x', side=1, line=3.5, adj=0.5, cex=1.6,  family='serif')
  }
  dev.off()
}

# time
{
  pdf('fig/flights_time.pdf', height = 5, width= 7.5)
  par(mai=c(0.5,0.6,0.3,0.3),omi=c(0.5,0.5,0,0),mfrow=c(1,1))
  
  idx=1:4
  plot(sub_streams[idx],tb[idx],
       ylim=c(0,max(tb[idx])),xlim=range(sub_streams[idx]),
       type='l',lwd=1.5,lty=2,
       xlab='',ylab = '',xaxt='n',font.main=1,
       cex.lab=2.5,cex.axis=2,cex.main=2.5,family='serif')
  lines(1:max(sub_streams[idx]),to[1:max(sub_streams[idx])])
  mtext('time(secs)', side=2, line=3.5, adj=0.5, cex=2, family='serif')
  mtext('blocks', side=1, line=3.5, adj=0.5, cex=2, family='serif')
  axis(side=1,at=c(1,seq(10,40,10)), family='serif',
       labels = c(1,seq(10,40,10)),cex.axis=2)
  dev.off()
}









