setwd('/home/yangy/OGAM_code_and_data')
library(latex2exp)

############################## organize data  ###########
d <- 4; Kmax <- 600
sub_streams <- c(1,seq(20,Kmax,20))
Nsim <- 100
L <- c(3,5,10)
lb_th <- 1/(1+0.1831/L+0.0032/L^2)

rssb<-matrix(0,d,length(sub_streams)); tb <- rep(0,length(sub_streams))
rss_L3<-matrix(0,d,Kmax); to <- matrix(0,length(L),Kmax)
rss_L5<-matrix(0,d,Kmax)
rss_L10<-matrix(0,d,Kmax)
for(i in 1:100){
  load(paste0('res/sim2/batch_gam_',i,'.Rdata'))
  rssb <- rssb+rss; tb <- tb+time
  load(paste0('res/sim2/online_gam_L3_',i,'.Rdata'))
  rss_L3 <- rss_L3+rss; to[1,] <- to[1,]+time
  load(paste0('res/sim2/online_gam_L5_',i,'.Rdata'))
  rss_L5 <- rss_L5+rss; to[2,] <- to[2,]+time
  load(paste0('res/sim2/online_gam_L10_',i,'.Rdata'))
  rss_L10 <- rss_L10+rss; to[3,] <- to[3,]+time
}
rsso <- array(0,dim = c(d,Kmax,3))
rsso[,,1] <- rss_L3; rsso[,,2] <- rss_L5; rsso[,,3] <- rss_L10
tb <- tb/Nsim; to <- to/Nsim
rm(rss_L3, rss_L5, rss_L10)
# rssb <- colSums(rssb); rsso <- colSums(rsso)

############################## plot figure XX-XX #################
# efficiency (figure XX)
{
  L <- c(3,5,10)
  lb_th <- 1/(1+0.1831/L+0.0032/L^2)
  titles <- c('(a)','(b)','(c)')
  
  pdf('fig/sim2_empirical_eff.pdf', height = 4.8 , width= 10)
  par(mai=c(0.4,0.4,0.5,0.5), omi=c(0.5,0.5,0,0),mfrow=c(2,3))
  for(i in 1:3){
    lb <- 0.85
    plot(sub_streams,rssb[1,]/rsso[1,sub_streams,i],type='l',
         ylim=c(lb,1.05),xaxt='n',yaxt='n',
         xlab='',  ylab = '',lwd=1.6, cex.main=1.7,
         font.main = 1,main=titles[i], family='serif')
    lines(sub_streams,rssb[2,]/rsso[2,sub_streams,i],lty=2,lwd=1.6)
    lines(sub_streams,rssb[3,]/rsso[3,sub_streams,i],lty=3,lwd=1.6)
    lines(sub_streams,rssb[4,]/rsso[4,sub_streams,i],lty=4,lwd=1.6)
    lines(sub_streams, rep(lb_th[i],length(sub_streams)),lty=3,lwd=3)
    axis(side=2,at=round(seq(lb,1,length.out = 3),3), font.axis = 1,
         labels =round(seq(lb,1,length.out = 3),3),cex.axis=1.6, family='serif')
    axis(side=1,at=c(1,seq(250,1000,250)),font.axis = 1,
         labels = c(1,seq(250,1000,250)),cex.axis=1.6, family='serif')
    if(i==1){
      mtext('efficiency', side=2, line=3.5, adj=0.5, cex=1.2, family='serif')
    }
  }
  rssb1 <- colSums(rssb); rsso1 <- colSums(rsso)
  titles <- c('(e)','(f)','(g)')
  for(i in 1:3){
    lb <- 0.85
    plot(sub_streams,rssb1/rsso1[sub_streams,i],type='l',
         ylim=c(lb,1.05),xaxt='n',yaxt='n',
         xlab='',  ylab = '',lwd=1.6, cex.main=1.7,
         font.main = 1,main=titles[i], family='serif')
    lines(sub_streams, rep(lb_th[i],length(sub_streams)),lty=3,lwd=3)
    axis(side=2,at=round(seq(lb,1,length.out = 3),3), font.axis = 1,
         labels =round(seq(lb,1,length.out = 3),3),cex.axis=1.6, family='serif')
    axis(side=1,at=c(1,seq(250,1000,250)),font.axis = 1,
         labels = c(1,seq(250,1000,250)),cex.axis=1.6, family='serif')
    mtext('blocks', side=1, line=3.5, adj=0.5, cex=1.2,  family='serif')
    if(i==1){
      mtext('efficiency', side=2, line=3.5, adj=0.5, cex=1.2, family='serif')
    }
  }
  dev.off()
}

# time (figure XX)
{
  pdf('fig/sim2_simutime.pdf', height = 4.5, width= 6)
  par(mai=c(0.5,0.6,0.3,0.3),omi=c(0.5,0.5,0,0),mfrow=c(1,1))
  
  idx=2:11; to <- t(to)
  plot(sub_streams[idx],tb[idx],
       ylim=c(0,max(tb[idx])),xlim=range(sub_streams[idx]),
       type='l',lwd=1.5,lty=2,
       xlab='',ylab = '',xaxt='n',font.main=1,
       cex.lab=2.5,cex.axis=1.6,cex.main=2.5,family='serif')
  for(i in 1:3){lines(sub_streams[idx],to[idx,i])}
  mtext('time(secs)', side=2, line=3.5, adj=0.5, cex=1.8, family='serif')
  mtext('blocks', side=1, line=3.5, adj=0.5, cex=1.8, family='serif')
  axis(side=1,at=c(20,100,200), family='serif',
       labels = c(20,100,200),cex.axis=1.6)
  
  dev.off()
}






















