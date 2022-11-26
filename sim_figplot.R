setwd('.../OGAM')
library(latex2exp)


############### Figure 4 ##########################
{
  lb_mu <- function(L,eps){-eps/(1-eps)*L^2+0.224*L+0.0017}
  lb_gam <- function(L,eps){-eps/(1-eps)*L^2+0.5924*L-0.3901}
  eps1 <- c(seq(0.16,0.03,-0.01),seq(0.027,0.01,-0.003),seq(0.01,0.005,-0.001),seq(0.0045,0.0025,-0.0001))
  L_mu <- sapply(eps1,function(x){uniroot(lb_mu,c(1,1000),eps=x)$root})
  # plot(L_mu, log(eps1))
  eps2 <- c(seq(0.16,0.03,-0.01),seq(0.027,0.01,-0.003),seq(0.01,0.005,-0.001),seq(0.0045,0.0025,-0.0001))
  L_gam <- sapply(eps2,function(x){uniroot(lb_gam,c(1,1000),eps=x)$root})
  # plot(L_gam, log(eps2))
  L <- 3:40
  lb_mu <- 1/(1+0.224/L+0.0017/L^2)
  lb_gam <- 1/(1+0.5924/L-0.3901/L^2)
  pdf('figures/theoretical_eff2.pdf',
      height = 4.32 , width= 11.52)
  par(mai=c(1,0.5,0.35,0.5),omi=c(0,0.5,0,0),mfrow=c(1,2))
  plot(L,lb_mu,xlab='',ylab=' ',ylim=c(0.85,1),
       type='l',lwd=1.5,family='serif',
       cex.lab=1.35,cex.axis=1.35)
  mtext('efficiency', side=2, line=3, adj=0.5,
        family='serif',cex=1.35)
  mtext('L', side=1, line=3, adj=0.5,
        family='serif',cex=1.35,font=3)
  lines(L,lb_gam,lty=2)
  plot(log(eps1),L_mu, xlab='',ylab=' ', 
       ylim=range(c(L_mu,L_gam)),xlim=log(c(0.0025,0.2)),
       type='l',lwd=1.5,family='serif',
       cex.lab=1.35,cex.axis=1.35)
  lines(log(eps2),L_gam,lty=2)
  mtext('L', side=2, line=3, adj=0.5,
        family='serif',cex=1.35,font=3)
  mtext(TeX('$\\log(\\delta)$'), side=1, line=3, adj=0.5,
        family='serif',cex=1.35)
  dev.off()
}


############### Figure 5 ##########################
ck <- function(L,d,r){
  0.183/L+0.001/L^2+L*(d+2)^3*r
}
memo <- function(c,d){
  0.183*10^3*4*(d+2)^3/c^2/2^30
}

pdf('fig/selection_L.pdf', height = 4 , width= 11)
par(mai=c(0.55,0.55,0.5,0.8), omi=c(0.5,0.5,0,0), mfrow=c(1,2))

Ls = 3:30
{
  ds = (7:2)
  ll = 1 # memory=1G
  for(d in ds){
    y = ck(Ls,d,3.73*1e-6)
    if(d==max(ds)){
      plot(Ls, y, type='l', ylim=c(0.01,0.09), xaxt='n',yaxt='n', 
           xlab='',  ylab = '',lty=ll)
    }else{
      lines(Ls, y, lty=ll)
    }
    points(Ls[which.min(y)], min(y), pch=4)
    ll = ll + 1
  }
  mtext('total cost', side=2, line=3.5, adj=0.5, cex=1.6, family='serif')
  mtext('L', side=1, line=3.5, adj=0.5, cex=1.6,  family='serif')
  axis(side=2,at=seq(0.01,0.09,0.02), font.axis = 1,
       labels=seq(0.01,0.09,0.02),cex.axis=1.5, family='serif')
  axis(side=1,at=c(3,seq(5,30,5)),font.axis = 1,
       labels=c(3,seq(5,30,5)),cex.axis=1.5, family='serif')
}

ds = 2:8
{
  y1 = memo(0.005, ds)
  y2 = memo(0.0075, ds)
  y3 = memo(0.01, ds)
  y4 = memo(0.02, ds)
  plot(ds, y1, type='l', ylim=c(0,28),#range(c(y1,y2,y3,y4)),
       xaxt='n',yaxt='n', 
       xlab='',  ylab = '',lty=1)
  lines(ds,y2,lty=2)
  lines(ds,y3,lty=3)
  lines(ds,y4,lty=4)
  mtext('required memory (GB)', side=2, line=3.5, adj=0.5, cex=1.6, family='serif')
  mtext('d', side=1, line=3.5, adj=0.5, cex=1.6,  family='serif')
  axis(side=2,at=seq(0,28,4), font.axis = 1,
       labels=seq(0,28,4),cex.axis=1.5, family='serif')
  axis(side=1,at=2:8,font.axis = 1,
       labels=2:8,cex.axis=1.5, family='serif')
}

dev.off()

############################## organize data for Sim1  ###########
Kmax <- 1000
sub_streams <- c(1,seq(20,Kmax,20))
Nsim <- 100
d <- 2

### batch
hb <- array(0, dim = c(length(sub_streams), d, Nsim))
rssb <- hb; tb <- matrix(0, length(sub_streams), Nsim)

for(i in 1:Nsim){
  
  load(paste('res/sim1/batch_gam_',i,'.Rdata',sep=''))
  hb[,,i] <- t(band); rssb[,,i] <- t(rss); tb[,i] <- time
  
}

hb <- rowMeans(hb, dims = 2)
rssb <- rowMeans(rssb, dims = 2)
tb <- rowMeans(tb)

### online 
Ls <- c(3,5,10)
ho <- array(0, dim = c(Kmax, d, length(Ls)))
rsso <- ho; to <- matrix(0, Kmax, length(Ls))

for(L in Ls){
  
  ho1 <- array(0, dim = c(Kmax, d, Nsim))
  rsso1 <- ho1; to1 <- matrix(0, Kmax, Nsim)
  
  for(i in 1:Nsim){
    
    load(paste('res/sim1/online_gam_L',L,'_',i,'.Rdata',sep=''))
    ho1[,,i] <- t(band); rsso1[,,i] <- t(rss); to1[,i] <- time
    
  }
  
  ho[,, which(Ls==L)] <- rowMeans(ho1, dims = 2)
  rsso[,, which(Ls==L)] <- rowMeans(rsso1, dims = 2)
  to[,which(Ls==L)] <- rowMeans(to1)
  
}

rm(ho1, rsso1, to1, band, time, rss)





############################## organize data for Sim2 ########
d <- 4; Kmax <- 600
sub_streams <- c(1,seq(20,Kmax,20))
Nsim <- 100
L <- c(3,5,10)
lb_th <- 1/(1+0.1831/L+0.0032/L^2)

rssb1 <-matrix(0,d,length(sub_streams)); tb1 <- rep(0,length(sub_streams))
rss_L3<-matrix(0,d,Kmax); to1 <- matrix(0,length(L),Kmax)
rss_L5<-matrix(0,d,Kmax)
rss_L10<-matrix(0,d,Kmax)
for(i in 1:100){
  load(paste0('res/sim2/batch_gam_',i,'.Rdata'))
  rssb1 <- rssb1+rss; tb1 <- tb1+time
  load(paste0('res/sim2/online_gam_L3_',i,'.Rdata'))
  rss_L3 <- rss_L3+rss; to1[1,] <- to1[1,]+time
  load(paste0('res/sim2/online_gam_L5_',i,'.Rdata'))
  rss_L5 <- rss_L5+rss; to1[2,] <- to1[2,]+time
  load(paste0('res/sim2/online_gam_L10_',i,'.Rdata'))
  rss_L10 <- rss_L10+rss; to1[3,] <- to1[3,]+time
}
rsso1 <- array(0,dim = c(d,Kmax,3))
rsso1[,,1] <- rss_L3; rsso1[,,2] <- rss_L5; rsso1[,,3] <- rss_L10
tb1 <- tb1/Nsim; to1 <- to1/Nsim
rm(rss_L3, rss_L5, rss_L10)


############################## Figure 6 #################
{
  L <- c(3,5,10)
  lb_th <- 1/(1+0.1831/L+0.0032/L^2)
  
  pdf('fig/sim_empirical_eff.pdf', height = 5 , width= 10)
  par(mai=c(0.4,0.4,0.5,0.5), omi=c(0.5,0.5,0,0),mfrow=c(2,3))
  
  sub_streams <- c(1,seq(20,1000,20))
  idx=1:51
  titles <- c('(a)','(b)','(c)')
  for(i in 1:3){
    lb <- 0.85
    plot(sub_streams[idx],rssb[idx,1]/rsso[sub_streams[idx],1,i],type='l',
         ylim=c(lb,1.05),xaxt='n',yaxt='n',col='red',
         xlab='',  ylab = '',lwd=1.6, cex.main=1.7,
         font.main = 1,main=titles[i], family='serif')
    lines(sub_streams[idx],rssb[idx,2]/rsso[sub_streams[idx],2,i],col='blue',lty=2,lwd=1.6)
    lines(sub_streams[idx], rep(lb_th[i],length(idx)),lty=2,lwd=1.6)
    axis(side=2,at=round(seq(lb,1,length.out = 3),3), font.axis = 1,
         labels =round(seq(lb,1,length.out = 3),3),cex.axis=1.6, family='serif')
    axis(side=1,at=c(1,seq(250,1000,250)),font.axis = 1,
         labels = c(1,seq(250,1000,250)),cex.axis=1.6, family='serif')
    mtext('blocks', side=1, line=2.8, adj=0.5, cex=1.2,  family='serif')
    if(i==1){
      mtext('efficiency', side=2, line=3.5, adj=0.5, cex=1.2, family='serif')
    }
  }
 
  sub_streams <- c(1,seq(20,600,20))
  titles <- c('(d)','(e)','(f)')
  for(i in 1:3){
    lb <- 0.85
    plot(sub_streams,rssb1[1,]/rsso1[1,sub_streams,i],type='l',
         ylim=c(lb,1.05),xaxt='n',yaxt='n',col='red',
         xlab='',  ylab = '',lwd=1.6, cex.main=1.7,
         font.main = 1,main=titles[i], family='serif')
    lines(sub_streams,rssb1[2,]/rsso1[2,sub_streams,i],col='blue',lty=2,lwd=1.6)
    lines(sub_streams,rssb1[3,]/rsso1[3,sub_streams,i],col='green',lty=3,lwd=1.6)
    lines(sub_streams,rssb1[4,]/rsso1[4,sub_streams,i],col='orange',lty=4,lwd=1.6)
    lines(sub_streams, rep(lb_th[i],length(sub_streams)),lty=2,lwd=1.6)
    axis(side=2,at=round(seq(lb,1,length.out = 3),3), font.axis = 1,
         labels =round(seq(lb,1,length.out = 3),3),cex.axis=1.6, family='serif')
    axis(side=1,at=c(1,seq(200, 600, 200)),font.axis = 1,
         labels = c(1,seq(200, 600, 200)),cex.axis=1.6, family='serif')
    mtext('blocks', side=1, line=2.8, adj=0.5, cex=1.2,  family='serif')
    if(i==1){
      mtext('efficiency', side=2, line=3.5, adj=0.5, cex=1.2, family='serif')
    }
  }
  
  dev.off()
}

############################## Figure 7 #################
{
  pdf('fig/sim_simutime.pdf', height = 4.5, width= 12)
  par(mai=c(0.5,0.6,0.5,0.4),omi=c(0.5,0.5,0,0),mfrow=c(1,2))
  
  idx=2:16
  sub_streams <- c(1,seq(20,1000,20))
  plot(sub_streams[idx],tb[idx],
       ylim=c(0,max(tb[idx])),xlim=range(sub_streams[idx]),
       type='l',lwd=1.5,lty=2,cex.main=1.7,
       font.main = 1,main='(a)',
       xlab='',ylab = '',xaxt='n',font.main=1,
       cex.lab=2.5,cex.axis=1.6,family='serif')
  for(i in 1:3){lines(sub_streams[idx],to[idx,i])}
  mtext('time(secs)', side=2, line=3.5, adj=0.5, cex=1.8, family='serif')
  mtext('blocks', side=1, line=3.5, adj=0.5, cex=1.8, family='serif')
  axis(side=1,at=c(20,100,200,300), family='serif',
       labels = c(20,100,200,300),cex.axis=1.6)
  
  idx=2:11; to1 <- t(to1)
  sub_streams <- c(1,seq(20,600,20))
  plot(sub_streams[idx],tb1[idx],
       ylim=c(0,max(tb1[idx])),xlim=range(sub_streams[idx]),
       type='l',lwd=1.5,lty=2,cex.main=1.7,
       font.main = 1,main='(b)',
       xlab='',ylab = '',xaxt='n',font.main=1,
       cex.lab=2.5,cex.axis=1.6,family='serif')
  for(i in 1:3){lines(sub_streams[idx],to1[idx,i])}
  mtext('blocks', side=1, line=3.5, adj=0.5, cex=1.8, family='serif')
  axis(side=1,at=c(20,100,200), family='serif',
       labels = c(20,100,200),cex.axis=1.6)
  
  
  dev.off()
}

