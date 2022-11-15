setwd('/home/yangy/OGAM_code_and_data')
library(latex2exp)

############################## organize data  ###########
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

############################## Figure XX-XX #################
# efficiency (figure XX)
{
  L <- c(3,5,10)
  lb_th <- 1/(1+0.1831/L+0.0032/L^2)
  titles <- c('(a)','(b)','(c)')
  
  idx=1:51
  pdf('fig/sim1_empirical_eff.pdf', height = 2.4 , width= 10)
  par(mai=c(0.4,0.4,0.5,0.5), omi=c(0.5,0.5,0,0),mfrow=c(1,3))
  for(i in 1:3){
    lb <- 0.85
    plot(sub_streams[idx],rssb[idx,1]/rsso[sub_streams[idx],1,i],type='l',
         ylim=c(lb,1.05),xaxt='n',yaxt='n',
         xlab='',  ylab = '',lwd=1.6, cex.main=1.7,
         font.main = 1,main=titles[i], family='serif')
    lines(sub_streams[idx],rssb[idx,2]/rsso[sub_streams[idx],2,i],lty=2,lwd=1.6)
    lines(sub_streams[idx], rep(lb_th[i],length(idx)),lty=3,lwd=3)
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
  pdf('fig/sim1_simutime.pdf', height = 4.5, width= 6)
  par(mai=c(0.5,0.6,0.3,0.3),omi=c(0.5,0.5,0,0),mfrow=c(1,1))
  
  idx=2:16
  plot(sub_streams[idx],tb[idx],
       ylim=c(0,max(tb[idx])),xlim=range(sub_streams[idx]),
       type='l',lwd=1.5,lty=2,
       xlab='',ylab = '',xaxt='n',font.main=1,
       cex.lab=2.5,cex.axis=1.6,cex.main=2.5,family='serif')
  for(i in 1:3){lines(sub_streams[idx],to[idx,i])}
  mtext('time(secs)', side=2, line=3.5, adj=0.5, cex=1.8, family='serif')
  mtext('blocks', side=1, line=3.5, adj=0.5, cex=1.8, family='serif')
  axis(side=1,at=c(20,100,200,300), family='serif',
       labels = c(20,100,200,300),cex.axis=1.6)
  
  dev.off()
}

# bandwidth (figure XX)
ho <- ho[,,3]
{
  
  band_plot <- function(p, comp, h1, h2){
    
    Kmax <- 1000
    d <- 2
    pd <- 1
    R <-100
    set.seed(2020)
    sds <- ceiling(runif(R)*1e6)
    n <- ceiling(rnorm(Kmax,100,10))
    n[1] <- 2000
    N <- cumsum(n)
    
    
    L <- 10
    K<-1000
    eta<-matrix(0,L,K)
    eta[1,]<-h1
    for(k in 1:K){
      eta[2:L,k]<-sapply(2:L, function(l) ((L-l+1)/L)^(1/p)*eta[1,k])
    }
    
    for(K in c(200,500,1000)){
      m_eta<-eta[,1]
      k<-2
      ord<-sapply(1:L,function(l) which.min(abs(eta[l,k]-m_eta)))
      m_eta<-(N[k-1]*m_eta[ord]+n[k]*eta[,k])/N[k]
      index<-cbind(ord,1:L)
      for(k in 3:K){
        ord<-sapply(1:L,function(l) which.min(abs(eta[l,k]-m_eta)))
        m_eta<-(N[k-1]*m_eta[ord]+n[k]*eta[,k])/N[k]
        index1<-matrix(0,L,k)
        index1[,k]<-1:L
        index1[,k-1]<-ord
        for(l in 1:L){
          index1[l,1:(k-2)]<-index[ord[l],1:(k-2)]
        }
        index<-index1
      }
      rm(index1)
      pseu.index<-index[1,]
      pseu<-sapply(1:K,function(k) eta[pseu.index[k],k]) 
      plot(rep(1,L),eta[,1],pch=16,cex=0.2,xlim=c(1,1000),ylim=range(eta),
           xlab='',ylab = '',family='serif', xaxt='n',yaxt='n',
           cex.lab=2.5,cex.axis=2,cex.main=2.5,font.main=1)
      for(k in seq(10,K,length.out = (K/10))){
        points(rep(k,L),eta[,k],pch=16,cex=0.2)
      }
      points(pseu,pch=16,cex=0.6)
      lines(eta[1,1:K])
      if(comp==1){title(paste('K=',K,sep=''))}
      if(comp==2){
        axis(side=1,at=c(1,seq(250,1000,250)),font.axis = 1,
             labels = c(1,seq(250,1000,250)),cex.axis=1.6, family='serif')
        mtext('blocks', side=1, line=3.5, adj=0.5, cex=1.6,  family='serif')
      }
      if(K==200){
        mtext(TeX("$\\eta_{l1}^{(k)}$"), side=2, line=3.5, adj=0.5, cex=1.6, 
              family='serif',las=1)
        axis(side=2,at=round(seq(min(hb[,comp]),max(hb[,comp]),length.out = 3),3),
             font.axis = 1,
             labels =round(seq(min(hb[,comp]),max(hb[,comp]),length.out = 3),3),
             cex.axis=1.6, family='serif')
      }
      
      if(comp == 2){
        axis(side=1,at=c(1,seq(250,1000,250)),font.axis = 1,
             labels = c(1,seq(250,1000,250)),cex.axis=1.6, family='serif')
        mtext('blocks', side=1, line=3.5, adj=0.5, cex=1.6,  family='serif')
      }
    }
    
  }
  
  
  pdf('fig/sim1_bandwidths.pdf', height = 5 , width= 11)
  par(mai=c(0.5,0.5,0.4,0.4),omi=c(0.5,0.5,0,0),mfrow=c(2,4))
  plot(sub_streams,hb[,1],type='l',lty=2,
       xaxt='n',yaxt='n', 
       xlab='',  ylab = '',lwd=1.6, cex.main=2,
       font.main = 1,family='serif')
  lines(sub_streams,ho[sub_streams,1])
  mtext(TeX("$h_{K1}$"), side=2, line=3.5, adj=0.5, cex=1.6, family='serif',las=1)
  axis(side=2,at=round(seq(min(hb[,1]),max(hb[,1]),length.out = 3),3),
       font.axis = 1,
       labels =round(seq(min(hb[,1]),max(hb[,1]),length.out = 3),3),
       cex.axis=1.6, family='serif')
  band_plot(p = 5,comp=1, ho[,1], hb[,1])
  
  plot(sub_streams,hb[,2],type='l',lty=2,
       xaxt='n',yaxt='n', 
       xlab='',  ylab = '',lwd=1.6, cex.main=2,
       font.main = 1,family='serif')
  lines(sub_streams,ho[sub_streams,2])
  mtext(TeX("$h_{K2}$"), side=2, line=3.5, adj=0.5, cex=1.6, family='serif',las=1)
  axis(side=2,at=round(seq(min(hb[,2]),max(hb[,2]),length.out = 3),3),
       font.axis = 1,
       labels =round(seq(min(hb[,2]),max(hb[,2]),length.out = 3),3),
       cex.axis=1.6, family='serif')
  axis(side=1,at=c(1,seq(250,1000,250)),font.axis = 1,
       labels = c(1,seq(250,1000,250)),cex.axis=1.6, family='serif')
  mtext('blocks', side=1, line=3.5, adj=0.5, cex=1.6,  family='serif')
  band_plot(p = 5,comp=2, ho[,2], hb[,2])
  
  dev.off()
  
}



############### Figure XXXX ##########################
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


