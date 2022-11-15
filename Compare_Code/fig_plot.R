setwd('/home/yangy/gam/Compare_Code')

######################### figure XX
onlinebasis <- 0
onlineshrink <- 0
onlinekernel <- 0
fullkernel <- 0
Nsim <- 100
for(i in 1:Nsim){
  
  load(paste0('Compare_Code/res/compare_',i,'.Rdata'))
  onlinebasis=onlinebasis+m1resultcom[[1]][[1]]/Nsim
  onlineshrink=onlineshrink+m1resultcom[[2]][[1]]/Nsim
  onlinekernel=onlinekernel+m1resultcom[[3]][[1]]/Nsim
  fullkernel=fullkernel+m1resultcom[[4]][[1]]/Nsim
  
}

# Note: the simulations of DPSR was performed by the same synthetic data and by 
#       the Python version developed by and obtained through personal communications with
#       the authors of DPSR; the Python code thus cannot be made publicly here.
load("Compare_Code/DPSR_results_from_python/DPSR_results.Rdata")
pdf('fig/compare.pdf', height = 4.5, width= 6)
par(mai=c(0.5,0.6,0.3,0.3),omi=c(0.5,0.5,0,0),mfrow=c(1,1))
plot(fullkernel/onlinekernel,col='red',type='l',lty=1,lwd=1.5,ylim=c(0.3,1.05),
     xlab='',ylab = '',xaxt='n',font.main=1,
     cex.lab=2.5,cex.axis=1.6,cex.main=2.5,family='serif')
lines(fullkernel/onlineshrink,col='green',lwd=1.5,lty=3)
lines(fullkernel/onlinebasis,col='blue',lwd=1.5,lty=2)
lines(fullkernel/rowMeans(mise),col='orange',lwd=1.5,lty=4)
mtext('efficiency', side=2, line=3.5, adj=0.5, cex=1.8, family='serif')
mtext('blocks', side=1, line=3.5, adj=0.5, cex=1.8, family='serif')
axis(side=1,at=c(1,seq(200,1000,200)), family='serif',
     labels =c(1,seq(200,1000,200)),cex.axis=1.6)
dev.off()


