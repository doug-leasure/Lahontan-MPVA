plotsim <- function(dat, popname, extinct.prob, outdir, datyears=NA, Nmax=1e9){
  # origwd <- getwd()
  # dir.create(outdir)
  # setwd(outdir)
  
  jpeg(paste(outdir, popname, '.jpg', sep=''))
  
  # Data
  yrs <- 1:ncol(dat)
  x.yrs <- gsub('X', '', names(dat))
  
  # Setup Plot
  layout(matrix(c(1,2,3), nrow=3, ncol=1), heights=c(0.45, 0.25, 0.30))
  
  # Top panel
  c <- 0.001
  
  par(mar=c(1,5.5,4,1))
  plot(x=x.yrs, y=dat['Nmean',yrs]+c, type='l', ylim=c(0, min(Nmax, max(dat['Nup',yrs], na.rm=T))), 
       xlab=NA, xaxt='n', ylab='Population Size (N)',#\nlog-scale', 
       main=paste('Population: ',popname,'\n Extinction Probability: ',extinct.prob, sep=''),
       cex.lab=1.1, cex.axis=1, cex.main=2)#, log='y')
  
  lines(x=x.yrs, y=dat['Nmean',yrs]+c, type='l', lty=1, lwd=1.5)
  lines(x=x.yrs, y=dat['Nup',yrs]+c, type='l', lty=2, lwd=1)
  lines(x=x.yrs, y=dat['Nlow',yrs]+c, type='l', lty=2, lwd=1)
  axis(side=1, labels=NA)
  
  #points(x=datyears, y=1)
  rug(datyears, lwd=3)
  
  # Bottom panels
  par(mar=c(1,5.5,1,1))  
  plot(x=x.yrs, y=dat['rmean',yrs], type='l', col='blue', lty=1, 
       ylim=c(min(dat['rlow',yrs],na.rm=T), max(dat['rup',yrs],na.rm=T)), ylab='Population\nGrowth Rate', xlab=NA, 
       xaxt='n', cex.lab=1.1, cex.axis=1, lwd=2)
  abline(h=0)
  lines(x=x.yrs, y=dat['rup',yrs], type='l', col='blue', lty=2)
  lines(x=x.yrs, y=dat['rlow',yrs], type='l', col='blue', lty=2)
  axis(side=1, labels=NA)
  
  par(mar=c(4.5,5.5,1,1))  
  plot(x=x.yrs, y=dat['phimean',yrs], type='l', col='red', lty=1, 
       ylab='Effect of Density', xlab='Year', ylim=c(min(dat['philow',yrs],na.rm=T), 0), 
       cex.lab=1.1, cex.axis=1, lwd=2)
  abline(h=0)
  lines(x=x.yrs, y=dat['phiup',yrs], type='l', col='red', lty=2)
  lines(x=x.yrs, y=dat['philow',yrs], type='l', col='red', lty=2)
  
  dev.off()
  
  #setwd(origwd)
}