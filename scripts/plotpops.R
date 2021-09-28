plotpop <- function(d, jd, i=1, name=popnames){
  # jd=jags.dat
  
  x <- jd$t_ti[i,1]:jd$t_ti[i,jd$nt[i]]
  y <- yup <- ylow <- c()
  for(t in x){
    y <- c(y, mean( d[, paste('N[',i,',',t,']',sep='')] ) )
    ylow <- c(ylow, quantile( d[, paste('N[',i,',',t,']',sep='')], probs=c(0.025) ) )
    yup <- c(yup, quantile( d[, paste('N[',i,',',t,']',sep='')], probs=c(0.975) ) )
    Yi <- apply(jd$Yj[i,,], 1, sum, na.rm=T)
  }
  
  if(length(y) == 1){
    plot(y~x, xlab='t', ylab='N', ylim=c(min(Yi, na.rm=T), max(yup)), main=name[i], pch=16)
    points(x=x, y=yup, pch=1)
    points(x=x, y=ylow, pch=1)
    points(x=jd$t_ti[i,], y=Yi, col='red', pch=16)
  } else{
    plot(y~x, type='l', xlab='t', ylab='N', ylim=c(min(Yi, na.rm=T), max(yup)), main=name[i])
    lines(x=x, y=yup, lty=2)
    lines(x=x, y=ylow, lty=2)
    points(x=jd$t_ti[i,], y=Yi)  
  }
  
}
