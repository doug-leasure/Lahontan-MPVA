# rm(list=ls())
# cat("\014") 
gc()

library(rjags);library(runjags);library(coda);library(lattice)

#setwd('D:/RESEARCH/2015 NASA Trout Project/wd/1.002')

load('output/objects/jags.dat.R')
attach(jags.dat, warn.conflicts=F)

thin <- 100

if (thin==1) { 
  d.thin <- d 
} else if (!file.exists(paste('output/objects/d.thin.',thin,'.R',sep=''))){
  d.thin <- d[seq(1, nrow(d), thin),]
  save(d.thin, file=paste('output/objects/d.thin.',thin,'.R',sep=''))
  gc()
} else {
  load('output/objects/d.thin.100.R')
}

# Expected y values from model
if(file.exists('output/objects/modpred.R')){
  load('output/objects/modpred.R')
  attach(modpred)
} else{
  modpred <- list()
  modpred$ypred <- modpred$exp0  <- modpred$yhat <- modpred$exp1 <- array(NA, dim=c(nrow(d.thin), npops, max(nt), max(nj, na.rm=T), max(npasses, na.rm=T)))
  modpred$Yjpred <- modpred$Yjhat <- modpred$exp2 <- array(NA, dim=c(nrow(d.thin), npops, max(nt), max(nj, na.rm=T)))
  modpred$Yihat <- modpred$exp3 <- array(NA, dim=c(nrow(d.thin), npops, max(nt)))
  for (i in 1:npops){
    print(paste('Pop:', i))
    for (t in 1:nt[i]){
      print(paste('  Year:', t))
      N <- d.thin[,paste('N[',i,',',t_ti[i,t],']',sep='')]
      Psite <- d.thin[,paste('Psite[',i,',',t,']',sep='')]
      
      modpred$exp3[,i,t] <- N * Psite
      modpred$Yihat[,i,t] <- rbinom(length(N), N, Psite)
      
      for(s in 1:nrow(d.thin)){
        PIsite <- t(d.thin[s, paste('PIsite[',i,',',t,',',1:nj[i,t],']',sep='')])
        
        modpred$exp2[s,i,t,1:nj[i,t]] <- Yi[i,t] * PIsite
        modpred$Yjhat[s,i,t,1:nj[i,t]] <- rmultinom(1, Yi[i,t], PIsite)
        modpred$Yjpred[s,i,t,1:nj[i,t]] <- rmultinom(1, modpred$Yihat[s,i,t], PIsite)
      }
      
      for (j in 1:nj[i,t]){
        print(paste('    Site:',j))
        
        for (s in 1:nrow(d.thin)){
          PI <- t(d.thin[s, paste('PI[',i,',',t,',',j,',',1:npasses[i,t,j],']',sep='')])
          
          modpred$exp1[s,i,t,j,1:npasses[i,t,j]] <- Yj[i,t,j] * PI
          modpred$yhat[s,i,t,j,1:npasses[i,t,j]] <- rmultinom(1, Yj[i,t,j], PI)
          
          modpred$exp0[s,i,t,j,1:npasses[i,t,j]] <- mean(modpred$Yjpred[,i,t,j]) * PI
          modpred$ypred[s,i,t,j,1:npasses[i,t,j]] <- rmultinom(1, modpred$Yjpred[s,i,t,j], PI)
        }
      }  
    }
  }
  save(modpred, file='output/objects/modpred.R')
  attach(modpred)
}


# Assess residuals
res <- bias <- prec <- acc <- c()
for (s in 1:nrow(d.thin)){
  resi <- ypred[s,,,,] - y
  bias[s] <- mean(resi, na.rm=T)
  prec[s] <- sd(resi, na.rm=T)
  acc[s] <- mean(abs(resi), na.rm=T)
  res <- c(res, resi)
}

pdf('output/plots/residuals.pdf')
boxplot(res, main='Residuals', ylab='(Observed - Expected) Fish caught per pass', outline=F, range=0)
boxplot(bias, main='Bias', ylab='Mean Residuals', outline=F, range=0)
boxplot(prec, main='Imprecision', ylab='St. Dev. of Residuals', outline=F, range=0)
boxplot(acc, main ='Inaccuracy', ylab='Mean Absolute Residuals', outline=F, range=0)
dev.off()


# Test quantities
ss <- function(obs, exp){sum((obs-exp)^2, na.rm=T)}
chi <- function(obs, exp){sum((obs-exp)^2/exp, na.rm=T)}
mad <- function(obs, exp){median(abs(obs-exp), na.rm=T)}
rmse <- function(obs, exp){sqrt(mean((obs-exp)^2, na.rm=T))}

Stat <- ss

# Calculate test quantities for observed data and for simulated data
T0 <- T0hat <- T1 <- T1hat <- T2 <- T2hat <- T3 <- T3hat <- c()
for (s in 1:nrow(d.thin)){
  T0[s] <- Stat(obs=y, exp=exp0[s,,,,])
  T0hat[s] <- Stat(obs=ypred[s,,,,], exp=exp0[s,,,,])
  T1[s] <- Stat(obs=y, exp=exp1[s,,,,])
  T1hat[s] <- Stat(obs=yhat[s,,,,], exp=exp1[s,,,,])
  T2[s] <- Stat(obs=Yj, exp=exp2[s,,,])
  T2hat[s] <- Stat(obs=Yjhat[s,,,], exp=exp2[s,,,])
  T3[s] <- Stat(obs=Yi, exp=exp3[s,,])
  T3hat[s] <- Stat(obs=Yihat[s,,], exp=exp3[s,,])
}

# Bayesian p-values
pvalue0 <- mean(T0hat > T0)
pvalue1 <- mean(T1hat > T1)
pvalue2 <- mean(T2hat > T2)
pvalue3 <- mean(T3hat > T3)

pvalue0
pvalue1
pvalue2
pvalue3



# Observed vs. predicted plot
jpeg('output/plots/fitplot.jpg')
plot(NA, xlim=c(min(y, na.rm=T), max(y, na.rm=T)), ylim=c(min(ypred, na.rm=T), max(ypred, na.rm=T)), ylab='Predicted ypred', xlab='Observed y', main=paste('Accuracy:', round(mean(acc),3)))
abline(0, 1, col='red')
for (i in 1:npops){
  for (t in 1:nt[i]){
    for (j in 1:nj[i,t]){
      for(m in 1:npasses[i,t,j]){
        points(x = y[i,t,j,m], y=mean(ypred[,i,t,j,m], na.rm=T))
        arrows(x0=y[i,t,j,m], y0=quantile(ypred[,i,t,j,m], probs=c(0.025), na.rm=T), y1=quantile(ypred[,i,t,j,m], probs=c(0.975), na.rm=T), length=0)
      }
    }
  }
}
dev.off()

jpeg('output/plots/fitplot1.jpg')
plot(NA, xlim=c(min(y, na.rm=T), max(y, na.rm=T)), ylim=c(min(yhat, na.rm=T), max(yhat, na.rm=T)), ylab='Predicted yhat', xlab='Observed y', main=paste('P =', round(pvalue1,3)))
abline(0, 1, col='red')
for (i in 1:npops){
  for (t in 1:nt[i]){
    for (j in 1:nj[i,t]){
      for(m in 1:npasses[i,t,j]){
        points(x = y[i,t,j,m], y=mean(yhat[,i,t,j,m], na.rm=T))
        arrows(x0=y[i,t,j,m], y0=quantile(yhat[,i,t,j,m], probs=c(0.025), na.rm=T), y1=quantile(yhat[,i,t,j,m], probs=c(0.975), na.rm=T), length=0)
      }
    }
  }
}
dev.off()

jpeg('output/plots/fitplot2.jpg')
plot(NA, xlim=c(min(Yj, na.rm=T), max(Yj, na.rm=T)), ylim=c(min(Yjhat, na.rm=T), max(Yjhat, na.rm=T)), ylab='Predicted Yjhat', xlab='Observed Yj', main=paste('P =', round(pvalue2,3)))
abline(0, 1, col='red')
for (i in 1:npops){
  for (t in 1:nt[i]){
    for (j in 1:nj[i,t]){
      points(x = Yj[i,t,j], y=mean(Yjhat[,i,t,j], na.rm=T))
      arrows(x0=Yj[i,t,j], y0=quantile(Yjhat[,i,t,j], probs=c(0.025), na.rm=T), y1=quantile(Yjhat[,i,t,j], probs=c(0.975), na.rm=T), length=0)
    }
  }
}
dev.off()

jpeg('output/plots/fitplot3.jpg')
plot(NA, xlim=c(min(Yi, na.rm=T), max(Yi, na.rm=T)), ylim=c(min(Yihat, na.rm=T), max(Yihat, na.rm=T)), ylab='Predicted Yihat', xlab='Observed Yi', main=paste('P =', round(pvalue3,3)))
abline(0, 1, col='red')
for (i in 1:npops){
  for (t in 1:nt[i]){
    points(x = Yi[i,t], y=mean(Yihat[,i,t], na.rm=T))
    arrows(x0=Yi[i,t], y0=quantile(Yihat[,i,t], probs=c(0.025), na.rm=T), y1=quantile(Yihat[,i,t], probs=c(0.975), na.rm=T), length=0)
  }
}
dev.off()

# Cleanup
rm(list=ls()[-which(ls() %in% c('d','d.thin','zm'))])
gc()
