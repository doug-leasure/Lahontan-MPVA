rm(list=ls())
gc()
cat("\014") 

# Working directory
setwd('D:/RESEARCH/2015 NASA Trout Project/wd/MPVA')

# CODA file
load('output/objects/d.yhat.R')

d <- d[,grepl('hat[',names(d),fixed=T)]

# JAGS data
load('output/objects/jags.dat.R')

# Residuals
resid <- data.frame(row.names=1:nrow(d))

for(i in 1:jags.dat$npops){
  for(t in 1:jags.dat$nt[i]){
    # resid[,paste('Yierr[',i,',',t,']',sep='')] <- d[,paste('Yihat[',i,',',t,']',sep='')] - jags.dat$Yi[i,t]
    
    for(j in 1:jags.dat$nj[i,t]){
      # resid[,paste('Yjerr[',i,',',t,',',j,']',sep='')] <- d[,paste('Yjhat[',i,',',t,',',j,']',sep='')] - jags.dat$Yj[i,t,j]
      
      if(jags.dat$npasses[i,t,j] > 1){
        for(m in 1:jags.dat$npasses[i,t,j]){
          resid[,paste('yerr[',i,',',t,',',j,',',m,']',sep='')] <- (d[,paste('yhat[',i,',',t,',',j,',',m,']',sep='')] - jags.dat$y[i,t,j,m]) / (jags.dat$sitelength[i,t,j]*1000/30)
        }  
      }
    }
  }
}

# Bias, imprecision, and inaccuracy
acc <- function(x){
  return(mean(abs(x)))
}

# resid[,'Yi.bias'] <- apply(resid[,grepl('Yierr[',names(resid),fixed=T)], 1, mean)
# resid[,'Yi.imp'] <- apply(resid[,grepl('Yierr[',names(resid),fixed=T)], 1, sd)
# resid[,'Yi.acc'] <- apply(resid[,grepl('Yierr[',names(resid),fixed=T)], 1, acc)
# 
# resid[,'Yj.bias'] <- apply(resid[,grepl('Yjerr[',names(resid),fixed=T)], 1, mean)
# resid[,'Yj.imp'] <- apply(resid[,grepl('Yjerr[',names(resid),fixed=T)], 1, sd)
# resid[,'Yj.acc'] <- apply(resid[,grepl('Yjerr[',names(resid),fixed=T)], 1, acc)

resid[,'y.bias'] <- apply(resid[,grepl('yerr[',names(resid),fixed=T)], 1, mean)
resid[,'y.imp'] <- apply(resid[,grepl('yerr[',names(resid),fixed=T)], 1, sd)
resid[,'y.acc'] <- apply(resid[,grepl('yerr[',names(resid),fixed=T)], 1, acc)

save(resid, file='output/objects/resid_alt.R')

## Plot for Yi
jpeg('output/plots/residuals.Yi.jpg',width=6, height=6, units='in', res=600)

layout(matrix(1:3, nrow=3, ncol=1), heights=c(1,1,1))

par(mar=c(4.5,4.5,4,1))

hist(resid[,'Yi.bias'], xlab='Bias', main=expression(Y['i,t']))

par(mar=c(4.5,4.5,1,1))

hist(resid[,'Yi.imp'], xlab='Imprecision', main=NA)

hist(resid[,'Yi.acc'], xlab='Inaccuracy', main=NA)

dev.off()

## Plot Yj
jpeg('output/plots/residuals.Yj.jpg',width=6, height=6, units='in', res=600)

layout(matrix(1:3, nrow=3, ncol=1), heights=c(1,1,1))

par(mar=c(4.5,4.5,4,1))

hist(resid[,'Yj.bias'], xlab='Bias', main=expression(Y['i,t,j']))

par(mar=c(4.5,4.5,1,1))

hist(resid[,'Yj.imp'], xlab='Imprecision', main=NA)

hist(resid[,'Yj.acc'], xlab='Inaccuracy', main=NA)

dev.off()

## Plot y
jpeg('output/plots/residuals.y.jpg',width=6, height=6, units='in', res=600)

layout(matrix(1:3, nrow=3, ncol=1), heights=c(1,1,1))

par(mar=c(4.5,4.5,4,1))

hist(resid[,'y.bias'], xlab='Bias', main=expression(y['i,t,j,m']))

par(mar=c(4.5,4.5,1,1))

hist(resid[,'y.imp'], xlab='Imprecision', main=NA)

hist(resid[,'y.acc'], xlab='Inaccuracy', main=NA)

dev.off()

## Summaries for table
quantile(resid[,'y.bias'], probs=c(0.025, 0.5, 0.975))
quantile(resid[,'y.imp'], probs=c(0.025, 0.5, 0.975))
quantile(resid[,'y.acc'], probs=c(0.025, 0.5, 0.975))

quantile(resid[,'Yj.bias'], probs=c(0.025, 0.5, 0.975))
quantile(resid[,'Yj.imp'], probs=c(0.025, 0.5, 0.975))
quantile(resid[,'Yj.acc'], probs=c(0.025, 0.5, 0.975))

quantile(resid[,'Yi.bias'], probs=c(0.025, 0.5, 0.975))
quantile(resid[,'Yi.imp'], probs=c(0.025, 0.5, 0.975))
quantile(resid[,'Yi.acc'], probs=c(0.025, 0.5, 0.975))

mean(resid[,'y.bias'])
mean(resid[,'y.imp'])
mean(resid[,'y.acc'])

mean(resid[,'Yj.bias'])
mean(resid[,'Yj.imp'])
mean(resid[,'Yj.acc'])

mean(resid[,'Yi.bias'])
mean(resid[,'Yi.imp'])
mean(resid[,'Yi.acc'])




