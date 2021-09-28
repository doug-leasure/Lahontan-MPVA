rm(list=ls())
gc()
cat("\014") 

library(coda);library(boot);library(mc2d)

setwd('D:/RESEARCH/2015 NASA Trout Project/wd/MPVA')

# JAGS data
load('output/objects/jags.dat.R')
npop <- jags.dat$npop
nt <- jags.dat$nt
nj <- jags.dat$nj
npasses <- jags.dat$npasses

Yi <- jags.dat$Yi
Yj <- jags.dat$Yj
y <- jags.dat$y

t_ti <- jags.dat$t_ti

drain <- jags.dat$drain

# JAGS model
load('output/objects/jm.R')

d <- jm$mcmc[[1]]
for (i in 2:length(jm$mcmc)) d <- rbind(d, jm$mcmc[[i]])
d <- data.frame(d, row.names=NULL, check.names=F)

nmcmc <- nrow(d)

rm(jm)
gc()

# Posterior predictions
for (i in 1:npop){
  print(paste(i,'/',npop, sep=''))
  for (t in 1:nt[i]){
    rho <- matrix(NA, nrow=nmcmc, ncol=nj[i,t])
    
    for (j in 1:nj[i,t]){
      
      # p = detection probability
      p <- matrix(NA, nrow=nmcmc, ncol=npasses[i,t,j])
      p[,1] <- inv.logit(d$b0p + d$b1p * drain[i,t,j])
      if(npasses[i,t,j] > 1) {
        for (m in 2:npasses[i,t,j]){
          p[,m] <- p[,1] * exp(d$delta*(m-1))
        }  
      }
      
      # Q = non-detection probability in all PRIOR passes
      Q <- matrix(NA, nrow=nmcmc, ncol=npasses[i,t,j]+1)
      Q[,1] <- 1
      if(npasses[i,t,j] > 1) {
        Q[,2] <- 1-p[,1]
      }
      if(npasses[i,t,j] > 2) {
        for (m in 3:npasses[i,t,j]){
          Q[,m] <- apply(1-p[,1:(m-1)], 1, prod, na.rm=T)
        }
      }
      if(npasses[i,t,j] > 1) {
        Q[,npasses[i,t,j]+1] <- apply(1-p[,1:npasses[i,t,j]], 1, prod, na.rm=T)
      } else {
        Q[,npasses[i,t,j]+1] <- 1-p[,1]
      }
      
      # Yjhat = posterior predictions for Yj
      Yjhat.col <- paste('Yjhat[',i,',',t,',',j,']', sep='')
      YjE.col <- paste('YjE[',i,',',t,',',j,']', sep='')
      n.col <- paste('n[',i,',',t,',',j,']', sep='')
      
      d[,Yjhat.col] <- rbinom(nmcmc, d[,n.col], 1-Q[,npasses[i,t,j]+1])
      d[,YjE.col] <- d[,n.col] * (1-Q[,npasses[i,t,j]+1])
      
      if( npasses[i,t,j] > 1 ){
        
        # pi = probability of detection in pass m AND non-detection in previous passes GIVEN detection occurred one of the passes
        pi <- matrix(NA, nrow=nmcmc, ncol=npasses[i,t,j])
        for(m in 1:npasses[i,t,j]){
          pi[,m] <- p[,m] * Q[,m] / (1-Q[,npasses[i,t,j]+1])
        }
        
        # yhat = posterior predictions for y
        yhat.cols <- paste('yhat[',i,',',t,',',j,',',1:npasses[i,t,j],']', sep='')
        yE.cols <- paste('yE[',i,',',t,',',j,',',1:npasses[i,t,j],']', sep='')
        
        d[,yhat.cols] <- rmultinomial(nmcmc, rep(Yj[i,t,j], nmcmc), pi[, 1:npasses[i,t,j]])
        d[,yE.cols] <- rep(Yj[i,t,j], nmcmc) * pi[, 1:npasses[i,t,j]]
      }
      
      # rho
      psi.col <- paste('psi[',i,',',t,',',j,']', sep='')
      rho[,j] <- d[,psi.col] * (1-Q[,npasses[i,t,j]+1])
      
    }
    
    Yihat.col <- paste('Yihat[',i,',',t,']', sep='')
    YiE.col <- paste('YiE[',i,',',t,']', sep='')
    N.col <- paste('N[',i,',',t_ti[i,t],']', sep='')
    
    d[,Yihat.col] <- rbinom( nmcmc, d[,N.col], 1-apply(1-rho, 1, prod))
    d[,YiE.col] <- d[,N.col] * rho
  }
}

save(d, file='output/objects/d.yhat.R')

# Vectorize data for plots and p-value calculations
y.obs <- y.pred <- y.up <- y.low <- Yj.obs <- Yj.pred <- Yj.up <- Yj.low <- Yi.obs <- Yi.pred <- Yi.up <- Yi.low <- c()

for(i in 1:npop){
  for(t in 1:nt[i]){
    Yi.obs <- c(Yi.obs, Yi[i,t])
    
    x <- d[,paste('Yihat[',i,',',t,']',sep='')]
    Yi.pred <- c(Yi.pred, median(x))
    Yi.up <- c(Yi.up, quantile(x, probs=c(0.975)))
    Yi.low <- c(Yi.low, quantile(x, probs=c(0.025)))
    
    for(j in 1:nj[i,t]){
      Yj.obs <- c(Yj.obs, Yj[i,t,j])
      
      x <- d[,paste('Yjhat[',i,',',t,',',j,']',sep='')]
      Yj.pred <- c(Yj.pred, median(x))
      Yj.up <- c(Yj.up, quantile(x, probs=c(0.975)))
      Yj.low <- c(Yj.low, quantile(x, probs=c(0.025)))
      
      if(npasses[i,t,j] > 1){
        for(m in 1:npasses[i,t,j]){
          y.obs <- c(y.obs, y[i,t,j,m])
          
          x <- d[,paste('yhat[',i,',',t,',',j,',',m,']',sep='')]
          y.pred <- c(y.pred, median(x))
          y.up <- c(y.up, quantile(x, probs=c(0.975)))
          y.low <- c(y.low, quantile(x, probs=c(0.025)))
        }  
      }
    }
  }
}

vectorize.dat <- list(y.obs=y.obs, y.pred=y.pred, y.up=y.up, y.low=y.low, Yj.obs=Yj.obs, Yj.pred=Yj.pred, Yj.up=Yj.up, Yj.low=Yj.low,  Yi.obs=Yi.obs, Yi.pred=Yi.pred, Yi.up=Yi.up, Yi.low=Yi.low)
save(vectorize.dat, file='output/objects/vectorize.dat.R')



#### Difference statistics functions ####

# Sum of squares
sumsq <- function(O, E){
  apply( (O-E)^2 , 1, sum) 
}

# Freeman-Tukey
freetuk <- function(O, E){
  apply( (sqrt(O)-sqrt(E))^2 , 1, sum) 
}

# Proportion zeros
propzero <- function(O, dat=F){
  O <- O==0
  if(dat){
    result <- mean(O)
  } else {
    result <- apply( O, 1, mean )
  }
  return(result)  
}

# Overdispersion
overdisp <- function(O, dat=F, NAlist=NULL){
  if(dat){
    if(length(NAlist)>0) O[NAlist] <- NA
    result <- var(O, na.rm=T) / mean(O, na.rm=T)
  } else {
    if(length(NAlist)>0) {
      O[matrix(NAlist, nrow=nrow(O), ncol=length(NAlist), byrow=T)] <- NA
    }
    result <- apply(O,1,var, na.rm=T) / apply(O,1,mean, na.rm=T)
  }
  return(result)  
}

# Skewness
skewness <- function(O, dat=F, NAlist=NULL){
  if(dat){
    if(length(NAlist)>0) O[NAlist] <- NA
    result <- mean( ((O - mean(O, na.rm=T)) / sd(O, na.rm=T))^3 , na.rm=T)
  } else {
    if(length(NAlist)>0) {
      O[matrix(NAlist, nrow=nrow(O), ncol=length(NAlist), byrow=T)] <- NA
    }
    mu <- apply(O,1,mean, na.rm=T)
    sigma <- apply(O,1,sd, na.rm=T)
    O <- ((O - mu) / sigma )^3
    result <- apply(O,1,mean, na.rm=T)
  }
  return(result)
}

## Calculate difference stats for each response variable

# Yi
hat.cols <- names(d)[grepl('Yihat', names(d))]
E.cols <- names(d)[grepl('YiE', names(d))]
Odat <- Yi.obs
Ohat <- d[,hat.cols]
NAlist <- Odat==0
NAlist <- NULL

Tdat1.Yi <- propzero(Odat, dat=T)
That1.Yi <- propzero(Ohat)

Tdat2.Yi <- overdisp(Odat, dat=T, NAlist=NAlist)
That2.Yi <- overdisp(Ohat, NAlist=NAlist)

Tdat3.Yi <- skewness(Odat, dat=T, NAlist=NAlist)
That3.Yi <- skewness(Ohat, NAlist=NAlist)

Tdat4.Yi <- freetuk(Odat, d[,E.cols])
That4.Yi <- freetuk(Ohat, d[hat.cols])


# Yj
hat.cols <- names(d)[grepl('Yjhat', names(d))]
E.cols <- names(d)[grepl('YjE', names(d))]
Odat <- Yj.obs
Ohat <- d[,hat.cols]
NAlist <- Odat==0
NAlist <- NULL

Tdat1.Yj <- propzero(Odat, dat=T)
That1.Yj <- propzero(Ohat)

Tdat2.Yj <- overdisp(Odat, dat=T, NAlist=NAlist)
That2.Yj <- overdisp(Ohat, NAlist=NAlist)

Tdat3.Yj <- skewness(Odat, dat=T, NAlist=NAlist)
That3.Yj <- skewness(Ohat, NAlist=NAlist)

Tdat4.Yj <- freetuk(Odat, d[,E.cols])
That4.Yj <- freetuk(Ohat, d[hat.cols])

# y
hat.cols <- names(d)[grepl('yhat', names(d))]
E.cols <- names(d)[grepl('yE', names(d))]
Odat <- y.obs
Ohat <- d[,hat.cols]
NAlist <- Odat==0
NAlist <- NULL

Tdat1.y <- propzero(Odat, dat=T)
That1.y <- propzero(Ohat)

Tdat2.y <- overdisp(Odat, dat=T, NAlist=NAlist)
That2.y <- overdisp(Ohat, NAlist=NAlist)

Tdat3.y <- skewness(Odat, dat=T, NAlist=NAlist)
That3.y <- skewness(Ohat, NAlist=NAlist)

Tdat4.y <- freetuk(Odat, d[,E.cols])
That4.y <- freetuk(Ohat, d[hat.cols])

# Save difference stats
T.dat <- list(Tdat1.y=Tdat1.y, That1.y=That1.y, Tdat2.y=Tdat2.y, That2.y=That2.y, Tdat3.y=Tdat3.y, That3.y=That3.y, Tdat4.y=Tdat4.y, That4.y=That4.y, 
              Tdat1.Yj=Tdat1.Yj, That1.Yj=That1.Yj, Tdat2.Yj=Tdat2.Yj, That2.Yj=That2.Yj, Tdat3.Yj=Tdat3.Yj, That3.Yj=That3.Yj, Tdat4.Yj=Tdat4.Yj, That4.Yj=That4.Yj, 
              Tdat1.Yi=Tdat1.Yi, That1.Yi=That1.Yi, Tdat2.Yi=Tdat2.Yi, That2.Yi=That2.Yi, Tdat3.Yi=Tdat3.Yi, That3.Yi=That3.Yi, Tdat4.Yi=Tdat4.Yi, That4.Yi=That4.Yi)
save(T.dat, file='output/objects/T.dat.R')

# Calculate Bayesian p-value
p1.Yi <- mean(That1.Yi > Tdat1.Yi, na.rm=T)
p2.Yi <- mean(That2.Yi > Tdat2.Yi, na.rm=T)
p3.Yi <- mean(That3.Yi > Tdat3.Yi, na.rm=T)
p4.Yi <- mean(That4.Yi > Tdat4.Yi, na.rm=T)

p1.Yj <- mean(That1.Yj > Tdat1.Yj, na.rm=T)
p2.Yj <- mean(That2.Yj > Tdat2.Yj, na.rm=T)
p3.Yj <- mean(That3.Yj > Tdat3.Yj, na.rm=T)
p4.Yj <- mean(That4.Yj > Tdat4.Yj, na.rm=T)

p1.y <- mean(That1.y > Tdat1.y, na.rm=T)
p2.y <- mean(That2.y > Tdat2.y, na.rm=T)
p3.y <- mean(That3.y > Tdat3.y, na.rm=T)
p4.y <- mean(That4.y > Tdat4.y, na.rm=T)

#### Yi ####
p1.Yi
p2.Yi
p3.Yi
p4.Yi

#### Yj ####
p1.Yj
p2.Yj
p3.Yj
p4.Yj

#### y ####
p1.y
p2.y
p3.y
p4.y

ptable <- data.frame(row.names=c('PropZero','Overdisp','Skewness'),
                     Yi=c(p1.Yi, p2.Yi, p3.Yi), 
                     Yj=c(p1.Yj, p2.Yj, p3.Yj),
                     y=c(p1.y, p2.y, p3.y))

ptable <- round(ptable,4)


#### Observed vs. predicted plots ####
library(plotrix)

jpeg('output/plots/obs_pred.jpg',width=6, height=6, units='in', res=600)

layout(matrix(1:6, nrow=3, ncol=2, byrow=T), widths=c(2,1), heights=c(1,1,1))


# mar0 <- c(5, 4, 4, 2) + 0.1

# Yi
par(mar=c(2,4.5,1,1))

plot(NA, xlim=c(min(Yi.obs), max(Yi.obs)), ylim=c(min(Yi.low), max(Yi.up)), xlab=NA, ylab='Predicted', main=NA, cex.lab=1.3)
for(i in 1:length(Yi.pred)){
  arrows(x0=Yi.obs[i], x1=Yi.obs[i], y0=Yi.low[i], y1=Yi.up[i], length=0, lwd=0.5)
}
points(y=Yi.pred, x=Yi.obs, cex=0.5)
abline(0, 1, col='gray', lty=2)

ptab <- data.frame(row.names=row.names(ptable), p_value=ptable[,'Yi'])
addtable2plot(x='topleft', table=ptab, display.rownames=T, display.colnames=T, cex=1)

mtext(expression(Y['i,t']), side=3, line=-2, cex=1.25)

#
par(mar=c(2,2,1,1))

cutoff <- 0.95

xlim <- c(min(Yi.obs), quantile(Yi.obs, probs=c(cutoff)))
ylim <- c(min(Yi.low), max(Yi.up[Yi.obs<xlim[2]]))

plot(NA, xlim=xlim, ylim=ylim, xlab=NA, ylab=NA, main=NA)

for(i in 1:length(Yi.pred)){
  arrows(x0=Yi.obs[i], x1=Yi.obs[i], y0=Yi.low[i], y1=Yi.up[i], length=0, lwd=0.5)
}
jit1 <- 0 # rnorm(length(Yi.pred),0,0.05)
jit2 <- 0 # rnorm(length(Yi.pred),0,0.05)
points(y=Yi.pred+jit1, x=Yi.obs+jit2, col=rgb(0,0,0,1), cex=0.5)
abline(0, 1, col='gray', lty=2)


## Yj ##
par(mar=c(2,4.5,1,1))

plot(NA, xlim=c(min(Yj.obs), max(Yj.obs)), ylim=c(min(Yj.low), max(Yj.up)), xlab=NA, ylab='Predicted', main=NA, cex.lab=1.3)
for(i in 1:length(Yj.pred)){
  arrows(x0=Yj.obs[i], x1=Yj.obs[i], y0=Yj.low[i], y1=Yj.up[i], length=0, lwd=0.5)
}
points(y=Yj.pred, x=Yj.obs, cex=0.5)
abline(0, 1, col='gray', lty=2)

ptab <- data.frame(row.names=row.names(ptable), p_value=ptable[,'Yj'])
addtable2plot(x='topleft', table=ptab, display.rownames=T, display.colnames=T, cex=1)

mtext(expression(Y['i,t,j']), side=3, line=-2, cex=1.25)


#
par(mar=c(2,2,1,1))

xlim <- c(min(Yj.obs), quantile(Yj.obs, probs=c(cutoff)))
ylim <- c(min(Yj.low), max(Yj.up[Yj.obs<xlim[2]]))

plot(NA, xlim=xlim, ylim=ylim, xlab=NA, ylab=NA, main=NA)

for(i in 1:length(Yj.pred)){
  arrows(x0=Yj.obs[i], x1=Yj.obs[i], y0=Yj.low[i], y1=Yj.up[i], length=0, lwd=0.5)
}
jit1 <- 0 # rnorm(length(Yj.pred),0,0.05)
jit2 <- 0 # rnorm(length(Yj.pred),0,0.05)
points(y=Yj.pred+jit1, x=Yj.obs+jit2, col=rgb(0,0,0,1), cex=0.5)
abline(0, 1, col='gray', lty=2)


## y ##
par(mar=c(4,4.5,1,1))

xlim <- c(min(y.obs), max(y.obs))
ylim <- c(min(y.low), max(y.up))

plot(NA, xlim=xlim, ylim=ylim, xlab='Observed', ylab='Predicted', main=NA, cex.lab=1.3)

for(i in 1:length(y.pred)){
  arrows(x0=y.obs[i], x1=y.obs[i], y0=y.low[i], y1=y.up[i], length=0, lwd=0.5)
}
points(y=y.pred, x=y.obs, cex=0.5)
abline(0, 1, col='gray', lty=2)

ptab <- data.frame(row.names=row.names(ptable), p_value=ptable[,'y'])
addtable2plot(x='topleft', table=ptab, display.rownames=T, display.colnames=T, cex=1)

mtext(expression(y['i,t,j,m']), side=3, line=-2, cex=1.25)


#
par(mar=c(4,2,1,1))

xlim <- c(min(y.obs), quantile(y.obs, probs=c(cutoff)))
ylim <- c(min(y.low), max(y.up[y.obs<xlim[2]]))

plot(NA, xlim=xlim, ylim=ylim, xlab=NA, ylab=NA, main=NA)

for(i in 1:length(y.pred)){
  arrows(x0=y.obs[i], x1=y.obs[i], y0=y.low[i], y1=y.up[i], length=0, lwd=0.5)
}
jit1 <- 0 # rnorm(length(y.pred),0,0.05)
jit2 <- 0 # rnorm(length(y.pred),0,0.05)
points(y=y.pred+jit1, x=y.obs+jit2, col=rgb(0,0,0,1), cex=0.5)
abline(0, 1, col='gray', lty=2)

dev.off()









