rm(list=ls())
gc()
cat("\014") 

# library(rjags);library(runjags);library(coda);library(lattice);library(boot);
library(runjags);library(boot);library(rjags)

dir.create('output/plots')
dir.create('output/tables')
dir.create('output/objects')

# Load JAGS data
load('output/objects/jags.dat.R')
attach(jags.dat, warn.conflicts=F)

# Parameters to monitor (Shiny app needs N[i,t])
par.monitor <- c('b0r','b1r','b2r','b0phi','b1phi','b2phi','musig','sigsig','tau','b0p','b1p','delta')

# JAGS setup
print(getwd())
set.seed(123)
load.module('lecuyer')
load.module('glm')

n.adapt <- 1e3
n.burn <- 0
n.iter <- 1e3
thin <- 100
n.chains <- 10

source('scripts/03 inits.R')
init <- inits(chains=n.chains)
save(init, file='output/objects/init.R')

# Run JAGS
jm <- run.jags(model='scripts/02 MPVA.JAGS.R', 
               monitor=par.monitor, 
               data=jags.dat, 
               n.chains=n.chains, 
               inits=init, 
               thin=thin,
               adapt=n.adapt,
               burnin=n.burn,
               sample=n.iter,
               summarise=F, 
               #keep.jags.files=TRUE,
               method='parallel' #'rjags' #'rjparallel'
               )

# Time taken in hours
as.numeric(jm$timetaken)/3600 

# rm(i, init, jags.dat, modelfile, n.adapt, n.burn, n.chains, par.monitor, t, thin, inits)
jm0 <- jm
save(jm0, file='output/objects/jm0.R')
rm(jm0)

# zm.init = A function that puts the initials in as the first row in the CODA file
source('scripts/zm.init.R')

# Make a CODA file
zm <- zm.init(jm$mcmc, init, burn=0, thinmore=1, addinit=T)

# Trace plots
dir.create('output/plots/trace.plots')
jpeg(paste('output/plots/trace.plots/trace.%02d.jpg', sep=''), quality=100, height=1000, width=1000)
vars <- c('b0r','b1r','b2r','b0phi','b1phi','b2phi[1]','b2phi[2]','musig','sigsig','tau','b0p','b1p','delta')
for (i in seq(1, length(vars), 3)){
  plot(zm[,vars[i:min(length(vars), i+2)]], cex=1.25)
}
dev.off()

# Gelman-Rubin convergence statistic
zm <- zm.init(jm$mcmc, init, burn=5e4, thinmore=1, addinit=F)
psrf <- gelman.diag(zm, multivariate=F)

# Add additional parameters to monitor after burnin
add.monitor <- paste('sigmaR[',1:npops,']',sep='')
for (i in 1:npops){
  for (t in (t_ti[i,1]-1):t_ti[i,nt[i]]){
    add.monitor <- c(add.monitor, paste('N[',i,',',t,']', sep=''))
  }
  for(t in 1:nt[i]){
    for(j in 1:nj[i,t]){
      add.monitor <- c(add.monitor, paste('psi[',i,',',t,',',j,']', sep=''))
      add.monitor <- c(add.monitor, paste('n[',i,',',t,',',j,']', sep=''))
    }
  }
}

# Re-run until convergence is achieved
iteration <- 1
n.iter <- 1000

while(max(psrf$psrf[,'Upper C.I.']) > 1.1){
  
  print(paste('Iteration:', iteration))
  
  if(sum(!add.monitor %in% jm$monitor) == 0) { 
    jm <- extend.jags(jm, sample=n.iter)
  } else {
    jm <- extend.jags(jm, sample=n.iter, add.monitor=add.monitor, thin=thin)
  }
  
  as.numeric(jm$timetaken)/3600 # time taken in hours
  
  save(jm, file='output/objects/jm.R')
  
  iteration <- iteration + 1
  
  thinmore <- 1 # round(nrow(jm$mcmc[[1]]) / 5000)
  zm <- zm.init(jm$mcmc, init, burn=5e5, thinmore=thinmore, addinit=F)
  
  # pdf('output/plots/trace.plots.pdf')
  # vars <- c('b0r','b1r','b2r','b0phi','b1phi','b2phi[1]','b2phi[2]','musig','sigsig','tau','b0p','b1p','delta',add.monitor[grepl('sigmaR',add.monitor,fixed=T)])
  # for (i in seq(1, length(vars), 4)){
  #   plot(zm[,vars[i:min(length(vars), i+3)]], cex=1.25, density=F, smooth=F)
  # }
  # dev.off()

  vars.psrf <- c('b0r','b1r','b2r','b0phi','b1phi','b2phi[1]','b2phi[2]','musig','sigsig','tau','b0p','b1p','delta')
  psrf <- gelman.diag(zm[,vars.psrf])
  write.csv(as.data.frame(psrf$psrf), 'output/tables/psrf.csv')

  rm(zm)
  gc()
}

## Save final model
save(jm, file='output/objects/jm.R')

zm <- zm.init(jm$mcmc, init, burn=1e6, thinmore=1, addinit=F, varlist=varnames(jm$mcmc))
gelman.diag(zm.thin[,c('b0r','b1r','b2r','b0phi','b1phi','b2phi[1]','b2phi[2]','musig','sigsig','tau','b0p','b1p','delta')])$psrf

zm.thin <- zm.init(jm$mcmc, init, burn=1e6, thinmore=5, addinit=F)
gelman.diag(zm.thin[,c('b0r','b1r','b2r','b0phi','b1phi','b2phi[1]','b2phi[2]','musig','sigsig','tau','b0p','b1p','delta')])$psrf

varlist <- varnames(jm$mcmc)
varlist <- varlist[!grepl('n[', varlist, fixed=T) & !grepl('psi[',varlist,fixed=T)]

zm.thin <- zm.init(jm$mcmc, init, burn=1e6, thinmore=5, addinit=F, varlist=varlist)
save(zm.thin, file='output/objects/zm.thin.R')

d <- zm.thin[[1]]
for (i in 2:length(zm.thin)) d <- rbind(d, zm.thin[[i]])
d <- data.frame(d, row.names=NULL, check.names=F)
save(d, file='output/objects/dthin.R')

pdf('output/plots/trace.plots.pdf')
vars <- c('b0r','b1r','b2r','b0phi','b1phi','b2phi[1]','b2phi[2]','musig','sigsig','tau','b0p','b1p','delta',add.monitor[grepl('sigmaR',add.monitor,fixed=T)])
for (i in seq(1, length(vars), 4)){
  plot(zm.thin[,vars[i:min(length(vars), i+3)]], cex=1.25, density=F, smooth=F)
}
dev.off()

vars.psrf <- c('b0r','b1r','b2r','b0phi','b1phi','b2phi[1]','b2phi[2]','musig','sigsig','tau','b0p','b1p','delta',paste('sigmaR[',1:npops,']',sep=''))
psrf <- gelman.diag(jm$mcmc[,vars.psrf])
psrf$psrf
write.csv(as.data.frame(psrf$psrf), 'output/tables/psrf.csv')

zm.summary <- summary(jm$mcmc, quantiles=c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975))
write.csv(zm.summary$statistics, file='output/tables/zm.statistics.csv')
write.csv(zm.summary$quantiles, file='output/tables/zm.quantiles.csv')







# Thin the model more and re-check convergence
newthin <- 10
zm.thin <- mcmc.list()

for(i in 1:length(zm)){
  zm.thin[[i]] <- as.mcmc(zm[[i]][seq(1, nrow(zm[[i]]), by=newthin),])
}
dim(zm.thin[[1]])

psrf.thin <- gelman.diag(zm)
psrf.thin

zm.thin.summary <- summary(zm.thin, quantiles=c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975))

zm.summary$quantiles[1:15,]
zm.thin.summary$quantiles[1:15,]

zm.valyear.thin <- zm.thin
save(zm.valyear.thin, file='output/objects/zm.valyear.thin.R')

write.csv(zm.thin.summary$statistics, file='output/tables/zm.valyear.thin.statistics.csv')
write.csv(zm.thin.summary$quantiles, file='output/tables/zm.valyear.thin.quantiles.csv')
write.csv(psrf.thin$psrf, file='output/tables/psrf.valyear.thin.csv')









#############################################

# Load original JAGS data
detach(jags.dat)
load('output/objects/jags.dat.R')
attach(jags.dat)

# Convert model to MCMC list
zm <- zm.thin

# Variables to assess
vars <- c('b0r','b1r','b2r','b3r',
          'b0phi','b1phi','b2phi[1]','b2phi[2]',
          'b0palpha','b1palpha','pbeta',
          'sigmaprop','sigmap',
          'alpha','beta',paste('sigmaR[',1:npops,']', sep=''),
)


# Convert MCMC list to data frame
d <- data.frame()
for (i in 1:length(zm)) d <- rbind(d, zm[[i]])
save(d, file='output/objects/d.val.year.R')
rm(zm)
gc()


# Observed vs predicted plots: y
ylim <- c(Inf, 0)
xlim <- c(Inf, 0)
for(i in 1:npops){
  if(nt[i] > 2){
    for(t in nt[i]:nt[i]){
      for(j in 1:nj[i,t]){
        for(m in 1:npasses[i,t,j]){
          if(y[i,t,j,m] < xlim[1]) xlim[1] <- y[i,t,j,m]
          if(y[i,t,j,m] > xlim[2]) xlim[2] <- y[i,t,j,m]
          col <- paste('y[',i,',',t,',',j,',',m,']',sep='')
          quantlow <- quantile(d[,col], probs=c(0.025))
          quanthigh <- quantile(d[,col], probs=c(0.975))
          if(quantlow < ylim[1]) ylim[1] <- quantlow
          if(quanthigh > ylim[2]) ylim[2] <- quanthigh
        }
      }
    }
  }
}

jpeg('output/plots/fitplot.y.valyear.jpg')
plot(NA, xlim=xlim, ylim=ylim, ylab='Predicted', xlab='Observed', main='y')
abline(0, 1, col='red')
for (i in 1:npops){
  if(nt[i] > 2){
    for (t in nt[i]:nt[i]){
      for (j in 1:nj[i,t]){
        for(m in 1:npasses[i,t,j]){
          col <- paste('y[',i,',',t,',',j,',',m,']',sep='')
          quantlow <- quantile(d[,col], probs=c(0.025))
          quanthigh <- quantile(d[,col], probs=c(0.975))
          
          points(x = y[i,t,j,m], y=mean(d[,col]))
          arrows(x0=y[i,t,j,m], y0=quantlow, y1=quanthigh, length=0)
        }
      }
    }
  }
}
dev.off()

# Observed vs predicted plots: Yj
ylim <- c(Inf, 0)
xlim <- c(Inf, 0)
for(i in 1:npops){
  if(nt[i] > 2){
    for(t in nt[i]:nt[i]){
      for(j in 1:nj[i,t]){
        if(Yj[i,t,j] < xlim[1]) xlim[1] <- Yj[i,t,j]
        if(Yj[i,t,j] > xlim[2]) xlim[2] <- Yj[i,t,j]
        col <- paste('Yj[',i,',',t,',',j,']',sep='')
        quantlow <- quantile(d[,col], probs=c(0.025))
        quanthigh <- quantile(d[,col], probs=c(0.975))
        if(quantlow < ylim[1]) ylim[1] <- quantlow
        if(quanthigh > ylim[2]) ylim[2] <- quanthigh
      }
    }
  }
}

jpeg('output/plots/fitplot.Yj.valyear.jpg')
plot(NA, xlim=xlim, ylim=ylim, ylab='Predicted', xlab='Observed', main='Yj')
abline(0, 1, col='red')
for (i in 1:npops){
  if(nt[i] > 2){
    for (t in nt[i]:nt[i]){
      for (j in 1:nj[i,t]){
        col <- paste('Yj[',i,',',t,',',j,']',sep='')
        quantlow <- quantile(d[,col], probs=c(0.025))
        quanthigh <- quantile(d[,col], probs=c(0.975))
        
        points(x = Yj[i,t,j], y=mean(d[,col]))
        arrows(x0=Yj[i,t,j], y0=quantlow, y1=quanthigh, length=0)
      }
    }
  }
}
dev.off()

# Observed vs predicted plots: Yi
ylim <- c(Inf, 0)
xlim <- c(Inf, 0)
for(i in 1:npops){
  if(nt[i] > 2){
    for(t in nt[i]:nt[i]){
      if(Yi[i,t] < xlim[1]) xlim[1] <- Yi[i,t]
      if(Yi[i,t] > xlim[2]) xlim[2] <- Yi[i,t]
      col <- paste('Yi[',i,',',t,']',sep='')
      quantlow <- quantile(d[,col], probs=c(0.025))
      quanthigh <- quantile(d[,col], probs=c(0.975))
      if(quantlow < ylim[1]) ylim[1] <- quantlow
      if(quanthigh > ylim[2]) ylim[2] <- quanthigh
    }
  }
}

jpeg('output/plots/fitplot.Yi.valyear.jpg')
plot(NA, xlim=xlim, ylim=ylim, ylab='Predicted', xlab='Observed', main='Yi')
abline(0, 1, col='red')
for (i in 1:npops){
  if(nt[i] > 2){
    for (t in nt[i]:nt[i]){
      col <- paste('Yi[',i,',',t,']',sep='')
      quantlow <- quantile(d[,col], probs=c(0.025))
      quanthigh <- quantile(d[,col], probs=c(0.975))
      
      points(x = Yi[i,t], y=mean(d[,col]))
      arrows(x0=Yi[i,t], y0=quantlow, y1=quanthigh, length=0)
    }
  }
}
dev.off()







# Cleanup
rm(list=ls()[-which(ls() %in% c('jm'))])
