rm(list=ls())
gc()
#cat("\014") 

library(rjags);library(runjags);library(coda);library(lattice);library(boot);library(beepr)

setwd('D:/RESEARCH/2015 NASA Trout Project/wd/MPVA')

dir.create('output/plots')
dir.create('output/tables')
dir.create('output/objects')
dir.create('output/setup')

# Load JAGS data
load('output/objects/jags.dat.R')

# JAGS setup
print(getwd())
set.seed(123)
load.module('lecuyer')
load.module('glm')

n.adapt <- 1e3
n.burn <- 0
n.iter <- 1e3
thin <- 500
n.chains <- 5

# JAGS initials
source('scripts/03 inits.R')
init <- inits(chains=n.chains)
save(init, file='output/objects/init.xvalyear.R')

# Remove last year of data for each population
for(i in 1:jags.dat$npops){
  if(jags.dat$nt[i] > 2){
    jags.dat$y[i,jags.dat$nt[i],,] <- NA
    jags.dat$Yj[i,jags.dat$nt[i],] <- NA
    jags.dat$Yi[i,jags.dat$nt[i]] <- NA
    # jags.dat$nt[i] <- jags.dat$nt[i] - 1
  }
}

# Parameters to monitor
par.monitor <- c('b0r','b1r','b2r','b0phi','b1phi','b2phi','musig','sigsig','tau','b0p','b1p','delta','sigmaR')


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

# rm(i, init, jags.dat, modelfile, n.adapt, n.burn, n.chains, par.monitor, t, thin, inits)
jm0 <- jm
save(jm0, file='output/objects/jm0.xvalyear.R')
rm(jm0)

# zm.init = A function that puts the initials in as the first row in the CODA file
source('scripts/zm.init.R')

# Make a CODA file
zm <- zm.init(jm$mcmc, init, burn=0, thinmore=1, addinit=T)

# Trace plots
vars <- c('b0r','b1r','b2r','b0phi','b1phi','b2phi[1]','b2phi[2]','musig','sigsig','tau','b0p','b1p','delta')
for (i in seq(1, length(vars), 3)){
  plot(zm[,vars[i:min(length(vars), i+2)]], cex=1.25)
}

# Gelman-Rubin convergence statistic
zm <- zm.init(jm$mcmc, init, burn=2e5, thinmore=1, addinit=F)
psrf <- gelman.diag(zm, multivariate=F)

# Add monitor
add.monitor <- c()
for (i in 1:jags.dat$npops){
  add.monitor <- c(add.monitor, paste('Yi[',i,',',jags.dat$nt[i],']', sep=''))
  
  if(nt[i] > 2){
    for (j in 1:jags.dat$nj[i,jags.dat$nt[i]]){
      add.monitor <- c(add.monitor, paste('Yj[',i,',',jags.dat$nt[i],',',j,']', sep=''))
      
      for (m in 1:npasses[i,jags.dat$nt[i],j]){
        add.monitor <- c(add.monitor, paste('y[',i,',',jags.dat$nt[i],',',j,',',m,']', sep=''))
      }
    }  
  }
}



# Re-run until convergence is achieved
iteration <- 1
n.iter <- n.iter
thin <- thin

while(max(psrf$psrf[,'Upper C.I.']) > 1.1){
  
  print(paste('Iteration:', iteration))
  
  if(sum(!add.monitor %in% jm$monitor) == 0) { 
    jm <- extend.jags(jm, sample=n.iter)
  } else {
    jm <- extend.jags(jm, sample=n.iter, add.monitor=add.monitor, thin=thin)
  }
  
  save(jm, file='output/objects/jm.xvalyear.R')
  
  iteration <- iteration + 1
  
  thinmore <- 1 # round(nrow(jm$mcmc[[1]]) / 5000)
  burn <- 1e6
  mcmc.count <- as.numeric(row.names(jm$mcmc[[1]])[nrow(jm$mcmc[[1]])])
  if (burn >= mcmc.count) {
    burn <- round(mcmc.count/2)
  }
  zm <- zm.init(jm$mcmc, init, burn=burn, thinmore=thinmore, addinit=F)
  
  vars.psrf <- c('b0r','b1r','b2r','b0phi','b1phi','b2phi[1]','b2phi[2]','musig','sigsig','tau','b0p','b1p','delta')
  psrf <- gelman.diag(zm[,vars.psrf])
  write.csv(as.data.frame(psrf$psrf), 'output/tables/psrf.xvalyear.csv')
  
  rm(zm)
  gc()
}




zm <- zm.init(jm$mcmc, init, burn=burn, addinit=F, varlist=varnames(jm$mcmc))

varlist <- varnames(jm$mcmc)
varlist <- varlist[!grepl('Yi[', varlist, fixed=T) & !grepl('Yj[',varlist,fixed=T) & !grepl('y[',varlist,fixed=T)]

zm.summary <- summary(zm[,varlist], quantiles=c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975))

psrf <- gelman.diag(zm[,varlist])

d <- jm$mcmc[[1]]
for (i in 2:length(jm$mcmc)) d <- rbind(d, jm$mcmc[[i]])
d <- data.frame(d, row.names=NULL, check.names=F)

# Save Results
save(d, file='output/objects/d.xvalyear.R')

jm.valyear <- jm
zm.valyear <- zm

save(jm.valyear, file='output/objects/jm.xvalyear.R')
save(zm.valyear, file='output/objects/zm.xvalyear.R')

write.csv(zm.summary$statistics, file='output/tables/zm.xvalyear.statistics.csv')
write.csv(zm.summary$quantiles, file='output/tables/zm.xvalyear.quantiles.csv')
write.csv(psrf$psrf, file='output/tables/psrf.xvalyear.csv')

