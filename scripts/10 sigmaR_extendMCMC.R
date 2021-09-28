rm(list=ls())
gc()
cat("\014") 

# library(rjags);library(runjags);library(coda);library(lattice);library(boot);
library(runjags);library(boot);library(rjags)

load('output/objects/jags.dat.R')
attach(jags.dat, warn.conflicts=F)

# Parameters to monitor (Shiny app needs N[i,t])
par.monitor <- c('musig','sigsig','sigmaR0')

# JAGS setup
print(getwd())
set.seed(123)
load.module('lecuyer')
load.module('glm')

n.adapt <- 1e3
n.burn <- 0
n.iter <- 1e4
thin <- 10
n.chains <- 10

load('output/objects/init.R')

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

source('scripts/zm.init.R')

# Make a CODA file
zm <- zm.init(jm$mcmc, init, burn=0, thinmore=1, addinit=F, varlist=par.monitor)

# Trace plots
vars <- par.monitor
for (i in seq(1, length(vars), 3)){
  plot(zm[,vars[i:min(length(vars), i+2)]], cex=1.25)
}

# If burn-in complete keep model object separate
jm0 <- jm

# Gelman-Rubin convergence statistic
zm <- zm.init(jm$mcmc, init, burn=5e4, thinmore=1, addinit=F, varlist=vars)
psrf <- gelman.diag(zm, multivariate=F)


# Change monitors
add.monitor <- c('sigmaR0')
# drop.monitor <- varnames(jm$mcmc)

# Setup JAGS extension
n.iter <- 1e5
thin <- 10

# Loop through extensions until convergence is achieved
iteration <- 1

while(max(psrf$psrf[,'Upper C.I.']) > 1.1){
  
  print(paste('Iteration:', iteration))
  
  if(sum(!add.monitor %in% jm$monitor) == 0) { 
    jm <- extend.jags(jm, sample=n.iter)
  } else {
    jm <- extend.jags(jm, sample=n.iter, add.monitor=add.monitor, drop.monitor=drop.monitor, thin=thin, model='scripts/02 MPVA.JAGS.R')
  }
  
  iteration <- iteration + 1
  
  zm <- zm.init(jm$mcmc, init, burn=1e4, addinit=F)
  
  vars.psrf <- c('sigmaR0')
  psrf <- gelman.diag(zm[,vars.psrf])
  
  rm(zm)
  gc()
}

