# Load JAGS model
load('output/objects/jm0.R')
jm <- jm0

# Make a CODA file
source('scripts/zm.init.R')
zm <- zm.init(jm$mcmc, init, burn=0, thinmore=1, addinit=T)

# Trace plots
dir.create('output/plots/trace.plots')
jpeg(paste('output/plots/trace.plots/trace.%02d.jpg', sep=''), quality=100, height=1000, width=1000)
vars <- c('b0r','b1r','b2r','b0phi','b1phi','b2phi[1]','b2phi[2]','musig','sigsig','tau','b0p','b1p','delta')
for (i in seq(1, length(vars), 3)){
  plot(zm[,vars[i:min(length(vars), i+2)]], cex=1.25)
}
dev.off()


# Create list of runjags objects
jm.list <- list()
for(iteration in 1:13){
  load(paste('output/objects/jm',iteration,'.R',sep=''))
  jm.list[[iteration]] <- jm
}

jm <- combine.mcmc(jm.list)

dir.create('output/plots/trace.plots')
jpeg(paste('output/plots/trace.plots/trace.%02d.jpg', sep=''), quality=100, height=1000, width=1000)
vars <- c('b0r','b1r','b2r','b0phi','b1phi','b2phi[1]','b2phi[2]','musig','sigsig','tau','b0p','b1p','delta',paste('sigmaR[',1:npops,']',sep=''))
for (i in seq(1, length(vars), 3)){
  plot(zm[,vars[i:min(length(vars), i+2)]], cex=1.25)
}
dev.off()


