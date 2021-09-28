#rm(list=ls())
#cat("\014") 
gc()
options(scipen=5) # stop R from converting to scientific notation so often

library(rjags);library(runjags);library(coda);library(lattice);library(boot)

dir.create('output')
dir.create('output/plots')
dir.create('output/plots/cov.effects')

# Load model
load('output/objects/jags.dat.R')
attach(jags.dat, warn.conflicts=F)

rescale <- read.csv('output/tables/scale.factors.csv', stringsAsFactors=F)

load('output/objects/d.R')


# r ~ temp
rng <- range(templag)
x <- seq(rng[1], rng[2], length=100)

p <- data.frame(matrix(NA, nrow=nrow(d), ncol=length(x)))
for (i in 1:length(x)) {
  p[,i] <- d$b0r + d$b1r * x[i] + d$b2r * x[i]^2
}

pmed <- apply(p, 2, median)
pup <- apply(p, 2, quantile, probs=c(0.975))
plow <- apply(p, 2, quantile, probs=c(0.025))

jpeg('output/plots/cov.effects/r.templag.jpg')
plot(x=x, y=pmed, type='l', ylim=c(min(plow), max(pup)), lwd=2, xlab='Temperature (t-1)', ylab='Pop Growth Rate (r)', xaxt='n')
lines(x=x, y=pup, lty=2)
lines(x=x, y=plow, lty=2)
at <- seq(round(min(x)), round(max(x)))
labels <- round(at * rescale[rescale$cov=='temp', 'sd'] + rescale[rescale$cov=='temp', 'mu'])
axis(1, at=at, labels=labels)
dev.off()

# r ~ hflow
rng <- range(hflowlag)
x <- seq(rng[1], rng[2], length=100)

p <- data.frame(matrix(NA, nrow=nrow(d), ncol=length(x)))
for (i in 1:length(x)) {
  p[,i] <- d$b0r + d$b3r * x[i]
}

pmed <- apply(p, 2, median)
pup <- apply(p, 2, quantile, probs=c(0.975))
plow <- apply(p, 2, quantile, probs=c(0.025))

jpeg('output/plots/cov.effects/r.hflowlag.jpg')
plot(x=x, y=pmed, type='l', ylim=c(min(plow), min(10,max(pup))), lwd=2, xlab='High flow (t-1) scaled per pop', ylab='Pop Growth Rate (r)')
lines(x=x, y=pup, lty=2)
lines(x=x, y=plow, lty=2)
dev.off()


# phi ~ bkt
rng <- range(bkt)
x <- seq(rng[1], rng[2], length=100)

p <- data.frame(matrix(NA, nrow=nrow(d), ncol=length(x)))
for (i in 1:length(x)) {
  p[,i] <- d$b0phi + d$b1phi * x[i]
}

pmed <- apply(p, 2, median)
pup <- apply(p, 2, quantile, probs=c(0.975))
plow <- apply(p, 2, quantile, probs=c(0.025))

jpeg('output/plots/cov.effects/phi.bkt.jpg')
plot(x=x, y=pmed, type='l', ylim=c(min(plow), max(pup)), lwd=2, xlab='Brook trout per meter', ylab='Strength of density-dependence (phi)', xaxt='n')
lines(x=x, y=pup, lty=2)
lines(x=x, y=plow, lty=2)
at <- seq(round(min(x)), round(max(x)))
labels <- round(at * rescale[rescale$cov=='bkt', 'sd'], 1)
axis(1, at=at, labels=labels)
dev.off()


# phi ~ pndvi
rng <- range(ndvi)
x <- seq(rng[1], rng[2], length=100)

p <- data.frame(matrix(NA, nrow=nrow(d), ncol=length(x)))
for (i in 1:length(x)) {
  p[,i] <- d$b0phi + d$b2phi * x[i]
}

pmed <- apply(p, 2, median)
pup <- apply(p, 2, quantile, probs=c(0.975))
plow <- apply(p, 2, quantile, probs=c(0.025))

jpeg('output/plots/cov.effects/phi.pndvi.jpg')
plot(x=x, y=pmed, type='l', ylim=c(min(plow), max(pup)), lwd=2, xlab='NDVI (scaled)', ylab='Strength of density-dependence (phi)')
lines(x=x, y=pup, lty=2)
lines(x=x, y=plow, lty=2)

# at <- seq(round(min(x)), round(max(x)))
# labels <- round(at * rescale[rescale$cov=='gndvi', 'sd'] + rescale[rescale$cov=='gndvi', 'mu'], 1)
# axis(1, at=at, labels=labels)

dev.off()

# phi ~ burn
rng <- range(burn)
x <- seq(rng[1], rng[2], length=100)

p <- data.frame(matrix(NA, nrow=nrow(d), ncol=length(x)))
for (i in 1:length(x)) {
  p[,i] <- d$b0phi + d$b3phi * x[i]
}

pmed <- apply(p, 2, median)
pup <- apply(p, 2, quantile, probs=c(0.975))
plow <- apply(p, 2, quantile, probs=c(0.025))

jpeg('output/plots/cov.effects/phi.burn.jpg')
plot(x=x, y=pmed, type='l', ylim=c(min(plow), max(pup)), lwd=2, xlab='Prop. drainage burned', ylab='Strength of density-dependence (phi)')
lines(x=x, y=pup, lty=2)
lines(x=x, y=plow, lty=2)
dev.off()


# p ~ drain
rng <- range(drain, na.rm=T)
x <- seq(rng[1], rng[2], length=100)

p <- data.frame(matrix(NA, nrow=nrow(d), ncol=length(x)))
for (i in 1:length(x)) {
  p[,i] <- inv.logit(d$b0palpha + d$b1palpha * x[i])
}

pmed <- apply(p, 2, median)
pup <- apply(p, 2, quantile, probs=c(0.975))
plow <- apply(p, 2, quantile, probs=c(0.025))

jpeg('output/plots/cov.effects/p.drain.jpg')
plot(x=x, y=pmed, type='l', ylim=c(min(plow), max(pup)), lwd=2, xlab='Drainage Area (sq km)', ylab='Detection Probability (p)', xaxt='n')
lines(x=x, y=pup, lty=2)
lines(x=x, y=plow, lty=2)
at <- seq(round(min(x)), round(max(x)))
labels <- round(at * rescale[rescale$cov=='drain', 'sd'] + rescale[rescale$cov=='drain', 'mu'], 1)
axis(1, at=at, labels=labels)
dev.off()

# p ~ lflow
rng <- range(lflow, na.rm=T)
x <- seq(rng[1], rng[2], length=100)

p <- data.frame(matrix(NA, nrow=nrow(d), ncol=length(x)))
for (i in 1:length(x)) {
  p[,i] <- inv.logit(d$b0palpha + d$b2palpha * x[i])
}

pmed <- apply(p, 2, median)
pup <- apply(p, 2, quantile, probs=c(0.975))
plow <- apply(p, 2, quantile, probs=c(0.025))

jpeg('output/plots/cov.effects/p.lflow.jpg')
plot(x=x, y=pmed, type='l', ylim=c(min(plow), max(pup)), lwd=2, xlab='Mean Summer Flow (scaled)', ylab='Detection Probability (p)')
lines(x=x, y=pup, lty=2)
lines(x=x, y=plow, lty=2)
dev.off()

# Cleanup
rm(list=ls()[-which(ls() %in% c('d','zm'))])
