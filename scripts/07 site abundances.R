rm(list=ls())
gc()
cat("\014") 

library(runjags);library(boot);library(rjags)

setwd('D:/RESEARCH/2015 NASA Trout Project/wd/MPVA')

# jm
load('output/objects/jm.R')

# zm
source('scripts/zm.init.R')
varlist <- varnames(jm$mcmc)[grepl('n[', varnames(jm$mcmc), fixed=T)]
zm <- zm.init(jm$mcmc, init, burn=1e6, thinmore=1, addinit=F, varlist=varlist)

# d
d <- zm[[1]]
for (i in 2:length(zm)) d <- rbind(d, zm[[i]])
d <- data.frame(d, row.names=NULL, check.names=F)

# jags.dat
load('output/objects/jags.dat.R')

# popnames
load('output/objects/popnames.R')

# metapop
metapop <- read.csv('input/LCT_MetaPop_Pop.20171101.csv', stringsAsFactors=F)

# lct.dat with original columns
lct.dat <- read.csv('output/tables/lct.dat.origcols.csv', stringsAsFactors=F)

# setup site abundance data frame
cols <- c('PopID','PopulationName','Year','SiteID','SiteLength_m','n.mean','n.median','n.low','n.up')
n <- data.frame(matrix(NA, nrow=0, ncol=length(cols)))
names(n) <- cols

# iterate through all sites
rowcount <- 0
for(i in 1:jags.dat$npops){
  popname <- popnames[i]
  popid <- metapop[metapop$PopulationName==popname, 'PopID']
  
  print(popname)
  
  for(t in 1:jags.dat$nt[i]){
    year <- 1983 + jags.dat$t_ti[i,t]
    
    for(j in 1:jags.dat$nj[i,t]){
      siteid <- lct.dat[lct.dat$pop==i & lct.dat$t==t & lct.dat$j==j, 'SiteID'][1]
      sitelength <- lct.dat[lct.dat$pop==i & lct.dat$t==t & lct.dat$j==j, 'SiteLength_m'][1]
      
      dat <- d[,paste('n[',i,',',t,',',j,']',sep='')]
      
      n.mean <- round(mean(dat), 1)
      n.median <- median(dat)
      n.up <- quantile(dat, probs=c(0.975))
      n.low <- quantile(dat, probs=c(0.025))
      
      rowcount <- rowcount + 1
      n[rowcount, cols] <- c(popid, popname, year, siteid, sitelength, n.mean, n.median, n.low, n.up)
    }
  }
}

# Write to file
write.csv(n, file='output/tables/site.abundances.csv', row.names=F)
