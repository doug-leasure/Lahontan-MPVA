rm(list=ls())
gc()
cat("\014") 

setwd('D:/RESEARCH/2015 NASA Trout Project/wd/MPVA')

library(coda)

load('output/objects/jm.R')

varnames(jm$mcmc)

d <- jm$mcmc[[1]]
for(i in 2:length(jm$mcmc)){
  print(i)
  d <- rbind(d, jm$mcmc[[i]])
}
d <- data.frame(d, row.names=NULL, check.names=F)

lct.mpva <- apply(d, 2, median)
save(lct.mpva, file='output/objects/lct.mpva.median.RData')

lct.mpva <- apply(d, 2, mean)
save(lct.mpva, file='output/objects/lct.mpva.mean.RData')
