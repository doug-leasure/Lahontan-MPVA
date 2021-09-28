rm(list=ls())
gc()
cat("\014") 

setwd('D:/RESEARCH/2015 NASA Trout Project/wd/MPVA')

dir.create('output')

source(file='scripts/00a fillBKT.R', echo=T)
source(file='scripts/00b extract.flowdat.R', echo=T)
source(file='scripts/01a initialize data.R', echo=T)
source(file='scripts/01b setup jags data.R', echo=T)

source(file='scripts/04 MPVA.R', echo=T)

