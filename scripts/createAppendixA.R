rm(list=ls())
gc()
cat("\014") 

setwd('D:/RESEARCH/2015 NASA Trout Project/wd/MPVA')

dir.create('output/appA')

# scripts
outdir <- 'output/appA/scripts/'
dir.create(outdir)

file.copy(from=paste('scripts/02 MPVA.JAGS.R',sep=''), to=outdir, overwrite = T, recursive = T, copy.mode = T)
file.copy(from=paste('scripts/03 inits.R',sep=''), to=outdir, overwrite = T, recursive = T, copy.mode = T)
file.copy(from=paste('scripts/04b MPVA.R',sep=''), to=outdir, overwrite = T, recursive = T, copy.mode = T)
file.copy(from=paste('scripts/zm.init.R',sep=''), to=outdir, overwrite = T, recursive = T, copy.mode = T)

file.rename(from=paste(outdir, '02 MPVA.JAGS.R', sep=''), to=paste(outdir, 'MPVA.JAGS.R',sep=''))
file.rename(from=paste(outdir, '03 inits.R', sep=''), to=paste(outdir, 'inits.R',sep=''))
file.rename(from=paste(outdir, '04b MPVA.R', sep=''), to=paste(outdir, 'MPVA.R',sep=''))

# objects
outdir <- 'output/appA/input/'
dir.create(outdir)

file.copy(from=paste('output/objects/jags.dat.R',sep=''), to=outdir, overwrite = T, recursive = T, copy.mode = T)
