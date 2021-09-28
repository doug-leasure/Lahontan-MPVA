install.packages('runjags', dependencies=T)
library('runjags')

setwd('C:\working_directory') # <-- set this to whatever folder has the files

load('jm.Robj')

jm.thin <- combine.mcmc(jm, thin=3)

save(jm.thin, file='jm.thin.Robj')

