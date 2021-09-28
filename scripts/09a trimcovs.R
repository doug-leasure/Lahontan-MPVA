#rm(list=ls())
#cat("\014") 
gc()

dir.create('output')
dir.create('output/tables')
dir.create('output/tables/trimcovs')

# Covariates
extent.dat <- read.csv('output/tables/covs/extent.csv', stringsAsFactors=F)
row.names(extent.dat) <- extent.dat$PopulationName

add.dat <- read.csv('output/tables/covs/add.csv', stringsAsFactors=F)
rem.dat <- read.csv('output/tables/covs/rem.csv', stringsAsFactors=F)

row.names(add.dat) <- add.dat$PopulationName
row.names(rem.dat) <- rem.dat$PopulationName

# r covs
templag.dat <- read.csv('output/tables/covs/templag.scale.csv', stringsAsFactors=F)
hflowlag.dat <- read.csv('output/tables/covs/qmax3lag.scale.csv', stringsAsFactors=F)
lflowlag.dat <- read.csv('output/tables/covs/mslag.scale.csv', stringsAsFactors=F)

row.names(templag.dat) <- templag.dat$PopulationName
row.names(hflowlag.dat) <- hflowlag.dat$PopulationName
row.names(lflowlag.dat) <- lflowlag.dat$PopulationName

# phi covs
bkt.dat <- read.csv('output/tables/covs/bkt.scale.csv', stringsAsFactors=F)
lflow.dat <- read.csv('output/tables/covs/ms.scale.csv', stringsAsFactors=F)
ndvi.dat <- read.csv('output/tables/covs/pndvi.scale.csv', stringsAsFactors=F)
#burn.dat <- read.csv('output/tables/covs/burn.csv', stringsAsFactors=F)

row.names(bkt.dat) <- bkt.dat$PopulationName
row.names(lflow.dat) <- lflow.dat$PopulationName
row.names(ndvi.dat) <- ndvi.dat$PopulationName
#row.names(burn.dat) <- burn.dat$PopulationName

# trim to pops with data for all covs
dats <- list(templag.dat, hflowlag.dat, lflowlag.dat, extent.dat, bkt.dat, lflow.dat, ndvi.dat, add.dat, rem.dat)#, burn.dat)
pops <- c()

for (dat in dats){
  if (length(pops)==0) {
    pops <- dat$PopulationName
  } else {
    pops <- intersect(pops, dat$PopulationName)
  }
}

templag.dat <- templag.dat[pops,]
hflowlag.dat <- hflowlag.dat[pops,]
lflowlag.dat <- lflowlag.dat[pops,]
extent.dat <- extent.dat[pops,]
bkt.dat <- bkt.dat[pops,]
lflow.dat <- lflow.dat[pops,]
ndvi.dat <- ndvi.dat[pops,]
#burn.dat <- burn.dat[pops,]
add.dat <- add.dat[pops,]
rem.dat <- rem.dat[pops,]


# trim to years with data for all covs
dats <- list(templag.dat, hflowlag.dat, lflowlag.dat, bkt.dat, lflow.dat, ndvi.dat)#, burn.dat)
years <- c()

for (dat in dats){
  if (length(years)==0) {
    years <- as.numeric(gsub('X', '', names(dat[,2:ncol(dat)]), fixed=T))
  } else {
    years <- intersect(years, as.numeric(gsub('X', '', names(dat[,2:ncol(dat)]), fixed=T)))
  }
}

templag.dat <- templag.dat[,c(1, which(names(templag.dat) %in% paste('X',years,sep='')))]
hflowlag.dat <- hflowlag.dat[,c(1, which(names(hflowlag.dat) %in% paste('X',years,sep='')))]
lflowlag.dat <- lflowlag.dat[,c(1, which(names(lflowlag.dat) %in% paste('X',years,sep='')))]
bkt.dat <- bkt.dat[,c(1, which(names(bkt.dat) %in% paste('X',years,sep='')))]
lflow.dat <- lflow.dat[,c(1, which(names(lflow.dat) %in% paste('X',years,sep='')))]
ndvi.dat <- ndvi.dat[,c(1, which(names(ndvi.dat) %in% paste('X',years,sep='')))]
#burn.dat <- burn.dat[,c(1, which(names(burn.dat) %in% paste('X',years,sep='')))]
add.dat <- add.dat[,c(1, which(names(add.dat) %in% paste('X',years,sep='')))]
rem.dat <- rem.dat[,c(1, which(names(rem.dat) %in% paste('X',years,sep='')))]

# save trimmed covs
write.csv(templag.dat, file='output/tables/trimcovs/templag.dat.csv', row.names=F)
write.csv(hflowlag.dat, file='output/tables/trimcovs/hflowlag.dat.csv', row.names=F)
write.csv(lflowlag.dat, file='output/tables/trimcovs/lflowlag.dat.csv', row.names=F)
write.csv(extent.dat, file='output/tables/trimcovs/extent.dat.csv', row.names=F)
write.csv(bkt.dat, file='output/tables/trimcovs/bkt.dat.csv', row.names=F)
write.csv(lflow.dat, file='output/tables/trimcovs/lflow.dat.csv', row.names=F)
write.csv(ndvi.dat, file='output/tables/trimcovs/ndvi.dat.csv', row.names=F)
#write.csv(burn.dat, file='output/tables/trimcovs/burn.dat.csv', row.names=F)
write.csv(add.dat, file='output/tables/trimcovs/add.dat.csv', row.names=F)
write.csv(rem.dat, file='output/tables/trimcovs/rem.dat.csv', row.names=F)

# Cleanup
rm(list=ls()[-which(ls() %in% c('d','zm'))])

