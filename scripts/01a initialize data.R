library(foreign)

rm(list=ls())
gc()
#cat("\014") 
options(scipen=5) # stop R from converting to scientific notation so often

#setwd('D:/RESEARCH/2015 NASA Trout Project/wd/1.002')

dir.create('output')
dir.create('output/objects')
dir.create('output/tables')
dir.create('output/tables/covs')

# Read original fish data
dat <- read.csv('input/LCT_data_20170518.csv', header=T, stringsAsFactors=F)

dat.orig <- dat

#### Make a map of PopulationName to COMID to extent
map <- read.csv('input/LCT_MetaPop_Pop.20170313.csv', stringsAsFactors=F)[,c('PopID','PopulationName','COMID','Pop_Extent_km')]
map <- map[!duplicated(map),]
map <- map[order(map$PopulationName),]
row.names(map) <- map$PopulationName

############# sTART:  REDUCE DATA #########################

## Remove problem records
dat <- dat[!is.na(dat$NumPasses),] # n=1
dat <- dat[!dat$NumPasses == 0,] # n=0
dat <- dat[!is.na(dat$PassNumber),] # n=25
dat <- dat[!is.na(dat$SiteLength_m),] # n=0

dat <- dat[!dat$SiteLength_m > 200,] # n=91

dat[!is.na(dat$Age1.Size_Definition) & dat$Age1.Size_Definition==152 & dat$N_Age0==0 & dat$N_Age1pl==0, 'Age1.Size_Definition'] <- NA
dat <- dat[is.na(dat$Age1.Size_Definition) | !dat$Age1.Size_Definition==152,]

dat <- dat[!dat$Lacustrine==1,] # n=49

dat <- dat[!dat$PopType=='Reconnect',]

## Remove pops that don't have covariate data
dat <- dat[dat$SiteID %in% read.dbf('input/SiteIDDrain.20170518.dbf')$SiteID, ] # n=41

row.names(map) <- map$PopulationName
pops.nocovs <- c()
pops.nocovs <- c(pops.nocovs, dat[!map[dat$PopulationName,'COMID'] %in% read.csv('output/tables/flow.csv',stringsAsFactors=F)$COMIDV2, 'PopulationName'])
pops.nocovs <- c(pops.nocovs, dat[!dat$PopulationName %in% read.csv('output/tables/bkt.fill.csv',stringsAsFactors=F)$PopulationName, 'PopulationName'])
pops.nocovs <- c(pops.nocovs, dat[!map[dat$PopulationName,'PopID'] %in% read.dbf('input/MeanAugC.20170313.dbf')$PopID, 'PopulationName'])

pops.nocovs <- c(pops.nocovs, dat[!map[dat$PopulationName,'PopID'] %in% read.csv('input/ndvi.20170313.csv',stringsAsFactors=F)$Pop_ID, 'PopulationName'])
#pops.nocovs <- c(pops.nocovs, dat[!map[dat$PopulationName,'PopID'] %in% read.csv('input/burn.20160509.csv',stringsAsFactors=F)$PopID, 'PopulationName'])
pops.nocovs <- unique(pops.nocovs)

# Pops Removed
pops.nocovs

# Pops remaining
unique(dat[,"PopulationName"])[!unique(dat[,"PopulationName"]) %in% pops.nocovs]

# Reduce dat to populations with covariates
dat <- dat[!dat$PopulationName %in% pops.nocovs, ]

#### Reduce pops for testing  ####

# dat <- dat[dat$PopulationName %in% c('Abel','Andorno','Pole_QRD','EMR&MRBC&CC&QC&Cutt&Short&WillBas','Washburn_up','Crowley'
# 'Coyote', 'Foreman_up','Frazer','Gance_RdCny_Warm','Humboldt_NF&ColeCyn', 'Indian',
# 'Maggie','MohawkCanyon','Threemile','Tierney','ToeJam','Cabin&DJ&Lye&Martin&Deep&RoundCorral&Alkali','Marysville','Pearl','KlondikeCan&Quinn_EF&Laca','Stewart&Stewart_NF&MF',
# 'Whitehorse_OR','Big&Big_NF', 'Boone', 'Brown', 'BuffaloCanyon'
# ),]

## Setup dat years
firstyear <- 1984
lastyear <- 2015
dat <- dat[dat$Year %in% c((firstyear+1):lastyear),]

# # Remove pops with only one year of data
# pops.1year <- c()
# for (pop in unique(dat$PopulationName)){
#   if (length(unique(dat[dat$PopulationName==pop, 'Year']))==1){
#     pops.1year <- c(pops.1year, pop)
#   }
# }
# dat <- dat[!dat$PopulationName %in% pops.1year, ]

# Re-order dat
dat <- dat[order(dat$PopulationName, dat$Year, dat$SiteID, dat$PassNumber), ]

print(paste('Rows Removed = ', nrow(dat.orig) - nrow(dat), ' (', round(100-(nrow(dat)/nrow(dat.orig)*100),1), '%)' , sep=''))
print(paste('Pops Removed = ', length(unique(dat.orig$PopID)) - length(unique(dat$PopID)), ' (', round(100-(length(unique(dat$PopID))/length(unique(dat.orig$PopID))*100),1), '%)' , sep=''))
unique(dat.orig$PopulationName)[!unique(dat.orig$PopulationName) %in% unique(dat$PopulationName)]

#### Save dat ####
write.csv(dat, file='output/tables/lct.reduced.dat.csv')

############### END: REDUCE DATA ###################



#### Save popnames object ####
popnames <- unique(dat$PopulationName)
save(popnames, file='output/objects/popnames.R')

allpopnames <- unique(map$PopulationName)
allpopnames <- allpopnames[allpopnames %in% read.csv('output/tables/bkt.fill.csv',stringsAsFactors=F)$PopulationName]
save(allpopnames, file='output/objects/allpopnames.R')

#### Initialize Covariates ###
row.names(map) <- map$PopID

## Additions/Removals
addrem <- read.csv('input/addrem.20170518.csv', stringsAsFactors=F)

popids <- list()
for (name in allpopnames) popids[name] <- map[which(map$PopulationName==name),'PopID']
  
rem <- add <- data.frame(PopulationName=allpopnames, PopID=NA)
for(i in 1:length(popnames)){
  rem[i,'PopID'] <- add[i,'PopID'] <- popids[[popnames[i]]]
}
rownames(rem) <- rownames(add) <- as.character(popids)
for (year in firstyear:lastyear) add[,as.character(year)] <- rem[,as.character(year)] <- 0

for(i in row.names(addrem)){
  if(addrem[i,'PopID'] %in% popids & addrem[i,'year'] %in% firstyear:lastyear){
    add[as.character(addrem[i,'PopID']), as.character(addrem[i,'year'])] <- addrem[i,'n']
    
    if (!is.na(addrem[i,'PopID_source']) & addrem[i,'PopID_source'] %in% popids) {
      rem[as.character(addrem[i,'PopID_source']), as.character(addrem[i,'year'])] <- addrem[i,'n']
    }
  }
}
write.csv(add[,-2], 'output/tables/covs/add.csv', row.names=F)
write.csv(rem[,-2], 'output/tables/covs/rem.csv', row.names=F)


# ndvi
ndvip <- read.csv('input/ndviP.20170313.csv', stringsAsFactors=F, check.names=F)
ndvip[,1] <- map[as.character(ndvip$Pop_ID), 'PopulationName']
names(ndvip)[1] <- 'PopulationName'
ndvip <- ndvip[,c('PopulationName',sort(names(ndvip[2:ncol(ndvip)])))]
row.names(ndvip) <- ndvip$PopulationName
ndvip <- ndvip[allpopnames,]

ndvi <- read.csv('input/ndvi.20170313.csv', stringsAsFactors=F, check.names=F)
ndvi[,1] <- map[as.character(ndvi$Pop_ID), 'PopulationName']
names(ndvi)[1] <- 'PopulationName'
ndvi <- ndvi[,c('PopulationName',sort(names(ndvi[2:ncol(ndvi)])))]
row.names(ndvi) <- ndvi$PopulationName
ndvi <- ndvi[allpopnames,]

# fill in data where proportion cloud cover was high
thresh <- 0.5
firstyear <- names(ndvi)[2]
lastyear <- names(ndvi)[ncol(ndvi)]
for (popname in ndvi$PopulationName) {
  for (year in names(ndvi[,2:ncol(ndvi)])) {
    if (ndvip[popname, year] > thresh) {
      if (year == firstyear) {
        ndvi[popname, year] <- ndvi[popname, as.character(as.numeric(year)+1)]
      } else if (year == lastyear) {
        ndvi[popname, year] <- ndvi[popname, as.character(as.numeric(year)-1)]
      } else{
        ndvi[popname, year] <- mean(as.matrix(ndvi[popname, as.character(c(as.numeric(year)-1, as.numeric(year)+1))]))
      }
    }
  }
}
write.csv(ndvi, 'output/tables/covs/ndvi.csv', row.names=F)

ndvilag <- ndvi
names(ndvilag)[2:ncol(ndvilag)] <- as.character(as.numeric(names(ndvilag)[2:dim(ndvilag)[2]])+1)
write.csv(ndvilag, 'output/tables/covs/ndvilag.csv', row.names=F)


# Temp
temp <- read.dbf('input/MeanAugC.20170313.dbf')

temp[,1] <- map[as.character(temp$PopID), 'PopulationName']
names(temp)[1] <- 'PopulationName'
row.names(temp) <- temp$PopulationName

names(temp)[2:dim(temp)[2]] <- gsub('nw', '', names(temp[2:dim(temp)[2]]))

temp <- temp[allpopnames, c('PopulationName', sort(names(temp)[2:ncol(temp)]))]

temp['WashONeal',2:ncol(temp)] <- temp['StoneHouse',2:ncol(temp)]
temp['MarySloan',2:ncol(temp)] <- temp['Happy',2:ncol(temp)]
temp['BuffaloCanyon',2:ncol(temp)] <- temp['HorseCanyon',2:ncol(temp)]
temp['Brown',2:ncol(temp)] <- temp['Pearl',2:ncol(temp)]

write.csv(temp, 'output/tables/covs/temp.csv', row.names=F)


templag <- temp
names(templag)[2:ncol(templag)] <- as.character(as.numeric(names(templag)[2:ncol(templag)])+1)
write.csv(templag, 'output/tables/covs/templag.csv', row.names=F)

rm(ndvip, firstyear, lastyear, popname, thresh, year)


#### Flow Data ####
flowdat <- read.csv('output/tables/flow.csv')
#mad.dat <- read.csv('MAD.csv')

comids <- map$COMID
years <- unique(flowdat$year)

qmax <- qmax3 <- ms <- data.frame(row.names=allpopnames, PopulationName=allpopnames)

ms$PopulationName <- as.character(ms$PopulationName)
qmax$PopulationName <- as.character(qmax$PopulationName)
qmax3$PopulationName <- as.character(qmax3$PopulationName)

for (popname in allpopnames){
  comid <- map[map$PopulationName==popname, 'COMID']
  for (year in years){
    ms[popname, as.character(year)] <- flowdat[flowdat$COMIDV2==comid & flowdat$year==year, 'MS']
    qmax[popname, as.character(year)] <- flowdat[flowdat$COMIDV2==comid & flowdat$year==year, 'Qmax']
    qmax3[popname, as.character(year)] <- flowdat[flowdat$COMIDV2==comid & flowdat$year==year, 'Qmax3']
  }
}

# # !!! LOG TRANSFORM !!! #
# ms[,-1] <- log(ms[,-1] + 1)
# qmax[,-1] <- log(qmax[,-1] + 1)
# qmax3[,-1] <- log(qmax3[,-1] + 1)

# write file
write.csv(ms, 'output/tables/covs/ms.csv', row.names=F)
write.csv(qmax, 'output/tables/covs/qmax.csv', row.names=F)
write.csv(qmax3, 'output/tables/covs/qmax3.csv', row.names=F)

# lag
mslag <- ms
names(mslag)[2:ncol(mslag)] <- as.character(as.numeric(names(mslag)[2:ncol(mslag)])+1)
write.csv(mslag, 'output/tables/covs/mslag.csv', row.names=F)

qmaxlag <- qmax
names(qmaxlag)[2:ncol(qmaxlag)] <- as.character(as.numeric(names(qmaxlag)[2:ncol(qmaxlag)])+1)
write.csv(qmaxlag, 'output/tables/covs/qmaxlag.csv', row.names=F)

qmax3lag <- qmax3
names(qmax3lag)[2:ncol(qmax3lag)] <- as.character(as.numeric(names(qmax3lag)[2:ncol(qmax3lag)])+1)
write.csv(qmax3lag, 'output/tables/covs/qmax3lag.csv', row.names=F)

rm(comid, comids, year, years, flowdat, popname)

### BKT
bktfill <- read.csv('output/tables/bkt.fill.csv', stringsAsFactors=F)
bkt <- bktfill[,-which(names(bktfill) %in% c('PopID', 'PopType'))]
write.csv(bkt, 'output/tables/covs/bkt.csv', row.names=F)

#### Get scale factors to scale covariates ####

## Scale globally (use only pops with data for the model)
load('output/objects/popnames.R')

# extent
extenti <- map$Pop_Extent_km[map$PopulationName %in% popnames]
scale.factors <- data.frame(cov='extent', mu=mean(extenti), sd=sd(extenti), PopulationName=NA)
scale.factors$cov <- as.character(scale.factors$cov)
scale.factors$PopulationName <- as.character(scale.factors$PopulationName)

# drain
drain.dat <- read.dbf('input/SiteIDDrain.20170518.dbf')
scale.factors <- rbind(scale.factors, c('drain', mean(drain.dat$CuDrainkm2), sd(drain.dat$CuDrainkm2), NA))

# temp
tempi <- as.vector(as.matrix(temp[temp$PopulationName %in% popnames,2:ncol(temp)]))
scale.factors <- rbind(scale.factors, c('temp', mean(tempi), sd(tempi), NA))

# gndvi
gndvii <- as.numeric(as.vector(as.matrix(ndvi[ndvi$PopulationName %in% popnames,2:ncol(ndvi)])))
scale.factors <- rbind(scale.factors, c('gndvi', mean(gndvii), sd(gndvii), NA))

# gms
gms <- as.numeric(as.vector(as.matrix(ms[ms$PopulationName %in% popnames,2:ncol(ms)])))
scale.factors <- rbind(scale.factors, c('gms', mean(gms), sd(gms), NA))

# bkt
bkt <- read.csv('output/tables/covs/bkt.csv', stringsAsFactors=F)
bkti <- as.numeric(as.vector(as.matrix(bkt[bkt$PopulationName %in% popnames, 2:ncol(bkt)])))
scale.factors <- rbind(scale.factors, c('bkt', mean(bkti), sd(bkti), NA))

rm(tempi, gndvii, bkti)

## Scale per pop

# pndvi
mu <- apply(as.matrix(ndvi[2:dim(ndvi)[2]]), 1, mean)
sigma <- apply(as.matrix(ndvi[2:dim(ndvi)[2]]), 1, sd)
for (popname in ndvi$PopulationName){
  scale.factors <- rbind(scale.factors, c('pndvi', mu[which(ndvi$PopulationName==popname)], sigma[which(ndvi$PopulationName==popname)], popname))
}
rm(mu, sigma)

# ms
mu <- apply(as.matrix(ms[2:dim(ms)[2]]), 1, mean)
sigma <- apply(as.matrix(ms[2:dim(ms)[2]]), 1, sd)
for (popname in ms$PopulationName){
  scale.factors <- rbind(scale.factors, c('ms', mu[which(ms$PopulationName==popname)], sigma[which(ms$PopulationName==popname)], popname))
}
rm(mu, sigma)

# qmax
mu <- apply(as.matrix(qmax[2:dim(qmax)[2]]), 1, mean)
sigma <- apply(as.matrix(qmax[2:dim(qmax)[2]]), 1, sd)
for (popname in qmax$PopulationName){
  scale.factors <- rbind(scale.factors, c('qmax', mu[which(qmax$PopulationName==popname)], sigma[which(qmax$PopulationName==popname)], popname))
}
rm(mu, sigma)

# qmax3
mu <- apply(as.matrix(qmax3[2:dim(qmax3)[2]]), 1, mean)
sigma <- apply(as.matrix(qmax3[2:dim(qmax3)[2]]), 1, sd)
for (popname in qmax3$PopulationName){
  scale.factors <- rbind(scale.factors, c('qmax3', mu[which(qmax3$PopulationName==popname)], sigma[which(qmax3$PopulationName==popname)], popname))
}
rm(mu, sigma)

## write scale factors to file
scale.factors$mu <- as.numeric(scale.factors$mu)
scale.factors$sd <- as.numeric(scale.factors$sd)
write.csv(scale.factors, file='output/tables/scale.factors.csv', row.names=F)

rm(popname)


#### Scale covariates

## Global scaling
covs <- c('temp', 'templag', 'gndvi', 'gndvilag', 'gms')
covs.scale <- c('temp', 'temp', 'gndvi', 'gndvi', 'gms')
dats <- list(temp, templag, ndvi, ndvilag, ms)

for (i in 1:length(covs)){
  dat <- dats[[i]][,2:ncol(dats[[i]])]
  dat.scale <- as.data.frame((dat - scale.factors[scale.factors$cov==covs.scale[i], 'mu']) / scale.factors[scale.factors$cov==covs.scale[i], 'sd'])
  scaled <- cbind(PopulationName=dats[[i]][,1], dat.scale)
  write.csv(scaled, file=paste('output/tables/covs/',covs[i],'.scale.csv',sep=''), row.names=F)
  rm(dat, dat.scale, scaled)
}
rm(i)

# Extent
extent <- map[,c('PopulationName', 'Pop_Extent_km')]
names(extent) <- c('PopulationName', 'extent')
extent.scale <- extent
extent.scale$extent <- (extent$extent - scale.factors[scale.factors$cov=='extent', 'mu']) / scale.factors[scale.factors$cov=='extent', 'sd']
write.csv(extent, 'output/tables/covs/extent.csv', row.names=F)
write.csv(extent.scale, 'output/tables/covs/extent.scale.csv', row.names=F)

# Drain
drain.dat$drain.scale <- (drain.dat$CuDrainkm2 - scale.factors[scale.factors$cov=='drain', 'mu']) / scale.factors[scale.factors$cov=='drain', 'sd']
write.csv(drain.dat, 'output/tables/covs/drain.csv', row.names=F)

# BKT
dat <- bkt[,2:ncol(bkt)]
dat.scale <- dat / scale.factors[scale.factors$cov=='bkt', 'sd']
bkt.scale <- cbind(PopulationName=bkt$PopulationName, dat.scale)
write.csv(bkt.scale, 'output/tables/covs/bkt.scale.csv', row.names=F)

## Per-pop scaling
covs <-       c('pndvi', 'pndvilag', 'ms', 'mslag', 'qmax', 'qmaxlag', 'qmax3', 'qmax3lag')
covs.scale <- c('pndvi', 'pndvi',    'ms', 'ms',    'qmax', 'qmax',    'qmax3', 'qmax3')
dats <-    list( ndvi,    ndvilag,    ms,   mslag,   qmax,   qmaxlag,   qmax3,   qmax3lag)

for (cov in 1:length(covs)){
  pops <- dats[[cov]]$PopulationName
  dat <- dats[[cov]][,2:ncol(dats[[cov]])]
  dat.scale <- dat
  
  for (pop in 1:nrow(dat.scale)){
    mu <- scale.factors[scale.factors$cov==covs.scale[cov] & scale.factors$PopulationName==pops[pop], 'mu'] 
    sigma <- scale.factors[scale.factors$cov==covs.scale[cov] & scale.factors$PopulationName==pops[pop], 'sd']
    dat.scale[pop,] <- as.data.frame((dat[pop,] - mu) / sigma)
    #rm(mu, sigma)
  }
  
  scaled <- cbind(PopulationName=row.names(dat.scale), dat.scale)
  write.csv(scaled, file=paste('output/tables/covs/',covs[cov],'.scale.csv',sep=''), row.names=F)
  #rm(dat, dat.scale, pops, scaled)
}

rm(covs, covs.scale, dats)

# Cleanup
rm(list=ls())
