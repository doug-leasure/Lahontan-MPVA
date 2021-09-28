rm(list=ls())
gc()
#cat("\014") 
options(scipen=5) # stop R from converting to scientific notation so often

#setwd('D:/RESEARCH/2015 NASA Trout Project/wd/1.002')

dir.create('output')
dir.create('output/objects')
dir.create('output/tables')

# Read fish data
dat <- read.csv('output/tables/lct.reduced.dat.csv', header=T, stringsAsFactors=F)
dat.orig <- dat

## Setup dat years
firstyear <- 1984
lastyear <- 2015
nyears <- lastyear - firstyear + 1

#dat <- dat[dat$Year %in% c((firstyear+1):lastyear),]
dat$year <- dat$Year - (firstyear - 1)

# Assign numeric population ID
load('output/objects/popnames.R')

for(i in 1:nrow(dat)){
  dat[i, 'pop'] <- which(popnames==dat[i,'PopulationName'])
}
npops <- length(unique(dat$pop))

# Prop: Calculate proportion of population's range reprsented by each "site" (i.e. reach of known length)
dat$prop <- round(dat$SiteLength_m/(dat$Pop_Extent_km*1000),5)

# Site drainage areas
drain.dat <- read.csv('output/tables/covs/drain.csv')

for (row in row.names(dat)){
  SiteID <- dat[row,'SiteID']
  dat[row,'drain'] <- drain.dat[drain.dat$SiteID==SiteID,'drain.scale']
}

## Pack data ##

# Determine max number of years with data for a population
ncol.t_ti <- 0
for (i in unique(dat$pop)) {
  ncol.new <- length(unique(dat[dat$pop==i, 'year']))
  if (ncol.new > ncol.t_ti) ncol.t_ti <- ncol.new
}

# Matrices to convert time index (ti) to time (t = years since 1984) and vice versa.
t_ti <- matrix(NA, nrow=length(unique(dat$pop)), ncol=ncol.t_ti)
ti_t <- matrix(NA, nrow=length(unique(dat$pop)), ncol=nyears)

# Vector to store number of years with data for each population
nti <- c()

# Matrix to store number of sites with data for each population x year
nji <- matrix(NA, nrow=npops, ncol=ncol(t_ti))

# Relate t to ti,  and add new columns to data for ti and ji
for (i in unique(dat$pop)) {
  t.list <- unique(dat[dat$pop==i, 'year'])
  nti[i] <- length(t.list)
  t_ti[i,] <- c(t.list, rep(NA, ncol.t_ti - nti[i]))
  
  for (t in t.list) {
    ti <- which(t.list==t)
    dat[dat$pop==i & dat$year==t, 't'] <- ti
    ti_t[i,t] <- ti
    
    j.list <- unique(dat[dat$pop==i & dat$year==t, 'SiteID'])
    nji[i,ti] <- length(j.list)
    
    # Assign numeric site IDs (consecutive from 1 for each pop x year)
    for (j in j.list) {
      dat[dat$pop==i & dat$year==t & dat$SiteID==j, 'j'] <- which(j.list==j)
    }
  }
}

nt <- nti
nj <- nji

rm(i, t, j, t.list, j.list, ncol.new, ncol.t_ti, nti, nji)

# Write csv: original columns + new columns
dat.orig <- dat
write.csv(dat.orig, file='output/tables/lct.dat.origcols.csv', row.names=F)



# Keep only columns needed for model
dat <- dat[,c('pop','t','j','PassNumber','N_Age1pl','prop','SiteLength_m','drain')]
names(dat) <- c('pop','year','site','pass','y','prop','sitelength','drain')
dat <- dat[order(dat$pop, dat$year, dat$site, dat$pass),]
write.csv(dat, file='output/tables/lct.dat.csv', row.names=F)

dat <- read.csv(file='output/tables/lct.dat.csv')

# site-specific drainage area
drain <- array(NA, dim=c(npops, max(nt,na.rm=T), max(nj,na.rm=T)))
for(i in 1:npops){
  for(t in 1:nt[i]){
    for(j in 1:nj[i,t]){
      drain[i,t,j] <- unique(dat[dat$pop==i & dat$year==t & dat$site==j, 'drain'])
    }
  }
}

# npasses
npasses <- array(NA, dim=c(npops, max(nt, na.rm=T), max(nj, na.rm=T)))
for (i in 1:npops){
  for (t in 1:nt[i]){
    for (j in 1:nj[i,t]){
      npasses[i,t,j] <- max(dat[dat$pop==i & dat$year==t & dat$site==j, 'pass'], na.rm=T)
    }
  }
}

# y, Yj, Yi, prop
y <- array(NA, dim=c(npops, max(nt, na.rm=T), max(nj, na.rm=T), max(npasses, na.rm=T)))
for (row in row.names(dat)) {
  x <- dat[row,]
  y[x$pop, x$year, x$site, x$pass] <- x$y
}

Yj <- array(NA, dim=c(npops, max(nt, na.rm=T), max(nj, na.rm=T)))
Yi <- matrix(NA, nrow=npops, ncol=max(nt, na.rm=T))
#jd <- array(NA, dim=c(npops, max(nt, na.rm=T), max(nj, na.rm=T)))
for (i in 1:npops[1]){
  for (t in 1:nt[i]){
    for (j in 1:nj[i,t]){
      Yj[i,t,j] <- sum(y[i,t,j,], na.rm=T)
      
    #   jdx <- unique(dat.orig[dat.orig$pop==i & dat.orig$t==t & dat.orig$j==j, 'Age1.Size_Definition'])
    #   if(is.na(jdx)){ jd[i,t,j] <- 0 
    #   } else { jd[i,t,j] <- as.numeric(jdx==152) }
    }
    Yi[i,t] <- sum(y[i,t,,], na.rm=T)
  }
}

prop <- array(NA, dim=c(npops, max(nt, na.rm=T), max(nj, na.rm=T)+1))
sitelength <- array(NA, dim=c(npops, max(nt, na.rm=T), max(nj, na.rm=T)+1))

for (i in 1:npops[1]){
  for (t in 1:nt[i]){
    for (j in 1:nj[i,t]){
      prop[i,t,j] <- unique(dat[dat$pop==i & dat$year==t & dat$site==j,'prop'])
      sitelength[i,t,j] <- unique(dat[dat$pop==i & dat$year==t & dat$site==j,'sitelength'])
    }
    prop[i,t,nj[i,t]+1] <- 1 - sum(prop[i,t,1:nj[i,t]])
  }
}

# Save Yi matrix index to actual years
Yi.yr <- matrix(NA, nrow=npops, ncol=lastyear-firstyear+1)
for (i in 1:npops[1]){
  for (ti in 1:nt[i]){
    Yi.yr[i,t_ti[i,ti]] <- Yi[i,ti]
  }
}
Yi.yr <- data.frame(Yi.yr)
names(Yi.yr) <- as.character(firstyear:lastyear)

comids <- c()
names <- c()
for (pop in unique(dat.orig$pop)){
  comids <- c(comids, unique(dat.orig[dat.orig$pop==pop,'COMID']))
  names <- c(names, unique(dat.orig[dat.orig$pop==pop,'PopulationName']))
}

Yi.yr <- data.frame(PopulationName=names, COMID=comids, pop=unique(dat.orig$pop), Yi.yr)

write.csv(Yi.yr, 'output/tables/Yi.yr.csv', row.names=F)





#### Covariates ####
load('output/objects/popnames.R')

metapop <- read.csv('input/LCT_MetaPop_Pop.20170313.csv', stringsAsFactors=F)

# Make a map of PopulationName to COMID to pop#
map <- metapop[metapop$PopulationName %in% popnames, c('PopulationName', 'COMID', 'PopID')]
for (popname in map$PopulationName){
  map[map$PopulationName==popname, 'pop'] <- which(popnames == popname)
}
map <- map[order(map$pop),]

# extent
extent <- read.csv('output/tables/covs/extent.csv', stringsAsFactors=F)
row.names(extent) <- extent$PopulationName
extent <- extent[popnames, 'extent']

# Read data:  temp, ndvi, ndviC, ndwi
temp.dat <- read.csv('output/tables/covs/temp.scale.csv', stringsAsFactors=F)
templag.dat <- read.csv('output/tables/covs/templag.scale.csv', stringsAsFactors=F)
gndvi.dat <- read.csv('output/tables/covs/gndvi.scale.csv', stringsAsFactors=F)
gndvilag.dat <- read.csv('output/tables/covs/gndvilag.scale.csv', stringsAsFactors=F)
pndvi.dat <- read.csv('output/tables/covs/pndvi.scale.csv', stringsAsFactors=F)
pndvilag.dat <- read.csv('output/tables/covs/pndvilag.scale.csv', stringsAsFactors=F)
ms.dat <- read.csv('output/tables/covs/ms.scale.csv', stringsAsFactors=F)
mslag.dat <- read.csv('output/tables/covs/mslag.scale.csv', stringsAsFactors=F)
gms.dat <- read.csv('output/tables/covs/gms.scale.csv', stringsAsFactors=F)
qmax.dat <- read.csv('output/tables/covs/qmax.scale.csv', stringsAsFactors=F)
qmaxlag.dat <- read.csv('output/tables/covs/qmaxlag.scale.csv', stringsAsFactors=F)
qmax3.dat <- read.csv('output/tables/covs/qmax3.scale.csv', stringsAsFactors=F)
qmax3lag.dat <- read.csv('output/tables/covs/qmax3lag.scale.csv', stringsAsFactors=F)

# Fill in data with pop#
for (popname in popnames){
  pop <- map[map$PopulationName==popname, 'pop']
  
  temp.dat[temp.dat$PopulationName==popname,'pop'] <- pop
  templag.dat[templag.dat$PopulationName==popname,'pop'] <- pop
  gndvi.dat[gndvi.dat$PopulationName==popname,'pop'] <- pop
  gndvilag.dat[gndvilag.dat$PopulationName==popname,'pop'] <- pop
  pndvi.dat[pndvi.dat$PopulationName==popname,'pop'] <- pop
  pndvilag.dat[pndvilag.dat$PopulationName==popname,'pop'] <- pop
  ms.dat[ms.dat$PopulationName==popname,'pop'] <- pop
  mslag.dat[mslag.dat$PopulationName==popname,'pop'] <- pop
  gms.dat[gms.dat$PopulationName==popname,'pop'] <- pop
  qmax.dat[qmax.dat$PopulationName==popname,'pop'] <- pop
  qmaxlag.dat[qmaxlag.dat$PopulationName==popname,'pop'] <- pop
  qmax3.dat[qmax3.dat$PopulationName==popname,'pop'] <- pop
  qmax3lag.dat[qmax3lag.dat$PopulationName==popname,'pop'] <- pop
}


# Remove rows with pops not in analysis, remove columns with years not in analysis, re-order by pop#, and convert to matrix
# For remote sensed covs, fill in NAs of early years with population means: ndvi, ndivC, and ndwi
temp <- temp.dat[!is.na(temp.dat$pop), ]
temp <- as.matrix(temp[order(temp$pop), paste('X',firstyear:lastyear,sep='')])

templag <- templag.dat[!is.na(templag.dat$pop), ]
templag <- as.matrix(templag[order(templag$pop), paste('X',firstyear:lastyear,sep='')])

gndvi <- gndvi.dat[!is.na(gndvi.dat$pop), ]
gndvi <- as.matrix(gndvi[order(gndvi$pop), paste('X',1985:lastyear,sep='')])
gndvi <- cbind(as.vector(apply(gndvi, 1, mean)), gndvi)

gndvilag <- gndvilag.dat[!is.na(gndvilag.dat$pop), ]
gndvilag <- as.matrix(gndvilag[order(gndvilag$pop), paste('X',1986:lastyear,sep='')])
gndvilag <- cbind(as.vector(apply(gndvilag, 1, mean)), as.vector(apply(gndvilag, 1, mean)), gndvilag)

pndvi <- pndvi.dat[!is.na(pndvi.dat$pop), ]
pndvi <- as.matrix(pndvi[order(pndvi$pop), paste('X',1985:lastyear,sep='')])
pndvi <- cbind(as.vector(apply(pndvi, 1, mean)), pndvi)

pndvilag <- pndvilag.dat[!is.na(pndvilag.dat$pop), ]
pndvilag <- as.matrix(pndvilag[order(pndvilag$pop), paste('X',1986:lastyear,sep='')])
pndvilag <- cbind(as.vector(apply(pndvilag, 1, mean)), as.vector(apply(pndvilag, 1, mean)), pndvilag)

ms <- ms.dat[!is.na(ms.dat$pop), ]
ms <- as.matrix(ms[order(ms$pop), paste('X',firstyear:lastyear,sep='')])

mslag <- mslag.dat[!is.na(mslag.dat$pop), ]
mslag <- as.matrix(mslag[order(mslag$pop), paste('X',firstyear:lastyear+1,sep='')])

gms <- gms.dat[!is.na(gms.dat$pop), ]
gms <- as.matrix(gms[order(gms$pop), paste('X',firstyear:lastyear,sep='')])

qmax <- qmax.dat[!is.na(qmax.dat$pop), ]
qmax <- as.matrix(qmax[order(qmax$pop), paste('X',firstyear:lastyear,sep='')])

qmaxlag <- qmaxlag.dat[!is.na(qmaxlag.dat$pop), ]
qmaxlag <- as.matrix(qmaxlag[order(qmaxlag$pop), paste('X',firstyear:lastyear+1,sep='')])

qmax3 <- qmax3.dat[!is.na(qmax3.dat$pop), ]
qmax3 <- as.matrix(qmax3[order(qmax3$pop), paste('X',firstyear:lastyear,sep='')])

qmax3lag <- qmax3lag.dat[!is.na(qmax3lag.dat$pop), ]
qmax3lag <- as.matrix(qmax3lag[order(qmax3lag$pop), paste('X',firstyear:lastyear+1,sep='')])


### BKT ###
bkt.dat <- read.csv('output/tables/covs/bkt.scale.csv')
names.bkt.dat <- names(bkt.dat)
for (popname in map$PopulationName){
  bkt.dat[bkt.dat$PopulationName==popname, 'pop'] <- map[map$PopulationName==popname,'pop']
}
bkt.dat <- bkt.dat[!is.na(bkt.dat$pop),]
bkt.dat <- bkt.dat[order(bkt.dat$pop),]
bkt.dat <- bkt.dat[,sort(names.bkt.dat[2:length(names.bkt.dat)])]
#bkt.dat <- bkt.dat * 10 # convert from BKT/m to BKT/10m
bkt <- as.matrix(bkt.dat[,paste('X',firstyear:lastyear,sep='')])

bktlag.dat <- bkt.dat
names(bktlag.dat) <- gsub('X', '', names(bktlag.dat))
names(bktlag.dat)[1:ncol(bktlag.dat)] <- as.character(as.numeric(names(bktlag.dat)[1:dim(bktlag.dat)[2]])+1)
bktlag <- as.matrix(bktlag.dat[,as.character(firstyear:lastyear)])

# # Burn proportion
# burn.dat <- read.dbf('input/burnprop.20170313.dbf')
# for (popname in popnames){
#   popid <- map[map$PopulationName==popname, 'PopID']
#   pop <- map[map$PopulationName==popname, 'pop']
#   burn.dat[burn.dat$PopID==popid,'pop'] <- pop
# }
# burn.dat <- burn.dat[!is.na(burn.dat$pop), ]
# burn.dat <- burn.dat[order(burn.dat$pop),]
# 
# names(burn.dat) <- gsub('wf', '', names(burn.dat))
# 
# burn <- as.matrix(burn.dat[,as.character(firstyear:lastyear)])
# write.csv(cbind(data.frame(PopulationName=popnames), burn), file='output/tables/covs/burn.csv', row.names=F)
# 
# burnlag <- cbind(matrix(0, ncol=1, nrow=nrow(burn)), burn[,1:(ncol(burn)-1)])
# colnames(burnlag) <- colnames(burn)
# write.csv(cbind(data.frame(PopulationName=popnames), burnlag), file='output/tables/covs/burnlag.csv', row.names=F)


# Addition/Removal Data
add.dat <- read.csv('output/tables/covs/add.csv')
row.names(add.dat) <- add.dat$PopulationName
add <- as.matrix(add.dat[popnames,paste('X',as.character(firstyear:lastyear),sep='')])

rem.dat <- read.csv('output/tables/covs/rem.csv')
row.names(rem.dat) <- rem.dat$PopulationName
rem <- as.matrix(rem.dat[popnames,paste('X',as.character(firstyear:lastyear),sep='')])

reintro <- add - rem


# ####  TEST  ####
# j.multi <- array(NA, dim=dim(npasses))
# nj.multi <- array(NA, dim=dim(nj))
# 
# for (i in 1:npops){
#   for (t in 1:nt[i]){
#     count.multi <- 0
#     for (j in 1:nj[i,t]){
#       if (npasses[i,t,j]>1) {
#         count.multi <- count.multi + 1
#         j.multi[i, t, count.multi] <- j
#       }
#     }
#     nj.multi[i,t] <- count.multi
#   }
# }
# 
# for (i in dim(j.multi)[3]:1){
#   if (sum(!is.na(j.multi[,,i])) == 0){
#     j.multi <- j.multi[,,-i]
#   }
# }
# 
# nj.tog <- nj.multi==0
# nj.tog <- matrix(as.numeric(nj.tog), ncol=ncol(nj.tog), nrow=nrow(nj.tog))
# ################


elev <- read.csv('input/elevation.20170615.csv')
for(row in row.names(elev)){
  popid <- as.numeric(elev[row, 'PopID'])
  if(popid %in% map$PopID){
    popname <- map[map$PopID==popid, 'PopulationName']
    elev[row, 'PopulationName'] <- popname  
  }
}
elev <- elev[!is.na(elev$PopulationName),]
row.names(elev) <- elev$PopulationName
elev <- elev[popnames,'elev_m']
elev <- as.numeric(elev > 2000) + 1


################
## EXPLORE NDVI EFFECT
# effort <- matrix(nrow=nrow(Yi), ncol=ncol(Yi))
# for(i in 1:npops){
#   for(t in 1:nt[i]){
#     effort[i,t] <- sum(prop[i,t,1:nj[i,t]]*extent[i]*1000*npasses[i,t,1:nj[i,t]], na.rm=T)
#   }
# }
# cpue <- Yi / effort
# ndvi <- pndvi
# plot(NA, 
#      ylab='CPUE', ylim=c(0,max(cpue,na.rm=T)),
#      xlab='NDVI',  xlim=c(min(ndvi), max(ndvi)))
# for(i in 1:npops){
#   for(t in 1:nt[i]){
#     if(elev[i]==1) points(y=cpue[i,t], x=ndvi[i,t_ti[i,t]], col='red')
#     if(elev[i]==2) points(y=cpue[i,t], x=ndvi[i,t_ti[i,t]], col='blue')
#   }
# }
# legend('topleft', legend=c('High elevation', 'Low elevation'), col=c('blue','red'), pch=1)
#################


### jags.dat ###
jags.dat <- list(y = y,
                 Yj = Yj,
                 # Yj2 = Yj,
                 Yi = Yi,
                 t_ti = t_ti,
                 # prop = prop,
                 sitelength = sitelength/1000,
                 npops = npops,
                 nt = nt,
                 nj = nj,     
                 npasses = npasses,
                 # jmult = j.multi,
                 # njmult = nj.multi,
                 # njtog = nj.tog,
                 reintro=reintro,
                 templag=templag,
                 hflowlag=qmax3lag,
                 # lflow=ms,
                 extent=extent,
                 bkt=bkt,
                 ndvi=pndvi,
                 elev=elev,
                 drain=drain
                 )

save(jags.dat, file='output/objects/jags.dat.R')

# Write csv: final data
write.csv(dat, file='output/tables/lct.dat.csv', row.names=F)

# Cleanup
rm(list=ls())



