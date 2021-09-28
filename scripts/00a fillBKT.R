rm(list=ls())
gc()
#cat("\014")

dir.create('output')
dir.create('output/plots')
dir.create('output/tables')

#bkt <- read.csv('input/BKT_e_mod.20170313.csv', stringsAsFactors=F)
bkt <- read.csv('input/BBRC_e_mod.20170518.csv', stringsAsFactors=F)
row.names(bkt) <- bkt$PopulationName
bkt.orig <- bkt

col1 <- which(names(bkt)=='X1983')

year0 <- as.numeric(gsub('X','',names(bkt)[col1])) - 1

erad <- read.csv('input/bkt.eradication.20160321.csv', stringsAsFactors=F)
erad <- erad[erad$Method==9,]

# # Delete false negatives (zeros that are after a non-zero value).  Assumes there are no natural BKT extinctions
# for (i in 1:nrow(bkt)){
#   j.extant1 <- which(bkt[i,] > 0)[2]
#   if (!is.na(j.extant1)) {
#     j.extant <- j.extant1:ncol(bkt)
#     j.zero <- which(bkt[i, ]==0)
#     
#     bkt[i, j.zero[j.zero %in% j.extant]] <- NA  
#   }
# }

# Fill in zeros after Rotenone treatments.  Assumes Rotenone = BKT local extinction
for (i in 1:nrow(erad)){
  pop <- erad[i,'PopulationName']
  year <- erad[i, 'Year']
  j.rote <- year-year0+1
  
  j.post <- (j.rote+1):ncol(bkt)
  
  for(j in j.post){
    if (is.na(bkt[pop,j])) {
      bkt[pop,j] <- 0
    } else {if (bkt[pop,j] > 0) break}
  }
}

# Fill in a population's time series with zeros if 1) all observations = 0, or 2) all observations = NA
for (i in 1:nrow(bkt)){
  js <- col1:ncol(bkt)
  if (sum(!is.na(bkt[i,js]))==0){ 
    bkt[i,js] <- 0
  } else {
    if (max(bkt[i,js], na.rm=T)==0){
      bkt[i,js] <- 0
    }
  }
}

bkt.test <- bkt
bkt <- bkt.test

# Fill remaining NAs with negative exponential weighted average
BKT <- bkt

lam <- 0.5
j1 <- col1
j2 <- ncol(BKT)

for (i in 1:nrow(bkt)) {
  j.dat <- which(!is.na(BKT[i,]))[-1]
  #for (j in which(is.na(BKT[i,]))) {
  for (j in col1:ncol(bkt)) {
    
    # negative exponential kernal function to fill in missing data
    jdist <- abs(j - (j1:j2))
    
    weights <- lam * exp(-lam * jdist)
    values <- as.numeric(BKT[i,j1:j2])
    weights <- weights / sum(weights[!is.na(values)])
    bkt[i,j] <- sum(values * weights, na.rm=T)
    
  }
}

# Fill in zeros after Rotenone treatments.  Assumes Rotenone = BKT local extinction
for (i in 1:nrow(erad)){
  pop <- erad[i,'PopulationName']
  year <- erad[i, 'Year']
  j.rote <- year - year0 + 1
  
  j.post <- (j.rote+1):ncol(bkt)
  
  for(j in j.post){
    if (is.na(bkt.orig[pop,j])) {
      bkt[pop,j] <- 0
    } else {if (bkt.orig[pop,j] > 0) break}
  }
}

write.csv(bkt, 'output/tables/bkt.fill.csv', row.names=F)


## ========= ##
## BKT Plots
## ========= ##
pdf('output/plots/bkt.pdf')
for(i in 1:nrow(bkt)){
  years <- as.numeric(gsub('X','',names(bkt)[col1:ncol(bkt)]))
  
  y = bkt[i,col1:ncol(bkt)] * 1000
  y.orig = bkt.orig[i,col1:ncol(bkt.orig)] * 1000
  
  plot(x=years, y=y, type='l', ylim=c(0,max(y, y.orig, na.rm=T)*1.05),
       xlab='Year', ylab='BKT / km',
       main=paste('Pop: ',i, ' (',bkt[i,'PopulationName'],')', sep=''))
  
  points(x=years, y=y.orig, pch=16)
  
  x <- erad[erad$PopulationName==bkt[i,'PopulationName'],'Year']
  points(x=x, y=rep(0, length(x)), pch=13, cex=2)
}
dev.off()

