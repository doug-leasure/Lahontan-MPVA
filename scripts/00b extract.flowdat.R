library(foreign)

rm(list=ls())
gc()
#cat("\014") 

dir.create('output')
dir.create('output/tables')

metapop <- read.csv('input/LCT_MetaPop_Pop.20170313.csv', stringsAsFactors=F)

comids <- metapop$COMID

years <- 1984:2015

PN <- TRUE

dat1.full <- read.dbf('input/GB16_flow_met_acc.20170302.dbf')

dat1 <- dat1.full[dat1.full$COMIDV2 %in% comids & dat1.full$year %in% years, ]
rm(dat1.full)
gc()

if(PN){
  dat2.full <- read.dbf('input/PN17_flow_met_acc_trim.20170314.dbf')
  dat2 <- dat2.full[dat2.full$COMIDV2 %in% comids & dat2.full$year %in% years, ]
  rm(dat2.full)
  gc()
}

cols <- c("COMIDV2", "year", "MA", "MS", "Qmax", "Qmax3")
if(PN){ 
  dat <- rbind(dat1[,cols], dat2[,cols]) 
  dat <- dat[!duplicated(dat),]
} else {
  dat <- dat1[,cols]
  dat <- dat[!duplicated(dat),]
}

na.count <- 0
for (comid in comids){
  for (year in years){
    na.count <- na.count + (sum(dat$COMIDV2==comid & dat$year==year) == 0)
  }
}
na.count

write.csv(dat, file='output/tables/flow.csv', row.names=F)
