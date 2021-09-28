#rm(list=ls())
#cat("\014") 
gc()

library(runjags)

# # Load model
load('output/objects/d.R')

# Load data
extent.dat <- read.csv('output/tables/trimcovs/extent.dat.csv', stringsAsFactors=F)
templag.dat <- read.csv('output/tables/trimcovs/templag.dat.csv', stringsAsFactors=F)
hflowlag.dat <- read.csv('output/tables/trimcovs/hflowlag.dat.csv', stringsAsFactors=F)
lflowlag.dat <- read.csv('output/tables/trimcovs/lflowlag.dat.csv', stringsAsFactors=F)
bkt.dat <- read.csv('output/tables/trimcovs/bkt.dat.csv', stringsAsFactors=F)
lflow.dat <- read.csv('output/tables/trimcovs/lflow.dat.csv', stringsAsFactors=F)
ndvi.dat <- read.csv('output/tables/trimcovs/ndvi.dat.csv', stringsAsFactors=F)
burn.dat <- read.csv('output/tables/trimcovs/burn.dat.csv', stringsAsFactors=F)
add.dat <- read.csv('output/tables/trimcovs/add.dat.csv', stringsAsFactors=F)
rem.dat <- read.csv('output/tables/trimcovs/rem.dat.csv', stringsAsFactors=F)

row.names(extent.dat) <- extent.dat$PopulationName
row.names(templag.dat) <- templag.dat$PopulationName
row.names(hflowlag.dat) <- hflowlag.dat$PopulationName
row.names(lflowlag.dat) <- lflowlag.dat$PopulationName
row.names(bkt.dat) <- bkt.dat$PopulationName
row.names(lflow.dat) <- lflow.dat$PopulationName
row.names(ndvi.dat) <- ndvi.dat$PopulationName
row.names(burn.dat) <- burn.dat$PopulationName
row.names(add.dat) <- add.dat$PopulationName
row.names(rem.dat) <- rem.dat$PopulationName

load('output/objects/popnames.R')
load('output/objects/jags.dat.R')
attach(jags.dat, warn.conflicts=F)

# Load functions
source('scripts/10a fun-stpvmforecast.R')
source('scripts/10b fun-plotsim.R')

# Run simulation
outdir <- paste(getwd(),'/output/forecast/', sep='')
dir.create(outdir)

popcompare <- c('EMR&MRBC&CC&QC&Cutt&Short&WillBas','WMR&CmpDraw&Gaws','MRBC_up','Frazer','Gance_RdCny_Warm',
                'Foreman','Abel','Indian','MohawkCanyon','Humboldt_NF&ColeCyn','T_Creek&Draw','Threemile','Tierney',
                'ToeJam','Boulder_4th')
popcompare <- sort(popcompare)

for (popname in popnames){
  popnum <- which(popnames==popname)
  covs <- list(extent = extent.dat[popname, 2:ncol(extent.dat)], 
               templag = templag.dat[popname, 2:ncol(templag.dat)], 
               hflowlag = hflowlag.dat[popname, 2:ncol(hflowlag.dat)], 
               lflowlag = lflowlag.dat[popname, 2:ncol(lflowlag.dat)], 
               bkt = bkt.dat[popname, 2:ncol(bkt.dat)], 
               lflow = lflow.dat[popname, 2:ncol(lflow.dat)], 
               ndvi = ndvi.dat[popname, 2:ncol(ndvi.dat)],
               burn = burn.dat[popname, 2:ncol(burn.dat)],
               reintro = add.dat[popname, 2:ncol(add.dat)] - rem.dat[popname, 2:ncol(rem.dat)]
               )
  
  extinct.row <- data.frame(row.names=c(popname), PopulationName=popname, ExtinctRisk=NA, ERYear=NA, ER2015=NA, ERForecastYear=NA)
  if (file.exists('output/tables/forecast.extinct.dat.csv')) {
	  extinct.dat <- read.csv('output/tables/forecast.extinct.dat.csv', stringsAsFactors=F)
	  row.names(extinct.dat) <- extinct.dat$PopulationName
	  if (!popname %in% extinct.dat$PopulationName) {
	    extinct.dat <- rbind(extinct.dat, extinct.row, stringsAsFactors=F)
	  }
  } else {
    extinct.dat <- extinct.row
  }
  stpvmforecast(forecastyear=2100, thin=100, nsim=50, popname=popname, popnum=popnum, 
                covs=covs, model=d, outdir=outdir, jags.dat=jags.dat, extinct.dat=extinct.dat,
                cut = list( templag=c(0,1), hflowlag=c(0,1), bkt=c(0,1), ndvi=c(0,1) ), 
                const = list( reintro=0, templag=NA, hflowlag=NA, bkt=NA, ndvi=NA, burn=NA )
                )
}

# Cleanup
rm(list=ls()[-which(ls() %in% c('d','zm'))])
