stpvmforecast <- function(forecastyear, popname, popnum, covs, model, outdir, jags.dat, extinct.dat, thin=100, nsim=100, Nmax=1e9, cut=list(templag=c(0,1), hflowlag=c(0,1), bkt=c(0,1), ndvi=c(0,1)), const=list(reintro=0, templag=NA, hflowlag=NA, bkt=NA, ndvi=NA)){
  # forecastyear=2115; thin=10; nsim=50; popname=popname; popnum=popnum; covs=covs; model=d; outdir=outdir; jags.dat=jags.dat; extinct.dat=extinct.dat; Nmax=1e9;
  # cut = list( templag=c(0,1), hflowlag=c(0,1), bkt=c(0,1), ndvi=c(0,1) ); 
  # const = list( reintro=0, templag=NA, hflowlag=NA, bkt=NA, ndvi=NA, burn=NA );
  print(popname)
  
  firstyear <- 1984
  
  # initialize simulation
  jags.dat <- jags.dat[c('t_ti','nt','sex.ratio')]
  attach(jags.dat, warn.conflicts=FALSE)
  
  # Setup model
  model <- model[seq(1,nrow(model),by=thin),]
  
  N1 <- model[,paste('N[',popnum,',',t_ti[popnum,nt[popnum]],']',sep='')]
  sigmaR <- model[,paste('sigmaR[',popnum,']',sep='')]
  
  attach(model[,c('b0phi','b1phi','b2phi','b3phi','b0r','b1r','b2r','b3r')], warn.conflicts=FALSE)  
  
  # Setup covariates
  covs.sims <- list()
  covnames <- names(covs)[which(!names(covs)=='extent')]
  year1.sim <- (firstyear - 1) + t_ti[popnum,nt[popnum]]

  # Identify years with data for all covariates
  resamp.years <- names(covs[[covnames[1]]])
  for(cov in covnames) resamp.years <- intersect(resamp.years, names(covs[[cov]]))
  resamp.years <- sort(resamp.years)
  
  # Identify last year with covariates
  yearendcovs <- as.numeric(sub('X','', resamp.years[length(resamp.years)]))
  
  # Identify years when covariates are within constraints of cut[]
  for(cov in names(cut)){
    q <- quantile(as.numeric(covs[[cov]]), probs=cut[[cov]])
    resamp.years <- resamp.years[ resamp.years %in% names(covs[[cov]])[as.numeric(covs[[cov]]) >= q[1] & as.numeric(covs[[cov]]) <= q[2]] ]
  }

  # Setup covariate time series for each simulation
  for (sim in 1:nsim){
    t <- sample(x=resamp.years, size=forecastyear-yearendcovs, replace=T)
    covs.sims[[sim]] <- list()
    
    for (cov in covnames) {
      x <- covs[[cov]]
      x <- x[,which(names(x)==paste('X', year1.sim, sep='')):ncol(x)]
      
      if (cov %in% names(const)[which(!is.na(const))]) { x <- cbind(x, matrix(rep(const[[cov]],length=length(t)), nrow=1, ncol=length(t)))
      } else { x <- cbind(x, covs[[cov]][t]) }
      
      names(x) <- paste('X',year1.sim:forecastyear, sep='')
      
      covs.sims[[sim]][[cov]] <- x
    }
  }
  
  # Prepare simulation run
  nyears <- forecastyear - year1.sim + 1 #length(templag)
  nmcmc <- nrow(model)
  
  N <- array(NA, dim=c(nsim, nmcmc, nyears))
  N[,,1] <- N1
  
  phi <- r <- Rmean <- R <- muN <- array(NA, dim=c(nsim, nmcmc, nyears))
  
  # Run simulation
  for (sim in 1:nsim){
    
    attach(covs.sims[[sim]], warn.conflicts=FALSE)
    
    for (t in 2:nyears){
      N[sim,,t-1] <-  N[sim,,t-1] + reintro[[t-1]]
      N[sim,,t-1][N[sim,,t-1] < 0] <-  0
      
      Nfem <- rbinom(nmcmc, N[sim,,t-1], sex.ratio)
      
      phi[sim,,t] <- b0phi + b1phi*bkt[[t]] + b2phi*ndvi[[t]] + b3phi*burn[[t]]
      r[sim,,t] <- b0r + b1r*templag[[t]] + b2r*templag[[t]]^2 + b3r*hflowlag[[t]]
      
      Rmean[sim,,t] <- r[sim,,t] + phi[sim,,t] * N[sim,,t-1] / extent
      R[sim,,t] <- rnorm(nmcmc, Rmean[sim,,t], sigmaR)  
      
      muN[sim,,t] <- N[sim,,t-1] * exp(R[sim,,t]) * (Nfem > 0) * (Nfem < N[sim, ,t-1]) 
      muN[sim,,t][muN[sim,,t] > Nmax] <- Nmax
      
      N[sim,,t] <- rpois(nmcmc, muN[sim,,t])
    }
    
    detach(covs.sims[[sim]])
  }
  
  # Output table
  rn <- c('Nmean','Nup','Nlow','rmean','rup','rlow','phimean','phiup','philow','extinct')
  out.df <- data.frame(row.names=rn, matrix(NA, nrow=length(rn), ncol=nyears))
  names(out.df) <-  as.character(year1.sim:forecastyear) #names(templag)
  
  out.df[c('Nlow','Nmean','Nup'),1] <- quantile(as.vector(N[,,1]), probs=c(0.025, 0.5, 0.975))
  out.df[c('Nmean'),1] <- mean(as.vector(N[,,1]))
  out.df['extinct',1] <- mean(as.vector(N[,,1]) == 0)
  
  for(t in 2:nyears){
    out.df[c('Nlow','Nmean','Nup'),t] <- quantile(as.vector(N[,,t]), probs=c(0.025, 0.5, 0.975))
    out.df['Nmean',t] <- mean(as.vector(N[,,t]))
    out.df[c('rlow','rmean','rup'),t] <- quantile(as.vector(r[,,t]), probs=c(0.025, 0.5, 0.975))
    out.df['rmean',t] <- mean(as.vector(r[,,t]))
    out.df[c('philow','phimean','phiup'),t] <- quantile(as.vector(phi[,,t]), probs=c(0.025, 0.5, 0.975))
    out.df['phimean',t] <- mean(as.vector(phi[,,t]))
    out.df['extinct',t] <- mean(as.vector(N[,,t]) == 0)
  }
  
  write.csv(out.df, file=paste(outdir, popname, '.csv', sep=''))
  
  # Output plot
  year1.wdat <- (firstyear-1) + t_ti[popnum,1]
  ncol.plot.df <- year1.sim - year1.wdat
  plot.df <- cbind(data.frame(matrix(NA, nrow=nrow(out.df), ncol=ncol.plot.df)), out.df )
  names(plot.df) <- year1.wdat:forecastyear

  for (col in 1:(ncol.plot.df+1)){
    t <- as.numeric(names(plot.df)[col]) - firstyear + 1

    plot.df[c('Nlow','Nmean','Nup'),col] <- quantile(model[,paste('N[',popnum,',',t,']',sep='')], probs=c(0.025, 0.5, 0.975))
    plot.df['Nmean',col] <- mean(model[,paste('N[',popnum,',',t,']',sep='')])
    plot.df[c('rlow','rmean','rup'),col] <- quantile(model[,paste('r[',popnum,',',t,']',sep='')], probs=c(0.025, 0.5, 0.975))
    plot.df['rmean',col] <- mean(model[,paste('r[',popnum,',',t,']',sep='')])
    plot.df[c('philow','phimean','phiup'),col] <- quantile(model[,paste('phi[',popnum,',',t,']',sep='')], probs=c(0.025, 0.5, 0.975))
    plot.df['phimean',col] <- mean(model[,paste('phi[',popnum,',',t,']',sep='')])
  }

  extinct.year <- (firstyear-1) + t_ti[popnum,nt[popnum]] + 25
  extinct.prob <- round(plot.df['extinct', as.character(extinct.year)], 4)

  extinct.dat[popname, c('ExtinctRisk','ERYear','ER2015','ERForecastYear')] <- c(extinct.prob, extinct.year, plot.df['extinct','2015'], plot.df['extinct',as.character(forecastyear)])
  write.csv(extinct.dat, file='output/tables/forecast.extinct.dat.csv', row.names=F)

  datyears <- (firstyear-1) + t_ti[popnum,]
  datyears <- datyears[!is.na(datyears)]

  plotsim(dat=plot.df, popname=popname, extinct.prob=plot.df['extinct',as.character(forecastyear)], datyears=datyears, outdir=outdir, Nmax=Nmax)

  
}


