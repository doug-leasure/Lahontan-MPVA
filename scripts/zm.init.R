zm.init <- function(zm, init, burn=0, thinmore=1, addinit=T, 
                    varlist=c('b0r','b1r','b2r','b0phi','b1phi','b2phi[1]','b2phi[2]','musig','sigsig','tau','b0p','b1p','delta',paste('sigmaR[',1:jags.dat$npops,']',sep=''))){
  
  # zm=jm$mcmc; burn=0; thinmore=1; addinit=T;
  # varlist=c('b0r','b1r','b2r','b0phi','b1phi','b2phi[1]','b2phi[2]','musig','sigsig','tau','b0p','b1p','delta',paste('sigmaR[',1:jags.dat$npops,']',sep=''))
  
  thin.orig <- thin(zm)
  
  varlist <- varlist[varlist %in% varnames(zm)]
  zm <- zm[,varlist]
  
  b <- row.names(zm[[1]])[which.min(abs(as.numeric(row.names(zm[[1]])) - burn ))]
  rn <- row.names(zm[[1]])
  
  zm <- window(zm, start=as.numeric(b), 
               end=as.numeric(rn[length(rn)]), 
               thin=thin.orig*thinmore)
  
  if(addinit){
    for(i in 1:length(zm)){
      newrow <- zm[[i]][1,]
      newrow[] <- NA
      for (var in varlist[!grepl('sigmaR', varlist, fixed=T)]){
        if(!grepl('b2phi', var, fixed=T)){
          newrow[var] <- init[[i]][[var]]
        } else {
          i1 <- gsub('b2phi[','',var,fixed=T)
          i1 <- as.numeric(gsub(']','',i1,fixed=T))
          newrow[var] <- init[[i]][['b2phi']][i1]
        }
      }
      for (var in varlist[grepl('sigmaR',varlist, fixed=T)]){
        popnum <- gsub('sigmaR[','',var, fixed=T)
        popnum <- as.numeric(gsub(']','',popnum, fixed=T))
        newrow[var] <- init[[i]]$sigmaR[popnum]
      }
      zm[[i]] <- rbind(newrow, zm[[i]])  
      row.names(zm[[i]])[1] <- '1'
      zm[[i]] <- mcmc(zm[[i]], start=as.numeric(row.names(zm[[i]])[2]) , end=as.numeric(row.names(zm[[i]])[nrow(zm[[i]])]) ,thin=thin.orig*thinmore)
    }
  }
  return(zm)
}
