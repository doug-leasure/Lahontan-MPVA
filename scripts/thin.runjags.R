thin.runjags <- function(jm, newthin=5){
  
  thin <- newthin*jm$thin
  
  for(c in 1:length(jm$mcmc)){
    jm$mcmc[[c]] <- window(jm$mcmc[[c]], thin=thin)
  }
  
  jm$thin <- jm$thin * newthin
  
  return(jm)  
}
