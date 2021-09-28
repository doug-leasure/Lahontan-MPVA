## Inits function
inits <- function(chains=1, dat=jags.dat){
  attach(dat,warn.conflicts=F)
  if (!'extinct.thresh' %in% names(dat)) extinct.thresh=10
  inits.out <- list()
  
  for (c in 1:chains){
    inits.i <- list()
    
    # Root node initials
    inits.i$b0r <- runif(1, 0, 1)
    inits.i$b1r <- runif(1, -1, 1)
    inits.i$b2r <- runif(1, -1, 1)
    
    inits.i$b0phi <- runif(1, -1e-2, -1e-5)
    inits.i$b1phi <- runif(1, -0.02, 0.02)
    inits.i$b2phi <- runif(2, -0.02, 0.02)
    
    inits.i$sigmaR <- runif(npops, 0.5, 1.5)
    inits.i$musig <- runif(1, 0.5, 1.5)
    inits.i$sigsig <- runif(1, 0.01, 1)
    
    inits.i$tau <- runif(1, 100, 1000)

    inits.i$b0p <- runif(1, -1, 1)
    inits.i$b1p <- runif(1, -3, 3)
    inits.i$delta <- runif(1, -1, 0)

    # Population size (N), site abundance (n), and occupancy probability (psi)
    inits.i$N <- matrix(NA, nrow=npops, ncol=max(t_ti, na.rm=T))
    inits.i$n <- array(NA, dim=dim(Yj))
    inits.i$psi <- array(NA, dim=dim(Yj))

    for (i in 1:npops) {
      
      # Create initials for N, n, and psi in years with data
      for(t in 1:nt[i]) {
        # gamma = probability of detection in at least one pass
        gamma <- c()
        
        for(j in 1:nj[i,t]){
          gamma[j] <- 1 - (1-inv.logit(inits.i$b0p))^npasses[i,t,j]
          
          inits.i$n[i,t,j] <- round( Yj[i,t,j] / gamma[j] )
          inits.i$psi[i,t,j] <- sitelength[i,t,j] / extent[i]
        }
        inits.i$N[i,t_ti[i,t]] <- round(Yi[i,t] / ( 1 - prod(1 - inits.i$psi[i,t,1:nj[i,t]] * gamma[1:nj[i,t]]) ))
      }
      
      # For pops where no LCT were observed, set N to 100 in all years
      if(sum(inits.i$N[i,t_ti[i,1:nt[i]]])==0) inits.i$N[i,t_ti[i,1:nt[i]]] <- 100
      
      # If N = NA or 0, set it to mean(N[i,])
      for (t in t_ti[i,1]:t_ti[i,nt[i]]){
        if(is.na(inits.i$N[i,t]) | inits.i$N[i,t]==0){
          inits.i$N[i,t] <- round(mean(inits.i$N[i,], na.rm=T))
        }
      }
      
      # If n = 0, set it to Binom(N, psi)
      for(t in 1:nt[i]) {
        for(j in 1:nj[i,t]){
          if(inits.i$n[i,t,j]==0){
            inits.i$n[i,t,j] <- rbinom(1, inits.i$N[i,t_ti[i,t]], inits.i$psi[i,t,j])
          } 
        }
      }
      
      # Initialize N in t=0
      inits.i$N1[i] <- inits.i$N[i,t_ti[i,1]]
    }
    # Define random number generator and seed
    inits.i$.RNG.name = "lecuyer::RngStream" #"base::Wichmann-Hill"
    inits.i$.RNG.seed = c
    inits.out[[c]] <- inits.i
  }
  detach(dat)
  return(inits.out)
}