model{
  
  ##=========##
  ## INDICES ##
  ##=========##
  # i = population ID
  # t (for process model) = years since 1984
  # t (for the rest of model) = year index referencing only years with data for each population
  # j = site
  # m = pass

  ##===============##
  ## PROCESS MODEL ##
  ##===============##
  for (i in 1:npops) {
    
    ## N ~ Initialize N first year
    N[i,t_ti[i,1]-1] ~ dpois( N1[i] )
    N1[i] ~ dgamma(0.001, 0.001)
    
    for (t in t_ti[i,1]:t_ti[i,nt[i]]){
      
      ## N ~ population size + demographic stochasticity
      N[i,t] ~ dpois( min(1e6, max(0, N[i,t-1] + reintro[i,t-1]) * exp(R[i,t])) )
      
      ## R ~ realized population growth rate with density-dependence
      R[i,t] ~ dnorm(r[i,t] + phi[i,t] * max(0, N[i,t-1] + reintro[i,t-1]) / extent[i], pow(sigmaR[i],-2))
      
      ## r = intrinsic population growth rate
      r[i,t] <- b0r + b1r*templag[i,t] + b2r*hflowlag[i,t]
      
      ## phi = strength of density dependence (r/k)
      phi[i,t] <- b0phi + b1phi*bkt[i,t] + b2phi[elev[i]]*ndvi[i,t]
    }
    
    ## sigmaR[i] ~ environmental stochasticity in population i (standard deviation in realized growth rate)
    sigmaR[i] ~ dt(musig, pow(sigsig,-2), 1) T(0.01,)
  }
  musig ~ dunif(0, 10)
  sigsig ~ dunif(0, 10)
  
  sigmaR0 ~ dt(musig, pow(sigsig,-2), 1) T(0.01,)
  
  ## Regression coefficients for r
  b0r ~ dnorm(0, pow(5, -2))
  b1r ~ dnorm(0, pow(5, -2))
  b2r ~ dnorm(0, pow(5, -2))
  
  ## Regression coefficients for phi
  b0phi ~ dnorm(0, pow(0.1, -2))
  b1phi ~ dnorm(0, pow(0.1, -2))
  b2phi[1] ~ dnorm(0, pow(0.1, -2))
  b2phi[2] ~ dnorm(0, pow(0.1, -2))
  
  ##================##
  ## SAMPLING MODEL ##
  ##================##
  for(i in 1:npops){
    for(t in 1:nt[i]){
      
      ## Yi  ~ total fish observed among sites at a population
      Yi[i,t] ~ dbin( 1-prod(1-rho[i,t,1:nj[i,t]]),  N[i,t_ti[i,t]] )
      
      for(j in 1:nj[i,t]){
        
        ## n ~ true site abundance
        n[i,t,j] ~ dbin( psi[i,t,j], N[i,t_ti[i,t]] )
        
        ## psi ~ probability that an individual from the population occupied a given site
        psi[i,t,j] ~ dbeta( (sitelength[i,t,j]/extent[i]) * tau , (1 - (sitelength[i,t,j]/extent[i])) * tau ) T(1e-10, 1-1e-10)
        
        ## rho = probability of an individual occupying a given site and being detected in one of the passes
        rho[i,t,j] <- psi[i,t,j] * ( 1 - prod( 1 - p[i,t,j,1:npasses[i,t,j]] ) )
      }
    }
  }
  ## tau ~ sampling precision
  tau ~ dunif(0, 1e4)
  
  ##===================##
  ## OBSERVATION MODEL ##
  ##===================##
  for(i in 1:npops){
    for(t in 1:nt[i]){
      for(j in 1:nj[i,t]){

        ## Yj ~ total fish observed among all passes at a site
        Yj[i,t,j] ~ dbin( 1-Q[i,t,j,npasses[i,t,j]+1] , n[i,t,j] )
        
        ## y ~ fish observed each pass
        y[i,t,j,1:npasses[i,t,j]] ~ dmulti( pi[i,t,j,1:npasses[i,t,j]], Yj[i,t,j] )

        for(m in 1:npasses[i,t,j]){

          ## PI[m] = probability of detection in pass m AND non-detection in previous passes GIVEN detection occurred one of the passes
          pi[i,t,j,m] <- p[i,t,j,m] * Q[i,t,j,m] / (1 - Q[i,t,j,npasses[i,t,j]+1])

          ## p = pass-specific detection probability
          p[i,t,j,m] <- max(1e-3, ilogit(b0p + b1p*drain[i,t,j]) * exp(delta * (m - 1)) )
        }

        ## Q = probability of non-detection in all passes prior to m
        Q[i,t,j,1] <- 1
        for(m in 2:(npasses[i,t,j]+1)){
          Q[i,t,j,m] <- prod( 1 - p[i,t,j,1:(m-1)] )
        }
  }}}
  
  # Regression coefficients for detection
  b0p ~ dnorm(0, pow(5,-2)) #  dnorm(-0.4, pow(0.7,-2))
  b1p ~ dnorm(0, pow(5,-2))

  ## delta ~ decrease in detection each pass (parameter of negative exponential function)
  delta ~ dnorm(0, pow(5,-2)) #I(,0) #  dnorm(-0.4, pow(0.7,-2)) I(,0)
}

