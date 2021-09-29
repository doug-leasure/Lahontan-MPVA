# cleanup
rm(list=ls()); gc(); cat("\014"); try(dev.off(), silent=T)

# working directory
setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path),'../wd'))

# libraries
library(runjags)
library(coda)

# load results
load('out/objects/jags.dat.R')
load('out/objects/jm.R')
load('out/objects/popnames.R')

# genetic years
dat <- read.csv('in/evolutionary_ecology_data.csv', stringsAsFactors=F)
dat <- dat[dat$include=="Y",c('creek_pva','year_ab','year_genetics')]

# coda data frame
x <- jm$mcmc[,which(grepl('N[',varnames(jm$mcmc), fixed=T)), drop=T]

d <- x[[1]]
for(i in 2:length(x)){
  d <- rbind(d, x[[i]])
}
d <- as.data.frame(d)

# harmonic mean function
mean_harmonic <- function(x){
  (sum(x^-1) / length(x))^-1
}

# calculations by population
result_harmonic <- data.frame(matrix(NA, nrow=0, ncol=9))
names(result_harmonic) <- c('popname','popnum','year_start','year_end','year_genetics','N_mean','N_median','N_lower','N_upper')

result_mean <- data.frame(matrix(NA, nrow=0, ncol=9))
names(result_mean) <- c('popname','popnum','year_start','year_end','year_genetics','N_mean','N_median','N_lower','N_upper')

for(i in 1:nrow(dat)){
  
  popname <- dat$creek_pva[i]
  popnum <- which(popnames == popname)
  
  cat(paste('row', i, ':', popname, popnum, '\n'))
  
  lastyear <- dat$year_ab[i]
  lastyearnum <- lastyear - 1983
  
  if(lastyearnum > jags.dat$t_ti[popnum, jags.dat$nt[popnum]]){
    warning(paste(popnum, popname, ': ab_year',lastyear,'is later than last data year',1983+jags.dat$t_ti[popnum, jags.dat$nt[popnum]]))
    lastyearnum <- jags.dat$t_ti[popnum, jags.dat$nt[popnum]]
    lastyear <- 1983 + lastyearnum
  }
  
  ts <- jags.dat$t_ti[popnum, 1]:lastyearnum
  cols <- paste0('N[',popnum,',',ts,']')
  result_row <- nrow(result_mean)+1
  
  # mean abundance
  if(length(cols) > 1){
    x <- apply(d[,cols], 1, mean)
  } else {
    x <- d[,cols]
  }
  result_mean[result_row, 'popname'] <- popnames[popnum]
  result_mean[result_row, 'popnum'] <- popnum
  result_mean[result_row, 'year_start'] <- min(ts+1983)
  result_mean[result_row, 'year_end'] <- max(ts+1983)
  result_mean[result_row, 'year_genetics'] <- dat$year_genetics[i]
  result_mean[result_row, 'N_mean'] <- mean(x)
  result_mean[result_row, 'N_median'] <- median(x)
  result_mean[result_row, 'N_lower'] <- quantile(x, probs=0.025)
  result_mean[result_row, 'N_upper'] <- quantile(x, probs=0.975)
  
  # harmonic mean abundance
  if(length(cols) > 1){
    x <- apply(d[,cols], 1, mean_harmonic)
  } else {
    x <- d[,cols]
  }
  result_harmonic[result_row, 'popname'] <- popnames[popnum]
  result_harmonic[result_row, 'popnum'] <- popnum
  result_harmonic[result_row, 'year_start'] <- min(ts+1983)
  result_harmonic[result_row, 'year_end'] <- max(ts+1983)
  result_harmonic[result_row, 'year_genetics'] <- dat$year_genetics[i]
  result_harmonic[result_row, 'N_mean'] <- mean(x)
  result_harmonic[result_row, 'N_median'] <- median(x)
  result_harmonic[result_row, 'N_lower'] <- quantile(x, probs=0.025)
  result_harmonic[result_row, 'N_upper'] <- quantile(x, probs=0.975)
}

# calculations by population x year
result_abundance <- data.frame(matrix(NA, nrow=0, ncol=8))
names(result_abundance) <- c('popname','year','popnum','yearnum','N_mean','N_median','N_lower','N_upper')

for(i in 1:jags.dat$npops){
  
  cat(paste(i, popnames[i], '\n'))
  
  ts <- jags.dat$t_ti[i,1]:jags.dat$t_ti[i,jags.dat$nt[i]]
  cols <- paste0('N[',i,',',ts,']')
  result_row <- nrow(result_mean)+1
  
  for(t in ts){
    result_row <- nrow(result_abundance) + 1
    x <- d[,paste0('N[',i,',',t,']')]
    
    result_abundance[result_row, 'popname'] <- popnames[i]
    result_abundance[result_row, 'popnum'] <- i
    result_abundance[result_row, 'year'] <- t + 1983
    result_abundance[result_row, 'yearnum'] <- t
    result_abundance[result_row, 'N_mean'] <- mean(x)
    result_abundance[result_row, 'N_median'] <- median(x)
    result_abundance[result_row, 'N_lower'] <- quantile(x, probs=0.025)
    result_abundance[result_row, 'N_upper'] <- quantile(x, probs=0.975)
  }
}

# save
write.csv(result_abundance, file='out/tables/abundances_annual.csv', row.names=F)
write.csv(result_mean, file='out/tables/abundances_mean.csv', row.names=F)
write.csv(result_harmonic, file='out/tables/abundances_harmonic.csv', row.names=F)


## TEST: Harmonic mean versus Bayesian harmonic mean

# number of random draws -- representing MCMC iterations
n <- 1e4

# data representing posterior distribution of abundance estimates for three years (a, b, c)
dat <- data.frame(a = rlnorm(n, log(500), 0.5),
                  b = rlnorm(n, log(1000), 0.5),
                  c = rlnorm(n, log(750), 0.5))

# function: harmonic mean
mean_harmonic <- function(x) (sum(x^-1)/length(x))^-1

# non-Bayesian harmonic mean calculated from the mean of Bayesian abundance estimates
result1 <- mean_harmonic(c(mean(dat$a),
                           mean(dat$b),
                           mean(dat$c)))
result1  # THIS RESULT IS ABOUT 786

# mean of Bayesian harmonic mean
result2 <- mean(apply(dat, 1, mean_harmonic))
result2  # THIS RESULT IS ABOUT 672


# demonstrate that the regular mean should come out the same either way

# non-Bayesian mean calculated from the mean of Bayesian abundance estimates
result3 <- mean(c(mean(dat$a),
                           mean(dat$b),
                           mean(dat$c)))
result3  # THIS RESULT IS ABOUT 851

# mean of Bayesian harmonic mean
result4 <- mean(apply(dat, 1, mean))
result4  # THIS RESULT IS ABOUT 851




