rm(list=ls())
cat("\014") 
gc()
#options(scipen=5) # stop R from converting to scientific notation so often

dat <- read.csv('input/jd.props.20160513.csv')$jd.props
dat <- dat[dat > 0 & dat < 1]

# Likelihood function
betalik <- function(x=dat, a, b){
  -sum(dbeta(x, a, b, log=T))
}

# MOM estimates
mu <- mean(dat)
sigma <- sd(dat)

alpha.mom <- (mu^2 - mu^3 - mu*sigma^2) / sigma^2
beta.mom <- (mu - 2*mu^2 + mu^3 - sigma^2 + mu*sigma^2) / sigma^2

nll.mom <- betalik(a=alpha.mom, b=beta.mom)

# MLE estimates
library(bbmle)
result.mle <- mle2(minuslogl=betalik, start=list(a=alpha.mom, b=beta.mom), hessian=T)

# Compare
alpha.mom
beta.mom
-nll.mom

result.mle




