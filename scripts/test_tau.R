par(mfrow=c(2,3))

# Population Extent (km)
extent <- 2

# x values
x <- seq(0,1,length=1000)

# site lengths
length <- c(0.03, 1, 1.8)

# Precision for plots 1 and 2
tau = 1000

#-------------------------------#
# Plot 1: 30 m site, static scale
mu <- length[1] / extent

y <- dbeta(x, mu*tau, (1-mu)*tau)
plot(y~x, type='l', main=paste(round(length[1]*1e3), 'm site, static precision'), ylab='Probability Density', xlab='Proportion of population')

#-------------------------------#
# Plot 2: 1000 m site, static scale
mu <- length[2] / extent

y <- dbeta(x, mu*tau, (1-mu)*tau)
plot(y~x, type='l', main=paste(round(length[2]*1e3), 'm site, static precision'), ylab='Probability Density', xlab='Proportion of population')

#-------------------------------#
# Plot 3: 1500 m site, static scale
mu <- length[3] / extent

y <- dbeta(x, mu*tau, (1-mu)*tau)
plot(y~x, type='l', main=paste(round(length[3]*1e3), 'm site, static precision'), ylab='Probability Density', xlab='Proportion of population')



# Precision for plots 3 and 4
tau = 25

#-------------------------------#
# Plot 1: 30 m site, static scale
mu <- length[1] / extent

y <- dbeta(x, mu*tau/length[1], (1-mu)*tau/length[1])
plot(y~x, type='l', main=paste(round(length[1]*1e3), 'm site, dynamic precision'), ylab='Probability Density', xlab='Proportion of population')

#-------------------------------#
# Plot 2: 1000 m site, static scale
mu <- length[2] / extent

y <- dbeta(x, mu*tau/length[2], (1-mu)*tau/length[2])
plot(y~x, type='l', main=paste(round(length[2]*1e3), 'm site, dynamic precision'), ylab='Probability Density', xlab='Proportion of population')

#-------------------------------#
# Plot 3: 1500 m site, static scale
mu <- length[3] / extent

y <- dbeta(x, mu*tau/length[3], (1-mu)*tau/length[3])
plot(y~x, type='l', main=paste(round(length[3]*1e3), 'm site, dynamic precision'), ylab='Probability Density', xlab='Proportion of population')
