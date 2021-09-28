#rm(list=ls())
#cat("\014") 
gc()

library(rjags)
library(runjags)
library(coda)
library(lattice)

#setwd('D:/RESEARCH/2015 NASA Trout Project/wd/1.002')

dir.create('output')

load('output/objects/popnames.R')

# Load model and data
load('output/objects/jags.dat.R')
attach(jags.dat, warn.conflicts=F)

#load('output/objects/jm.R')
zm <- as.mcmc.list(jm)
save(zm, file='output/objects/zm.R')

rm(jm)
gc()

# Variables to assess
vars <- c('b0r','b1r','b2r','sigmabr',
          'b0phi','b1phi','b2phi','sigmabphi',
          'b0palpha','b1palpha','pbeta','sigmabpalpha',
          'tauprop','taup',#'sigmaprop','sigmap',
          'alpha','beta', 
          paste('sigmaR[',1:npops,']', sep='')
          )

# Trace Plots
dir.create('output/plots/trace.plots')
jpeg(paste('output/plots/trace.plots/trace.%02d.jpg', sep=''), quality=100, height=720, width=720)
for (i in seq(1, length(vars), 4)){  #seq(1, 20, 4)
  print(vars[i:(i+3)])
  plot(zm[,vars[seq(i,min(i+3, length(vars)))]])
}
dev.off()
gc()

# Convert model to data frame
d <- data.frame()
for (i in 1:length(zm)) d <- rbind(d, zm[[i]])

save(d, file='output/objects/d.R')

rm(zm)
gc()

# Parameter quantiles
q <- rbind(apply(d, 2, mean, na.rm=T), apply(d, 2, quantile, probs=c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975), na.rm=T))
parms <- t(q)
colnames(parms)[1] <- 'mean'

wbracket <- grepl(pattern='[', rownames(parms), fixed=T)
wcomma <- grepl(pattern=',', rownames(parms), fixed=T)

parms1 <- parms[!wbracket & !wcomma, ]
parms2 <- parms[wbracket & !wcomma, ]
parms3 <- parms[wbracket & wcomma, ]
parms3 <- parms3[order(row.names(parms3)),]

parms <- rbind(parms1, parms2, parms3)
write.csv(parms, file='output/tables/parms.csv')

rm(parms1, parms2, parms3, parms, wbracket, wcomma)
gc()

# # Plot detection across passes
# pdf('output/plots/p.pdf')
# 
# for(i in 1:npops){
#   for(t in 1:nt[i]){
#     for(j in 1:nj[i,t]){
#       plot(x=1:npasses[i,t,j], y=q['50%',paste('p[',i,',',t,',',j,',',1:npasses[i,t,j],']', sep='')], type='p', ylim=c(0,1), 
#            xlab='Pass (m)', ylab='Detection Probability (p)',xaxt='n',
#            main=paste('Pop: ',i, ' (',popnames[i],')','\nYear: ',t_ti[i,t], ', Site: ',j, sep=''))
#       axis(1,at=1:npasses[i,t,j])
#       
#       for (m in 1:npasses[i,t,j]) {
#         arrows(y0=q['2.5%',paste('p[',i,',',t,',',j,',',m,']', sep='')], y1=q['97.5%',paste('p[',i,',',t,',',j,',',m,']', sep='')], x0=m, x1=m, length=0)
#       }
#     }
#   }
# }
# dev.off()


# # Plot density of rprop
# library(lattice)
# jpeg('output/plots/rprop.jpg')
# 
# tau <- q['50%', 'tauprop']
# prop <- 0.05
# 
# alpha <- prop * tau
# beta <- (1-prop) * tau
# 
# x <- seq(0, 0.4, length=100)
# y <- dbeta(x, alpha, beta)
# 
# plot(x=x, y=y, type='l',lty=1, ylim=c(0,max(y)), xlab='Proportion of Population Sampled', ylab='Probability Density')
# points(x=prop, y=0, cex=2, pch=16)
# 
# prop <- 0.1
# alpha <- prop * tau
# beta <- (1-prop) * tau
# y <- dbeta(x, alpha, beta)
# lines(x=x, y=y, lty=2)
# points(x=prop, y=0, cex=2, pch=1)
# 
# prop <- 0.2
# alpha <- prop * tau
# beta <- (1-prop) * tau
# y <- dbeta(x, alpha, beta)
# lines(x=x, y=y, lty=3)
# points(x=prop, y=0, cex=2, pch=8)
# 
# 
# dev.off()

# Cleanup
rm(list=ls()[-which(ls() %in% c('d','zm'))])

