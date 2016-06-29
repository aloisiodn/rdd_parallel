#install.packages('rdd')
#install.packages('rdrobust')
#install.packages('doMC')
#require(rdrobust)
require(rdd)
require(doMC)
#x11()
library(doMC)
nr_cores <- detectCores() -1
nr_cores <-4
registerDoMC(nr_cores)

simulador=function(n, cutPoint=0, impact=-10000, fix.stddev=20000, var.stdev=40000, xmin=-100, xmax=100,coef3=0,coef2=-4,coef1=-300,coef0=200000)
{
  set.seed(171)
  sim.data = matrix(c(NA), nrow= n, ncol=5 )
  sim.data = data.frame(sim.data)
  
  sim.data[,1] = sample(c(xmin:xmax), n, replace = TRUE)
  sim.data[,3] = ifelse(sim.data[,1] > cutPoint,1,0) #Treatment
  sim.data[,4] = round(runif(n, min =1, max=5))
  sim.data[,5] = round(runif(n, min =1, max=5))
  sim.data[,2] = (sim.data[,1]^3)*coef3 +(sim.data[,1]^2)*coef2 + sim.data[,1]*coef1 + coef0 + 
    rnorm(n,0,fix.stddev)  + 
    ((sim.data[,1] - cutPoint)/(xmax-xmin)) * rnorm(n,0,var.stdev) +
    impact*sim.data[,3] + 
    impact*sim.data[,4]/2 +
    impact*sim.data[,5]/3
  
  names(sim.data) <- c('x', 'y', 'T','Z1','Z2')
  
  sim.data
}

sim = simulador(8000000)


amostra <- sim[sample(nrow(sim),2000),]

y = sim$y
x = sim$x
#z1 = sim$Z1
#z2 = sim$Z2

system.time(RDestimate(y~x))
system.time(RDestimate.par(y~x))

rd <- RDestimate(y~x)
rd$bw
rd$est
rd.par <- RDestimate.par(y~x)
rd.par$bw
rd.par$est

#

#par(mfrow=c(1,2))
#plot(amostra$x,amostra$y)
#rdplot(sim$y,sim$x)

