require(rdd)
library(doMC)
require(sm)
require(rdrobust)
load("~/rdd_parallel/par_pessoa.rda")
#save(pess, file="par_pessoa_1M.rda")


source('~/rdd_parallel/rdestimate_par.R')
source('~/rdd_parallel/rdestimate.R')
source('~/rdd_parallel/IKbandwidth_par.R')
source('~/rdd_parallel/hatvalues_par2.R')
source('~/rdd_parallel/vcovHC_par.R')

nr_cores <-9
registerDoMC(nr_cores)

n = 20

rd_res = data.frame(size = 1:n, time= 1:n, est=1:n, bw= 1:n,ci1=1:n,  ci2=1:n)
rdp_res = data.frame(size = 1:n, time= 1:n, est=1:n, bw= 1:n,ci1=1:n,  ci2=1:n)

for (i in 1:20) {
    if (i==20) {
      amostra <- pess
    } else {
      amostra <- pess[sample(nrow(pess),i*nrow(pess)/20),]
    }
    print(paste ("Sample Size:",nrow(amostra) ))

    y = amostra$VLR_RAIS
    x = amostra$DIAS18A
    z = ifelse(amostra$VLR_VARIAVEL_017_2014>0,0,1)
    
    rd.t <- system.time(rd<-RDestimate(y~x+z))
    rd_res$size[i] = nrow(amostra)
    rd_res$time[i] = rd.t[3]
    rd_res$est[i]  = rd$est[1]
    rd_res$bw[i]   = rd$bw[1]
    rd_res$ci1[i]  = rd$ci[1,1]
    rd_res$ci2[i]  = rd$ci[1,2]
    print(paste ("Seq Time:",rd.t[3] ))
    
    rdp.t <- system.time(rdp<-RDestimate.par(y~x+z))
    rdp_res$size[i] = nrow(amostra)
    rdp_res$time[i] = rdp.t[3]
    rdp_res$est[i]  = rdp$est[1]
    rdp_res$bw[i]   = rdp$bw[1]
    rdp_res$ci1[i]  = rdp$ci[1,1]
    rdp_res$ci2[i]  = rdp$ci[1,2]
    print(paste ("Par Time:",rdp.t[3] ))
}

#save(rd_res, file="bf_rd_res.rda")
#save(rdp_res, file="bf_rdp_res.rda")

load(file="bf_rd_res.rda")
load(file="bf_rdp_res.rda")

par(mfrow=c(1,3))

#time
plot(x=rd_res$size,y=rd_res$time, type="o", 
     col = "black",cex = .5, lty = 2, lwd=2, 
     main = "(a)" , ylab="Execution Time (s)",xlab=" Sample Size" )
lines(x=rdp_res$size,y=rdp_res$time,type="o",col = "black",cex = .5, lty = 1, lwd=2)
legend(x="topleft",legend=c("seq.", "parall."), lty= c(2,1), lwd=c(2,2), col = c("black","black"))


#Bandwidths
plot(x=rd_res$size,y=rd_res$bw, type="o", 
     col = "black",cex = .5, lty = 1, lwd=2,
     ylim=c(min(c(rd_res$bw, rdp_res$bw)), max(c(rd_res$bw, rdp_res$bw))),
     main = "(b)" , ylab="Bandwith",xlab=" Sample Size" )

#Estimates
plot(x=rd_res$size,y=rd_res$est, type="o", 
         col = "black",cex = .5, lty = 1, lwd=2,
     ylim=c(-75, 125),
     main = "(c)" , ylab="Estimates & Confidence Intervals",xlab="Sample Size" )
lines(x=c(0,max(rdp_res$size+1000)),y=c(0,0),type="l",col = "black",cex = .5, lty = 3, lwd=1)

lines(x=rd_res$size,y=rd_res$ci1, type="o", 
      col = "black",cex = .5, lty = 2, lwd=1 )
lines(x=rd_res$size,y=rd_res$ci2, type="o", 
      col = "black",cex = .5, lty = 2, lwd=1 )
lines(x=c(min(rdp_res$size),max(rdp_res$size)),y=c(0,0),type="o",col = "gray",cex = .5, lty = 3, lwd=1)
legend(x="topright",legend=c("estim.", "c. int."), lty= c(1,2), lwd=c(2,1), col = c("black","black"))



#######################################################
# pré análise rdd
#######################################################

ps <- pess[pess$DIAS18A>=-1000&pess$VLR_RAIS<=2000&pess$DIAS18A<=400,]

amostra <- ps[sample(nrow(ps),15000),]


y = amostra$VLR_RAIS
x = amostra$DIAS18A
z = ifelse(amostra$VLR_VARIAVEL_017_2014>0,0,1)

par(mfrow=c(1,3))
plot(x,y, main="(a)", xlab="Age (days)", ylab="Formal Salary")

set.seed(173)
amostra <- ps[sample(nrow(ps),200000),]


y = amostra$VLR_RAIS
x = amostra$DIAS18A
z = ifelse(amostra$VLR_VARIAVEL_017_2014>0,0,1)
rdplot(x=x,y=y, title="(b)", x.label="Age (days)", y.label="Formal Salary",y.lim=c(0,400), col.dots = "gray", col.lines = "black", type.dots=1 )
sm_l <- sm.regression(x[x <0], z[x<0], display="none"  )  
sm_r <- sm.regression(x[x >=0], z[x>=0], h=30, display="none"  )  
plot(c(sm_l$eval.points,sm_r$eval.points),c(sm_l$estimate,sm_r$estimate),type="l", main="(c)", ylab="Treatment Assigment (W)", xlab="Age (days)")


?rdplot

#nr_cores <-8
#registerDoMC(nr_cores)

#system.time(ff<-foreach(genero = c("Masculino","Feminino"), .combine=rbind  ) %:%
#  foreach(local = c("Urbano","Rural"), .combine=rbind )  %dopar% {
#      rd<-RDestimate.par(y~x+z,subset=(pess$genero==genero) &(pess$local==local)  )
#      data.frame(Genero=genero, Local=local,  est=rd$est[1], bw=rd$bw[1], ci_a=rd$ci[1,1], ci_b=rd$ci[1,2], obs=rd$obs[1])
#})
#ff  

