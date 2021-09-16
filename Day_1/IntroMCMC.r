set.seed(2231) # For reproducible results
u<-runif(10000,0,1)
plot(density(u))

################################################ Done  
######## Inverse transform method ##############
################################################

# Example 1 to calculate density: fx=3x^2, 0<x<1, CDF Fx(x)=x^3, Inv-CDF(u)=u^(1/3)
#F^(-1)x(u)=Fx(x)
n=10000
u=runif(n)
x=u^(1/3)
hist(x,prob=T)
#plot(x,u)

y=seq(0,1,0.01)
lines(y,3*y^2)
#plot(u,3*u^2)

# Monte carlo integration
mean(3*u^2)
1^3

# Example 2: Exponential distribution with lambda=2
# FY(x)= P(Y<x)= int 2*exp(-2*y)dy= 1-exp(-2x)
# FY^(-1)(y)=-ln(1-y)/2
lambda      <-  2
num.samples <-  10000
U           <-  runif(num.samples)
X           <- -log(1-U)/lambda
hist(X,freq=F)
curve(dexp(x, rate=lambda),add=T)

# Monte carlo integration
mean(2*exp(-2*u))
1-exp(-2)

################################################  Done
######## Simple Monte Carlo  ###################
################################################
sims=10000
x=seq(0,1,length=sims)
plot(x,exp(-x),type="l")


u<-runif(sims,0,1)
mean(exp(-u))*1 # Expected value from simulation
(1-0)^2/sims *var(exp(-u)) # Variance/efficiency of simulation
1-exp(-1) # Expected value from integration
plot(density(exp(-u))) 

Monte_Basic<-function(sims,a,b){
  stopifnot(b>a)
u<-runif(sims,a,b)
val1<-mean(exp(-u))*(b-a)
val2<-exp(-a)-exp(-b)
varMC<-(b-a)^2/sims *var(exp(-u))
return(c("MCmean",val1,"INTmean",val2,"MCvar",varMC))
}
Monte_Basic(10000,1,5)
Monte_Basic(10000,9,10)
x=seq(9,10,length=1000)
plot(x,exp(-x),type="l")


################################################   Done
######## Acceptrance/Rejection #################
################################################

triangle.pdf = function(x){
  ifelse((0 < x) & (x < 1), x, 
         ifelse((1 <= x) & (x < 2), 2 - x, 0)) 
}
x=seq(0, 2, by = 0.01)
y= triangle.pdf(x)
plot(x,y,type="l",ylab="f(x)",xlab="x")
runif(20,0,2)
triangle_plot.dat = data.frame(x = seq(-0.5, 2.5, by = 0.01),
                               y = triangle.pdf(seq(-0.5, 2.5, by = 0.01)),
                               x_sample = runif(301, 0, 2),
                               y_sample = runif(301, 0, 1))
library(ggplot2)
jpeg(filename = "AR.jpg",width = 1300,height = 1000,res = 300)
ggplot(triangle_plot.dat, aes(x = x, y = y))+
  geom_line() + xlab("x") + ylab("Probability density function")+
  geom_point(aes(x_sample, y_sample), shape = 1,size=0.9)+
  geom_line(data=data.frame(x=0:2,y=1), linetype="dotted",size=1.2)+
  geom_line(data=data.frame(x=0,y=0:1), linetype="dotted",size=1.2)+
  geom_line(data=data.frame(x=2,y=0:1), linetype="dotted",size=1.2)
dev.off()

accept_reject = function(fx, n = 100, prop.val=1) {
  #simulates a sample of size n from the pdf fx using the accept reject algorithm
  x = numeric(n)
  count = 0
  
  while(count < n) {
    temp.x <- runif(1, 0, 2)
    u <- runif(1, 0, 1)
    if (u < fx(temp.x)) {
      count = count + 1
      x[count] <- temp.x
    }
  }
  
  return(x)
}

triangle.sample = accept_reject(triangle.pdf, 90000,prop.val = 1)


ggplot(triangle_plot.dat, aes(x = x, y = y)) +
  geom_line(color = "blue", size = 1.5) + xlab("x") + ylab("pdf") +
  geom_histogram(data = data.frame(x = triangle.sample), aes(x = x, y = ..density..), col = "gray" )


plot(density(triangle.sample,adjust = 0.5),main = "")

length(triangle.sample)/50000 # Probability of acceptance

time<-proc.time()
triangle.sample = accept_reject(triangle.pdf, 50000, prop.val = 1 ) # Optimal
time-proc.time()


time<-proc.time()
triangle.sample = accept_reject(triangle.pdf, 50000, prop.val = 10) # Less optimal
time-proc.time()


################################################  Done
######## Metropolis (random-walk) ##############
################################################ M Denwood tweaked example

jpeg(filename = "Prior_Metropolis.jpg",width = 1000,height = 1000,res=150)
plot(density(rbeta(10000,0.1,1)),main = "",xlab = "Prevalence")
dev.off()

round(table(rbinom(n = 5*10^7,size = 1,prob = 50/10^4))[2]/
        table(rbinom(n = 5*10^7,size = 1,prob = 50/10^4))[1],7)


library('coda')
metropolis <- function(burnin = 0, sample = 1000000, sigma = 0.01,
                       initial_value = 0.01, plot=TRUE){
  stopifnot(initial_value > 0, initial_value < 1)
  stopifnot(sigma > 0)
  burnin <- as.integer(burnin)
  sample <- as.integer(sample)
  stopifnot(burnin >= 0)
  stopifnot(sample > 0)
  # Redefine these to work on the log scale:
  llikelihood_fun <- function(prevalence)
    dbinom(50, 10^4, prevalence, log=TRUE)
  lprior_fun <- function(prevalence)
    dbeta(prevalence, 0.1, 1, log=TRUE)
  parameters <- numeric(burnin+sample)
  parameters[1] <- initial_value
  current <- initial_value
  post <- llikelihood_fun(current) + lprior_fun(current)
  for(i in 2:(burnin+sample)){
    proposal <- rnorm(1, current, sigma)
    if(proposal > 0 && proposal < 1){
      U=log(runif(1,0,1))
      newpost <- llikelihood_fun(proposal) + lprior_fun(proposal)
      accept <- U <= newpost-post      
      if(accept){
        current <- proposal
        post <- newpost
      }
    }
    parameters[i] <- current
  }
  
  if(plot && burnin > 0){
    plot(1:burnin, parameters[1:burnin], type='l', col='red',
         xlim=c(1,burnin+sample), #  ylim=c(0,1)
         main='Parameter values (red:burnin, blue:sample)',
         ylab='prevalence', xlab='Iteration')
    lines((burnin+1):(burnin+sample), parameters[-(1:burnin)], col='blue')
  }else if(plot){
    plot(1:sample, parameters, type='l', col='blue',
         xlim=c(1,burnin+sample),# ylim=c(0,1),
         main='Parameter values (red:burnin, blue:sample)',
         ylab='prevalence', xlab='Iteration')
  }
  parameters_coda <- window(coda::as.mcmc(parameters), start=burnin+1)
  varnames(parameters_coda) <- 'prevalence'
  return(list(parameters_coda=parameters_coda,parameters=parameters))
}

jpeg(filename = "Metropolis_Results.jpg",width = 1000,height = 1000,res=150)
set.seed(2231)
samples1 <- metropolis(burnin = 5000, sample = 10000, initial_value=0.01)
dev.off()

samples2 <- metropolis(burnin = 5000, sample = 10000, initial_value=0.0001)

jpeg(filename = "Metropolis_Results_Mix.jpg",width = 1000,height = 1000,res=150)
plot(samples1$parameters,type = 'l',col="blue",ylab = "Prevalence",xlab="Iteration",
     main = "Multiple MCMC chains")
lines(samples2$parameters,type = 'l',col="red")
dev.off()


set.seed(2231)
samples1 <- metropolis(burnin = 5000, sample = 10000, initial_value=0.01)
samples2 <- metropolis(burnin = 5000, sample = 10000, initial_value=0.0001)
mean(samples1)
effectiveSize(samples1$parameters_coda) # 10000/647.67=15.43 iterations to get an independent sample
HPDinterval(samples1$parameters_coda) # True posterior 95% credible interval
autocorr(samples$parameters_coda, lags=1)
autocorr.plot(samples$parameters_coda)


# Changing the variance of the candidate value.
par(mfrow=c(1,3))
samples_SmallS <- metropolis(sigma=0.00001, sample = 10000) # Moves only close
samples_NormalS <- metropolis(sigma=0.01, sample = 10000) # Decent movement
samples_largeS <- metropolis(sigma=1, sample = 10000) # Difficulty to identify move
par(mfrow=c(1,3))
geweke.plot(samples_SmallS$parameters_coda) # geweke.diag
geweke.plot(samples_NormalS$parameters_coda)
geweke.plot(samples_largeS$parameters_coda)

1-pnorm(abs(geweke.diag(samples_largeS$parameters_coda)[1]$z)) # p-value for convergence.

autocorr.plot(samples_SmallS$parameters_coda, lags=1) # autocorr
autocorr.plot(samples_NormalS$parameters_coda, lags=1)
autocorr.plot(samples_largeS$parameters_coda, lags=1)


gelman.diag() # For multiple chains




################################################  
######## Gibbs sampling ########################
################################################ I Ntzoufras tweaked example
# Dataa$temperature
burnin=450
y<-c(32,36,37,34,38,36,33,36,37,35,32,35); bary<-mean(y); n<-length(y) 
Iterations<-20000
muO<-0; s0<-100; a0<-0.001; b0<-0.001

theta <- matrix(nrow=Iterations, ncol=2) 
cur.mu<-0; cur.tau<-1; cur.s<-sqrt(1/cur.tau)

for (t in 1:Iterations){ 
  w<- s0^2/( cur.s^2/n+ s0^2)
  m <- w*bary + (1-w)*muO 
  s <- sqrt( w * cur.s^2/n)
  cur.mu <- rnorm(1, m, s)
  a <- a0 + 0.5*n
  b <- b0 + 0.5 * sum( (y-cur.mu)^2)
  cur.tau <- rgamma( 1, a, b)
  cur.s <- sqrt(1/cur.tau) 
  theta[t,]<-c( cur.mu, cur.s)
}

jpeg(filename = "Gibbs_Results_Uni.jpg",width = 1000,height = 700,res=150)
par(mfrow=c(1,2))
plot(theta[450:Iterations,1],type="l",main = "mu", ylab="")
plot(theta[450:Iterations,2],type="l",main = "sigma", ylab="")
dev.off()

jpeg(filename = "Gibbs_Results_both.jpg",width = 1000,height = 1000,res=150)
par(mfrow=c(2,2))
plot(theta[450:459,1],theta[450:459,2],
     type="l",main = "[10]", ylab="sigma", xlab="mu",lwd=0.1)
plot(theta[450:474,1],theta[450:474,2],
     type="l",main = "[25]", ylab="sigma", xlab="mu",lwd=0.1)
plot(theta[450:649,1],theta[450:649,2],
     type="l",main = "[200]", ylab="sigma", xlab="mu",lwd=0.1)
plot(theta[450:Iterations,1],theta[450:Iterations,2],
     type="p",main = "[19550]", ylab="sigma", xlab="mu",lwd=0.1)
dev.off()

jpeg(filename = "Gibbs_Results_both_stepbystep.jpg",width = 1000,height = 1000,res=150)
plot(rep(theta[450:459,2],each=2),c(NA,rep(theta[450:459,1],each=2)[-20]),
     type="l",main = "[10]", ylab="sigma", xlab="mu",lwd=0.1)

dev.off()


parameters_coda <- window(coda::as.mcmc(theta), start=burnin+1)
varnames(parameters_coda) <- c('mu','sigma')


effectiveSize(parameters_coda) # XX iterations to get an independent sample
HPDinterval(parameters_coda) # True posterior 95% credible interval
autocorr(parameters_coda, lags=1)
autocorr.plot(parameters_coda)
geweke.plot(parameters_coda)
geweke.diag(parameters_coda) # z-score
1-pnorm(abs(geweke.diag(parameters_coda)[1]$z)) # p-value for convergence.
