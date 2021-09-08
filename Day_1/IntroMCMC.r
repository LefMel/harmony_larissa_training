################################################  
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
num.samples <-  10000
U           <-  runif(num.samples)
X           <- -log(1-U)/2
hist(X,freq=F)
curve(dexp(x, rate=2),add=T)

# Monte carlo integration
mean(2*exp(-2*u))
1-exp(-2)



################################################  
######## Acceptrance/Rejection #################
################################################

triangle.pdf = function(x){
  ifelse((0 < x) & (x < 1), x, 
         ifelse((1 <= x) & (x < 2), 2 - x, 0)) 
}

triangle.pdf(seq(0, 2, by = 0.1))
runif(20,0,2)
triangle_plot.dat = data.frame(x = seq(-0.5, 2.5, by = 0.01),
                               y = triangle.pdf(seq(-0.5, 2.5, by = 0.01)),
                               x_sample = runif(301, 0, 2),
                               y_sample = runif(301, 0, 1))
library(ggplot2)
jpeg(filename = "AR.jpg",width = 1300,height = 1200,res = 300)
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


time<-proc.time()
triangle.sample = accept_reject(triangle.pdf, 50000,prop.val = 1 ) # Optimal
time-proc.time()


time<-proc.time()
triangle.sample = accept_reject(triangle.pdf, 50000,prop.val = 10) # Less optimal
time-proc.time()


################################################  
######## Metropolis (random-walk) ##############
################################################

jpeg(filename = "Prior_Metropolis.jpg",width = 1000,height = 1000,res=150)
plot(density(rbeta(10000,0.1,1)),main = "",xlab = "Prevalence")
dev.off()

round(table(rbinom(n = 5000000,size = 1,prob = 20/1000000))[2]/
        table(rbinom(n = 5000000,size = 1,prob = 20/1000000))[1],7)


library('coda')
metropolis <- function(burnin = 0, sample = 10000, sigma = 0.05,
                       initial_value = 0.05, plot=TRUE){
  
  stopifnot(initial_value > 0, initial_value < 1)
  stopifnot(sigma > 0)
  burnin <- as.integer(burnin); sample <- as.integer(sample)
  stopifnot(burnin >= 0)
  stopifnot(sample > 0)
  
  
  # Redefine these to work on the log scale:
  llikelihood_fun <- function(prevalence) dbinom(7, 10, prevalence, log=TRUE)
  lprior_fun <- function(prevalence) dbeta(prevalence, 2, 2, log=TRUE)
  
  parameters <- numeric(burnin+sample)
  parameters[1] <- initial_value
  current <- initial_value
  post <- llikelihood_fun(current) + lprior_fun(current)
  
  for(i in 2:(burnin+sample)){
    proposal <- rnorm(1, current, sigma)
    
    if(proposal > 0 && proposal < 1){
      newpost <- llikelihood_fun(proposal) + lprior_fun(proposal)
      accept <- newpost > post || rbinom(1, 1, exp(newpost-post))
      
      if(accept){
        current <- proposal
        post <- newpost
      }
    }
    parameters[i] <- current
    
    if(plot && burnin > 0){
      plot(1:burnin, parameters[1:burnin], type='l', col='red',
           xlim=c(1,burnin+sample), ylim=c(0,0.01),
           main='Parameter values (red:burnin, blue:sample)',
           ylab='prevalence', xlab='Iteration')
      lines((burnin+1):(burnin+sample), parameters[-(1:burnin)], col='blue')
    }else if(plot){
      plot(1:sample, parameters, type='l', col='blue',
           xlim=c(1,burnin+sample), ylim=c(0,0.01),
           main='Parameter values (red:burnin, blue:sample)',
           ylab='prevalence', xlab='Iteration')
    }
    parameters <- window(coda::as.mcmc(parameters), start=burnin+1)
    varnames(parameters) <- 'prevalence'
    return(parameters)
  }
  
}

# Simple Diagnostics (EFF, Gelmans, Mix chains, Autocor)


################################################  
######## Gibbs sampling ########################
################################################
# Dataa$temperature
y<-rnorm(1000,35,10); bary<-mean(y); n<-length(y) 
Iterations<-500
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

plot(theta[100:Iterations,1],type="l")
plot(theta[100:Iterations,2],type="l")

# Simple Diagnostics (EFF, Gelmans, Mix chains, Autocor)


