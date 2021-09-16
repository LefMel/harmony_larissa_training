"model{
    y ~ dbinom(main.ap, n)                                                                                   
    pre.y ~ dbinom(main.ap, m)                                                                                
    #Uniform (non-informative) prior for apparent prevalence (ap)                                           
    main.ap  ~ dbeta(a,b)
    plessthanSetvalue <- step(perVal-main.ap)                                                                 
  }" 

load("~/GitHub/harmony_larissa_training/Day_3/Example.mcmc.RData")
# Combine multiple chains to one
library(runjags)
fitbined <- combine.mcmc(Model1.mcmc)

library(rjags)
fitbined<-mcmc(do.call(rbind, Model1.mcmc))


# Further analysis 

traceplot(fitbined)

autocorr.plot(fitbined)

geweke.plot(fitbined)

require(ggmcmc)
S2<-ggs(fitbined)

ggs_histogram(S2) # Histograms (Main)

ggs_density(S2,family = "main") # Density plots for main parameters

ggs_density(S2,family = "sub") # Density plots for study prevalences

ggs_density(S2,family = "pre") # Density plots for predictive parameters

ggs_traceplot(S2,family = "main") # Trace plots (Main)

ggs_running(S2,family = "main") # Running means plots (Main)

ggs_compare_partial(S2,family = "main") # Partial chain comparison plots (Main)

ggs_autocorrelation(S2,family = "main") # Autocorrelation plots (Main)

#ggs_pairs(S2,lower = list(continuous = "density")) #crosscorrelation plot, at least 2 continuous inputs

# Bayesplot - rstan 
#https://cran.r-project.org/web/packages/bayesplot/vignettes/visual-mcmc-diagnostics.html
