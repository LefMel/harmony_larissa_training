
################### Session 2021 - 09 - 21 ###################

###   Multi-population Hui-Walter models ###

# load necessary libraries
library(tidyverse)
library(runjags)
library(rjags)
library(PriorGen)
runjags.options(silent.jags=TRUE, silent.runjags=TRUE)



# To collect temporary filenames:
cleanup <- character(0)

## Recap

#- Fitting models using MCMC is easy with JAGS / runjags

#- But we must **never forget** to check convergence and effective sample size!  -- especially when we know the model in non-identifiable
  
#  - More complex models become easy to implement

#* For example imperfect diagnostic tests, and Hui-Walter models
#* But remember to be realistic about what is possible with your data
#* Also carefully consider the influence of your priors

###   Two-population Hui-Walter models ###   


#Simulate data for the two_test-two_population setting.
  
# Set a random seed so that the data are reproducible:
set.seed(2021-09-21)

sensitivity <- c(0.9, 0.6)
specificity <- c(0.95, 0.9)
N <- 1000

# Change the number of populations here:
Populations <- 2
# Change the variation in prevalence here:
prevalence <- runif(Populations, min=0.1, max=0.9)
prevalence <- c(0.2, 0.6)

# replace = TRUE allows element to occurs twice
data <- tibble(Population = sample(seq_len(Populations), N, replace=TRUE)) %>%
  mutate(Status = rbinom(N, 1, prevalence[Population])) %>%
  mutate(Test1 = rbinom(N, 1, sensitivity[1]*Status + (1-specificity[1])*(1-Status))) %>%
  mutate(Test2 = rbinom(N, 1, sensitivity[2]*Status + (1-specificity[2])*(1-Status)))

twoXtwoXpop <- with(data, table(Test1, Test2, Population))
Tally <- matrix(twoXtwoXpop, ncol=Populations)

#Population_Size
TotalTests <- apply(Tally, 2, sum)


two_test_two_pop <- "
model{
  for(p in 1:2){
    Tally[1:4, p] ~ dmulti(prob[1:4, p], TotalTests[p])
  
    # Test1- Test2-
	  prob[1,p] <- (prev[p] * ((1-se[1])*(1-se[2]))) + ((1-prev[p]) * ((sp[1])*(sp[2])))

    # Test1+ Test2-
  	prob[2,p] <- (prev[p] * ((se[1])*(1-se[2]))) + ((1-prev[p]) * ((1-sp[1])*(sp[2])))

    # Test1- Test2+
  	prob[3,p] <- (prev[p] * ((1-se[1])*(se[2]))) + ((1-prev[p]) * ((sp[1])*(1-sp[2])))

  	 # Test1+ Test2+
	  prob[4,p] <- (prev[p] * ((se[1])*(se[2]))) + ((1-prev[p]) * ((1-sp[1])*(1-sp[2])))
	  
    prev[p] ~ dbeta(1, 1)
  }

  se[1] ~ dbeta(se_prior[1,1], se_prior[1,2])T(1-sp[1], )
  sp[1] ~ dbeta(sp_prior[1,1], sp_prior[1,2])
  se[2] ~ dbeta(se_prior[2,1], se_prior[2,2])T(1-sp[2], )
  sp[2] ~ dbeta(sp_prior[2,1], sp_prior[2,2])

  #data# Tally, TotalTests, se_prior, sp_prior
  #monitor# prev, se, sp
  #inits# prev, se, sp
}"
cat(two_test_two_pop)
cat(two_test_two_pop, file="two_test_two_pop.txt")
cleanup <- c(cleanup, "two_test_two_pop.txt")

se_prior <- matrix(1, nrow=2, ncol=2)
sp_prior <- matrix(1, nrow=2, ncol=2)

# Set up initial values for 2 populations:
se <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))
sp <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))
prev <- list(chain1=c(0.1, 0.9), chain2=c(0.9, 0.1))

# And run the model:
results_2p <- run.jags("two_test_two_pop.txt", n.chains=2)

pt_2p <- plot(results_2p)

print(pt_2p[["prev[1].plot1"]])

print(pt_2p[["se[1].plot1"]])

print(pt_2p[["sp[1].plot1"]])

print(pt_2p[["prev[2].plot1"]])

print(pt_2p[["se[2].plot1"]])

print(pt_2p[["sp[2].plot1"]])


# Multi-population Hui-Walter models

## Hui-Walter models with multiple populations

#- Basically an extension of the two-population model

# - Works best with multiple populations each with differing prevalence

# * BUT be aware of assumptions regarding constant sensitivity/specificity across populations with very different types of infections


## Different prevalence in different populations

#Simulate data for the two_test-two_population setting.

# Set a random seed so that the data are reproducible:
set.seed(2021-09-21)

sensitivity <- c(0.9, 0.6)
specificity <- c(0.95, 0.9)
N <- 1000

# Change the number of populations here:
Populations <- 5
# Change the variation in prevalence here:
#prevalence <- runif(Populations, min=0.1, max=0.9)
prevalence <- c(0.1, 0.2, 0.4, 0.6, 0.8)

# replace = TRUE allows element to occurs twice
data <- tibble(Population = sample(seq_len(Populations), N, replace=TRUE)) %>%
  mutate(Status = rbinom(N, 1, prevalence[Population])) %>%
  mutate(Test1 = rbinom(N, 1, sensitivity[1]*Status + (1-specificity[1])*(1-Status))) %>%
  mutate(Test2 = rbinom(N, 1, sensitivity[2]*Status + (1-specificity[2])*(1-Status)))

twoXtwoXpop <- with(data, table(Test1, Test2, Population))
Tally <- matrix(twoXtwoXpop, ncol=Populations)

#Population_Size
TotalTests <- apply(Tally, 2, sum)


two_test_five_pop <- "
model{
  for(p in 1:Populations){
    Tally[1:4, p] ~ dmulti(prob[1:4, p], TotalTests[p])
  
    # Test1- Test2-
	  prob[1,p] <- (prev[p] * ((1-se[1])*(1-se[2]))) + ((1-prev[p]) * ((sp[1])*(sp[2])))

    # Test1+ Test2-
  	prob[2,p] <- (prev[p] * ((se[1])*(1-se[2]))) + ((1-prev[p]) * ((1-sp[1])*(sp[2])))

    # Test1- Test2+
  	prob[3,p] <- (prev[p] * ((1-se[1])*(se[2]))) + ((1-prev[p]) * ((sp[1])*(1-sp[2])))

  	 # Test1+ Test2+
	  prob[4,p] <- (prev[p] * ((se[1])*(se[2]))) + ((1-prev[p]) * ((1-sp[1])*(1-sp[2])))
	  
    prev[p] ~ dbeta(1, 1)
  }

  se[1] ~ dbeta(se_prior[1,1], se_prior[1,2])T(1-sp[1], )
  sp[1] ~ dbeta(sp_prior[1,1], sp_prior[1,2])
  se[2] ~ dbeta(se_prior[2,1], se_prior[2,2])T(1-sp[2], )
  sp[2] ~ dbeta(sp_prior[2,1], sp_prior[2,2])

  #data# Tally, TotalTests, Populations, se_prior, sp_prior
  #monitor# prev, se, sp
  #inits# prev, se, sp
}"
cat(two_test_five_pop)
cat(two_test_five_pop, file="two_test_five_pop.txt")
cleanup <- c(cleanup, "two_test_five_pop.txt")

se_prior <- matrix(1, nrow=2, ncol=2)
sp_prior <- matrix(1, nrow=2, ncol=2)

# Set up initial values for 2 populations:
se <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))
sp <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))

#The length of initial values for `prev` in each chain is equal to the number of populations
prev <- list(chain1=c(0.1, 0.1, 0.1, 0.9, 0.9), chain2=c(0.9, 0.9, 0.9, 0.1, 0.1))

# And run the model:
results_5p <- run.jags("two_test_five_pop.txt", n.chains=2)

results_5p

pt_5p <- plot(results_2p)

print(pt_5p[["prev[1].plot1"]])

print(pt_5p[["prev[2].plot1"]])

print(pt_5p[["prev[3].plot1"]])

print(pt_5p[["prev[4].plot1"]])

print(pt_5p[["prev[5].plot1"]])

print(pt_5p[["se[1].plot1"]])

print(pt_5p[["sp[1].plot1"]])

print(pt_5p[["prev[2].plot1"]])

print(pt_5p[["se[2].plot1"]])

print(pt_2p[["sp[2].plot1"]])



## Multiple populations: assumptions

# We typically assume that the sensitivity and specificity *must* be consistent between populations
#  -- Do you have an endemic and epidemic population?
#  -- Or vaccinated and unvaccinated?
#  * If so then the assumptions might not hold!
  
#- The populations can be artificial (e.g. age groups) but must not be decided based on the diagnostic test results
#* It helps if the prevalence differs between the populations



## Multiple populations: special cases

# A small disease-free group is extremely helpful
# Contains strong data regarding specificity
# As long as specificity can be assumed to be the same in the other populations

# Check assumptions - leave one population out

## Exercise 1

#Simulate data using multiple populations. 

# First assume we have 3 and 10 populations

    #- How does this affect the confidence intervals for the diagnostic test parameters?
  
#Now change the simulated prevalence so that it varies between 0.4-0.6 instead of 0.1-0.9

    #- How does this affect the confidence intervals for the diagnostic test parameters?
  

# Set a random seed so that the data are reproducible:
set.seed(2021-09-21)

sensitivity <- c(0.9, 0.6)
specificity <- c(0.95, 0.9)
N <- 1000

# Change the number of populations here:
Populations <- 5
# Change the variation in prevalence here:
(prevalence <- runif(Populations, min=0.1, max=0.9))
(prevalence <- runif(Populations, min=0.4, max=0.6))

data <- tibble(Population = sample(seq_len(Populations), N, replace=TRUE)) %>%
  mutate(Status = rbinom(N, 1, prevalence[Population])) %>%
  mutate(Test1 = rbinom(N, 1, sensitivity[1]*Status + (1-specificity[1])*(1-Status))) %>%
  mutate(Test2 = rbinom(N, 1, sensitivity[2]*Status + (1-specificity[2])*(1-Status)))

(twoXtwoXpop <- with(data, table(Test1, Test2, Population)))
(Tally <- matrix(twoXtwoXpop, ncol=Populations))
(TotalTests <- apply(Tally, 2, sum))

multipop <- "
model{
  for(p in 1:Populations){
    Tally[1:4, p] ~ dmulti(prob[1:4, p], TotalTests[p])
  
    # Test1- Test2-
	  prob[1,p] <- (prev[p] * ((1-se[1])*(1-se[2]))) + ((1-prev[p]) * ((sp[1])*(sp[2])))

    # Test1+ Test2-
  	prob[2,p] <- (prev[p] * ((se[1])*(1-se[2]))) + ((1-prev[p]) * ((1-sp[1])*(sp[2])))

    # Test1- Test2+
  	prob[3,p] <- (prev[p] * ((1-se[1])*(se[2]))) + ((1-prev[p]) * ((sp[1])*(1-sp[2])))

  	 # Test1+ Test2+
	  prob[4,p] <- (prev[p] * ((se[1])*(se[2]))) + ((1-prev[p]) * ((1-sp[1])*(1-sp[2])))
	  
    prev[p] ~ dbeta(1, 1)
  }

  se[1] ~ dbeta(se_prior[1,1], se_prior[1,2])T(1-sp[1], )
  sp[1] ~ dbeta(sp_prior[1,1], sp_prior[1,2])
  se[2] ~ dbeta(se_prior[2,1], se_prior[2,2])T(1-sp[2], )
  sp[2] ~ dbeta(sp_prior[2,1], sp_prior[2,2])

  #data# Tally, TotalTests, Populations, se_prior, sp_prior
  #monitor# prev, se, sp
  #inits# prev, se, sp
}"
cat(multipop)
cat(multipop, file="multipopulation.txt")
cleanup <- c(cleanup, "multipopulation.txt")

se <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))
sp <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))
prev <- list(chain1=c(0.1, 0.1, 0.1, 0.9, 0.9), chain2=c(0.9, 0.9, 0.9, 0.1, 0.1))

se_prior = matrix(1, ncol=2, nrow=2)
sp_prior = matrix(1, ncol=2, nrow=2)

# And run the model:
results_prev_range <- run.jags("multipopulation.txt", n.chains=2)

results_prev_range

plot(results_prev_range)

# Change the number of populations here:
Populations <- 3
# Change the variation in prevalence here:
(prevalence <- runif(Populations, min=0.4, max=0.6))

data <- tibble(Population = sample(seq_len(Populations), N, replace=TRUE)) %>%
  mutate(Status = rbinom(N, 1, prevalence[Population])) %>%
  mutate(Test1 = rbinom(N, 1, sensitivity[1]*Status + (1-specificity[1])*(1-Status))) %>%
  mutate(Test2 = rbinom(N, 1, sensitivity[2]*Status + (1-specificity[2])*(1-Status)))

(twoXtwoXpop <- with(data, table(Test1, Test2, Population)))
(Tally <- matrix(twoXtwoXpop, ncol=Populations))
(TotalTests <- apply(Tally, 2, sum))

# Adjust initial values for 3 populations:
se <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))
sp <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))
prev <- list(chain1=c(0.1, 0.1, 0.9), chain2=c(0.9, 0.9, 0.1))

se_prior = matrix(1, ncol=2, nrow=2)
sp_prior = matrix(1, ncol=2, nrow=2)

# And run the model:
results_3p <- run.jags("multipopulation.txt", n.chains=2)

# Remember to check convergence!
results_3p

plot(results_3p)



#Note that the effective sample size is not enough here - you either need to run the model for longer in the first place, or extend it to get more samples:

# Extend the model:
#results_3p <- extend.jags(results_3p, sample=50000)

#results_3p
# Remember to check convergence!
# plot(results_3p)



#As a general rule, the more populations you have, and the more the prevalence varies between them, the better.  However, this is conditional on having a consistent sensitivity and specificity between your populations!!!

## Exercise 2

#There are 3 populations:
  
# - The CertifiedFree population is individuals that are known to be free of disease due to their geographical location. 
#You can assume that none of these animals are truly infected.

# - The NaturalInfection population is the main population of interest, 
#where some animals will be infected and some will be uninfected.

# - The ExperimentalInfection population have been experimentally infected with an extremely large infectious dose of the causative agent of your disease, 
#whilst being kept isolated from any other possible disease. You can assume that all of these individuals are infected.

##########
#To do
# Fit a 2-test 3-population model to this data using minimally informative Beta(1,1) priors.  Interpret the results.  
#How well is sensitivity and specificity estimated relative to the other datsaets we have looked at?  Why might there be a difference here?
  
#  Now exclude individuals in the ExperimentalInfection population and re-run the model.  Do the estimated results change?  Why might that be?
set.seed(2021-06-21)

N <- 1000
Populations <- 3
se1 <- c(0.8, 0.8, 0.95)
sp1 <- c(0.95, 0.95, 0.99)
se2 <- c(0.6, 0.6, 0.9)
sp2 <- c(0.9, 0.9, 0.99)
prevalence <- c(0, 0.33, 1)

multi_pop_data <- tibble(Population = sample(seq_len(Populations), N, prob=c(0.1,0.8,0.1), replace=TRUE)) %>%
  mutate(Status = rbinom(N, 1, prevalence[Population])) %>%
  mutate(Test1 = rbinom(N, 1, se1[Population]*Status + (1-sp1[Population])*(1-Status))) %>%
  mutate(Test2 = rbinom(N, 1, se1[Population]*Status + (1-sp2[Population])*(1-Status))) %>%
  mutate(Population = factor(Population, levels=1:3, labels=c("CertifiedFree", "NaturalInfection", "ExperimentalInfection"))) %>%
  select(-Status) %>%
  arrange(Population)

save(multi_pop_data, file="multi_pop_data.Rdata")

multi_pop_data %>% count(Population, Test1, Test2)

#First we need to make a minor change to the model to include #data# prev - only 1 line of this model has changed compared to before:

multipop <- "
model{
  for(p in 1:Populations){
    Tally[1:4, p] ~ dmulti(prob[1:4, p], TotalTests[p])
  
    # Test1- Test2-
	  prob[1,p] <- (prev[p] * ((1-se[1])*(1-se[2]))) + ((1-prev[p]) * ((sp[1])*(sp[2])))

    # Test1+ Test2-
  	prob[2,p] <- (prev[p] * ((se[1])*(1-se[2]))) + ((1-prev[p]) * ((1-sp[1])*(sp[2])))

    # Test1- Test2+
  	prob[3,p] <- (prev[p] * ((1-se[1])*(se[2]))) + ((1-prev[p]) * ((sp[1])*(1-sp[2])))

  	 # Test1+ Test2+
	  prob[4,p] <- (prev[p] * ((se[1])*(se[2]))) + ((1-prev[p]) * ((1-sp[1])*(1-sp[2])))
	  
    prev[p] ~ dbeta(1, 1)
  }

  se[1] ~ dbeta(se_prior[1,1], se_prior[1,2])T(1-sp[1], )
  sp[1] ~ dbeta(sp_prior[1,1], sp_prior[1,2])
  se[2] ~ dbeta(se_prior[2,1], se_prior[2,2])T(1-sp[2], )
  sp[2] ~ dbeta(sp_prior[2,1], sp_prior[2,2])

  #data# Tally, TotalTests, Populations, se_prior, sp_prior, prev
  #monitor# prev, se, sp
  #inits# prev, se, sp
}"
cat(multipop)
cat(multipop, file="multipopulation_fixprev.txt")
cleanup <- c(cleanup, "multipopulation_fixprev.txt")

(twoXtwoXpop <- with(multi_pop_data, table(Test1, Test2, Population)))
(Tally <- matrix(twoXtwoXpop, ncol=dim(twoXtwoXpop)[3]))

se_prior = matrix(1, ncol=2, nrow=2)
sp_prior = matrix(1, ncol=2, nrow=2)

## Gather all data
data <- list(
  Tally = Tally,
  TotalTests = apply(Tally, 2, sum),
  Populations = dim(twoXtwoXpop)[3],
  prev = rep(NA, dim(twoXtwoXpop)[3]),
  se_prior = matrix(1, ncol=2, nrow=2),
  sp_prior = matrix(1, ncol=2, nrow=2)
)

#Then we need to add data corresponding to the known prevalence in the first and third population, but leave it missing in the second:
data$prev <- c(0.0, NA, 1.0)
names(data$prev) <- levels(multi_pop_data$Population)
data


#Then we can set appropriate initial values and run the model:

inits <- list(
  chain1 = list(
    prev = c(NA, 0.1, NA),
    se = c(0.5, 0.99),
    sp = c(0.5, 0.99)
  ),
  chain2 = list(
    prev = c(NA, 0.9, NA),
    se = c(0.99, 0.5),
    sp = c(0.99, 0.5)
  )
)

results_ex_2 <- run.jags("multipopulation_fixprev.txt", n.chains=2, data = data, inits = inits)

# Note: this is only commented out to save space in the exercise file!
# plot(results_3p)
# check convergence and effective sample size, and then interpret results:
results_ex_2


#The 95% CI for sensitivity and specificity of both tests are extremely narrow compared to the previous dataset 
#with the same minimally informative priors.  This is because of the 1st and 3rd population:  
#because we know the prevalence is 0 or 1, the latent state is no longer latent 
#so the model is *much* better identified. 

#This has a knock-on effect to population 2 because we assume the sensitivty and specificity are the same.

#Re-running the model for the first 2 populations is quite straightforward:


(twoXtwoXpop <- with(multi_pop_data, table(Test1, Test2, Population, exclude="ExperimentalInfection")))

Populations <- dim(twoXtwoXpop)[3]
stopifnot(Populations==2)

(Tally <- matrix(twoXtwoXpop, ncol=Populations))

data <- list(
  Tally = Tally,
  TotalTests = apply(Tally, 2, sum),
  Populations = dim(twoXtwoXpop)[3],
  prev = c(0.0, NA),
  se_prior = matrix(1, ncol=2, nrow=2),
  sp_prior = matrix(1, ncol=2, nrow=2)
)

inits <- list(
  chain1 = list(
    prev = c(NA, 0.1),
    se = c(0.5, 0.99),
    sp = c(0.5, 0.99)
  ),
  chain2 = list(
    prev = c(NA, 0.9),
    se = c(0.99, 0.5),
    sp = c(0.99, 0.5)
  )
)

results_2p_ex_2 <- run.jags("multipopulation_fixprev.txt", n.chains=2, data = data, inits = inits)

# Note: this is only commented out to save space in the exercise file!
# plot(results_2p_ex_2)
# check convergence and effective sample size, and then interpret results:
results_2p_ex_2


#The sensitivity of both tests is now much more poorly estimated because we don't have a known positive population on which to base this estimate.  
#Furthermore, the posteriors have generally shifted downwards compared to the 3 population model.  We can illustrate this using a graph:
  

## Summary

#- Multiple populations setting help to estimate Se and Sp
#- Particularly if the prevalences differ

#- Populations may be artificial
#- But cannot be based on the result of either test

#- But if Se / Sp are inconsistent then we will get misleading results
#- In practice, groups with widely varying prevalence rarely have consistent Se / Sp
#- It is possible to allow Se / Sp to differ between populations, but then there is no benefit of combining the data

