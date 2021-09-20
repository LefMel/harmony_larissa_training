

# Code for Hui - Walter Model

library(tidyverse)
library(runjags)
library(rjags)
runjags.options(silent.jags=TRUE, silent.runjags=TRUE)
set.seed(2021-09-20)

cleanup <- character(0)


## Model Specification


hw_definition <- c("model{
  Cross_Classified_Data ~ dmulti(prob, N)
  
  # Test1+ Test2+
	prob[1] <- (prev * ((se[1])*(se[2]))) + ((1-prev) * ((1-sp[1])*(1-sp[2])))
  
  # Test1+ Test2-
	prob[2] <- (prev * ((se[1])*(1-se[2]))) + ((1-prev) * ((1-sp[1])*(sp[2])))

  # Test1- Test2+
	prob[3] <- (prev * ((1-se[1])*(se[2]))) + ((1-prev) * ((sp[1])*(1-sp[2])))

  # Test1- Test2-
	prob[4] <- (prev * ((1-se[1])*(1-se[2]))) + ((1-prev) * ((sp[1])*(sp[2])))

  
  prev ~ dbeta(1, 1)
  se[1] ~ dbeta(1, 1)
  sp[1] ~ dbeta(1, 1)
  se[2] ~ dbeta(1, 1)
  sp[2] ~ dbeta(1, 1)

  #data# Cross_Classified_Data, N
  #monitor# prev, prob, se, sp, deviance
  #inits# prev, se, sp
}
")
cat(hw_definition, sep='', file='basic_hw.txt')

cleanup <- c(cleanup, 'basic_hw.txt')
cat(hw_definition, sep='\n')


twoXtwo <- matrix(c(36, 4, 12, 48), ncol=2, nrow=2)
twoXtwo


library(runjags)

Cross_Classified_Data <- as.numeric(twoXtwo)
N <- sum(Cross_Classified_Data)

prev <- list(chain1=0.05, chain2=0.95)
se <- list(chain1=c(0.01,0.99), chain2=c(0.99,0.01))
sp <- list(chain1=c(0.01,0.99), chain2=c(0.99,0.01))

results <- run.jags('basic_hw.txt', n.chains=2)


#[Remember to check convergence and effective sample size!]

results


pt <- plot(results)


res <- summary(results)[,c(1:3,9,11)]
res[] <- round(res, 3)
res

print(pt[["prev.plot1"]])

print(pt[["se[1].plot1"]])

print(pt[["sp[1].plot1"]])

print(pt[["deviance.plot1"]])

print(pt[["crosscorr"]])

## Label Switching

#We can force Se+Sp >= 1:
#Jouden's Index Se + Sp -1 >0:
  

#se[1] ~ dbeta(1, 1)
#sp[1] ~ dbeta(1, 1)T(1-se[1], )
#se[1] ~ dbeta(1, 1)T(1-sp[1], )
#sp[1] ~ dbeta(1, 1)


#This allows the test to be useless, but not worse than useless.

# Alternatively we can have the weakly informative priors:
  
#se[1] ~ dbeta(2, 1)
#sp[1] ~ dbeta(2, 1)

#To give the model some information that we expect the test characteristics to be closer to 100% than 0%.

#Or we can use stronger priors for one or both tests.


## Priors

#A quick way to see the distribution of a prior:
  
curve(dbeta(x, 1, 1), from=0, to=1)
qbeta(c(0.025,0.975), shape1=1, shape2=1)



#  This was minimally informative, but how does that compare to a weakly informative prior for e.g. sensitivity?
  

curve(dbeta(x, 2, 1), from=0, to=1)
qbeta(c(0.025,0.975), shape1=2, shape2=1)



  
qbeta(c(0.025,0.975), shape1=2, shape2=1)

#Or more accurately:
  

library(TeachingDemos)
hpd(qbeta, shape1=2, shape2=1)

#Credible vs confidence intervals:
  
#  - For MCMC these are usually calculated using highest posterior density (HPD) intervals
#- Therefore there is a difference between:
#  - `qbeta(c(0.025,0.975), ...)`
#- `hpd(qbeta, ...)`
#- Technically HPD intervals are credible intervals...


  
#  What about a more informative prior?
  

curve(dbeta(x, 20, 2), from=0, to=1)
qbeta(c(0.025,0.975), shape1=20, shape2=2)
hpd(qbeta, shape1=20, shape2=2)


## Choosing a prior

#What we want is e.g. Beta(20,1)

#But typically we have median and 95% confidence intervals from a paper, e.g.:
  
#  "The median (95% CI) estimates of the sensitivity and specificity of the shiny new test were 94% (92-96%) and 99% (97-100%) respectively"


#How can we generate a Beta( , ) prior from this?
  
  ## The PriorGen package
  
#  "The median (95% CI) estimates of the sensitivity and specificity of the shiny new test were 94% (92-96%) and 99% (97-100%)"


library(PriorGen)
findbeta(themedian = 0.94, percentile=0.95, percentile.value = 0.92)
hpd(qbeta, shape1=429.95, shape2=27.76)

  
curve(dbeta(x, shape1=429.95, shape2=27.76))

## Initial values

#Part of the problem before was also that we were specifying extreme initial values:
  
se <- list(chain1=c(0.01,0.99), chain2=c(0.99,0.01))
sp <- list(chain1=c(0.01,0.99), chain2=c(0.99,0.01))

#Let's change these to:


se <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))
sp <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))


## Analysing simulated data

#This is useful to check that we can recover parameter values!

# Set a random seed so that the data are reproducible:
set.seed(2021-09-20)

sensitivity <- c(0.9, 0.6)
specificity <- c(0.95, 0.9)
N <- 1000
prevalence <- 0.5

data <- tibble(Status = rbinom(N, 1, prevalence)) %>%
  mutate(Test1 = rbinom(N, 1, sensitivity[1]*Status + (1-specificity[1])*(1-Status))) %>%
  mutate(Test2 = rbinom(N, 1, sensitivity[2]*Status + (1-specificity[2])*(1-Status)))

twoXtwo <- with(data, table(Test1, Test2))
Cross_Classified_Data <- as.numeric(twoXtwo) #str(Cross_Classified_Data) = Vector



#We know that e.g. the first test has Sensitivity of 90% and Specificity of 95% - so the model *should* be able to tell us that...

# Hands - On Training Session

## Point to consider {.fragile}

#How does changing the prior distributions for the se and sp of one test affect the inference for the other test parameters?




# Exercise 1

#Simulate some data using the code given above (under "Analysing simulated data"), and run it using the following Hui-Walter model with truncated Beta(1,1) priors for sensitivity and specificity of both tests:


hw_definition <- c("model{
  Cross_Classified_Data ~ dmulti(prob, N)
  
  # Test1- Test2-
	prob[1] <- (prev * ((1-se[1])*(1-se[2]))) + ((1-prev) * ((sp[1])*(sp[2])))
	
  # Test1+ Test2-
	prob[2] <- (prev * ((se[1])*(1-se[2]))) + ((1-prev) * ((1-sp[1])*(sp[2])))

  # Test1- Test2+
	prob[3] <- (prev * ((1-se[1])*(se[2]))) + ((1-prev) * ((sp[1])*(1-sp[2])))
	
	# Test1+ Test2+
	prob[4] <- (prev * ((se[1])*(se[2]))) + ((1-prev) * ((1-sp[1])*(1-sp[2])))


  prev ~ dbeta(1, 1)
  se[1] ~ dbeta(1, 1)T(1-sp[1], )
  sp[1] ~ dbeta(1, 1)
  se[2] ~ dbeta(1, 1)T(1-sp[2], )
  sp[2] ~ dbeta(1, 1)

  #data# Cross_Classified_Data, N
  #monitor# prev, prob, se, sp, deviance
  #inits# prev, se, sp
}
")
cat(hw_definition, sep='', file='hw_truncated.txt')
cleanup <- c(cleanup, 'hw_truncated.txt')
cat(hw_definition)


#What are the results?



## Exercise 2

#- Find beta distribution priors for:

#  * Sensitivity: mean estimate = 0.9 (95% CI: 0.85 - 0.95)
#  * Specificity: mean estimate = 0.95 (95%CI: 0.92-0.97)

#- Look at these distributions using curve and hpd

#- Modify your model from exercise 1 using these priors for test 1 (leave the priors for test 2 unchanged)
#  - Make sure to name your new model something different, so that you can easily run it using either set of priors for test 1!

#- How does this affect the inference for test 2?


## Exercise 3

#Now adjust the sample size so that you have N=100, re-simulate the data, and re-run the models with both sets of priors.

#  - What do you notice about the results compared to N=1000?
  
#Also change the prevalence from 50% to 10% or 90%

#  - How does this affect your ability to estimate the sensitivity and specificity of test 2 (using strong priors for test 1)?



## Optional exercise A

#Adapt the model so that you can specify the 'hyper-priors' of the sensitivity and specificity for both tests as data

#Now pretend that the manufacturer of the test told you that Test 1 actually has these characteristics:

#  * Sensitivity = 0.95 (95% CI: 0.92 - 0.98)
#  * Specificity = 0.999 (95%CI: 0.99 - 1.00)

#Re-estimate the values you would need to use for the priors

#Now run your adapted model using these values instead (using the original dataset with N=1000 and prevalence=0.5)
#  - What effect does the change to Test 1 priors have on the posterior for Test 2?
#  - Other than comparing to the simulation parameters (which you would not know in real life!) is there any way that you can tell the priors for test 1 are not realistic?
  

# In general the best way to assess the effect of your priors is by sensitivity analysis:  change the priors slightly (for example make them less informative) and see if your posterior changes substantially.  In this case we can see that the stronger priors for test 1 move the estimates for test 2 substantially, so we would have to be extremely sure that these priors are correct in order to believe the model results.


