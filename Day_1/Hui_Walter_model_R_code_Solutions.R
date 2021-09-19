## Solution to Exercises 1 -3 & Optional Exercise A


## Solution 1


prev <- list(chain1=0.05, chain2=0.95)
se <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))
sp <- list(chain1=c(0.5,0.99), chain2=c(0.99,0.5))

results_tr_1000 <- run.jags('hw_truncated.txt', n.chains=2, sample=10000)

# Note: this is only commented out to save space in the exercise file!
# plot(results_tr_1000)
# check convergence and effective sample size, and then interpret results:
results_tr_1000


#Note that model does converge, but the effective sample size is NOT high enough with 10000 samples - we need to run for longer to get reliable results:


results_tr_1000 <- run.jags('hw_truncated.txt', n.chains=2, sample=75000)

# Note: this is only commented out to save space in the exercise file!
# plot(results_tr_1000)
# check convergence and effective sample size, and then interpret results:
results_tr_1000


#Now we can see that the 95% confidence intervals for prev, se and sp are all quite wide, but at least they do contain the simulation values!


## Solution 2

#Parameters for Sensitivity: mean estimate = 0.9 (95% CI: 0.85 - 0.95):

# [Note: this is `themean` rather than `themedian`!]


PriorGen::findbeta(themean=0.9, percentile = 0.975, percentile.value = 0.85)
hpd(qbeta, shape1=148.43, shape2=16.49)
curve(dbeta(x, 148.43, 16.49), from=0, to=1)


#Parameters for Specificity: mean estimate = 0.95 (95%CI: 0.92-0.97):


PriorGen::findbeta(themean=0.95, percentile = 0.975, percentile.value = 0.92)
hpd(qbeta, shape1=240.03, shape2=12.63)
curve(dbeta(x, 240.03, 12.63), from=0, to=1)


#Here is the updated model using the new prior values:
  
  
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
  se[1] ~ dbeta(148.43, 16.49)T(1-sp[1], )
  sp[1] ~ dbeta(240.03, 12.63)
  se[2] ~ dbeta(1, 1)T(1-sp[2], )
  sp[2] ~ dbeta(1, 1)

  #data# Cross_Classified_Data, N
  #monitor# prev, prob, se, sp, deviance
  #inits# prev, se, sp
}
")
cat(hw_definition, sep='', file='hw_stronginf.txt')
cleanup <- c(cleanup, 'hw_stronginf.txt')
cat(hw_definition)



results_si_1000 <- run.jags('hw_stronginf.txt', n.chains=2)

# Note: this is only commented out to save space in the exercise file!
# plot(results_si_1000)
# check convergence and effective sample size, and then interpret results:
results_si_1000


#Note that the 95% confidence intervals are much narrower now, including for test 2!!!  The effective sample size is also much higher because the model is more identifiable (i.e. better behaved).


## Solution 3

#We can change the sample size like so:


# Set a random seed so that the data are reproducible:
set.seed(2021-09-20)

se <- c(0.9, 0.6)
sp <- c(0.95, 0.9)
N <- 100
prevalence <- 0.5

data <- tibble(Status = rbinom(N, 1, prevalence)) %>%
  mutate(Test1 = rbinom(N, 1, se[1]*Status + (1-sp[1])*(1-Status))) %>%
  mutate(Test2 = rbinom(N, 1, se[2]*Status + (1-sp[2])*(1-Status)))

twoXtwo <- with(data, table(Test1, Test2))
Cross_Classified_Data <- as.numeric(twoXtwo)

results_si_100 <- run.jags('hw_stronginf.txt', n.chains=2)
results_tr_100 <- run.jags('hw_truncated.txt', n.chains=2)

# Remember to check convergence!
# plot(results_si_100)
# plot(results_tr_100)

# Comparison to larger dataset:
results_si_100
results_tr_100
results_si_1000
results_tr_1000


#Note that the posteriors have wider confidence intervals for the smaller dataset, particularly with the weakly informative (truncated) prior for test 1.

#With a very low prevalence:
  
  # Set a random seed so that the data are reproducible:
  set.seed(2021-09-20)

se <- c(0.9, 0.6)
sp <- c(0.95, 0.9)
N <- 1000
prevalence <- 0.1

data <- tibble(Status = rbinom(N, 1, prevalence)) %>%
  mutate(Test1 = rbinom(N, 1, se[1]*Status + (1-sp[1])*(1-Status))) %>%
  mutate(Test2 = rbinom(N, 1, se[2]*Status + (1-sp[2])*(1-Status)))

twoXtwo <- with(data, table(Test1, Test2))
Cross_Classified_Data <- as.numeric(twoXtwo)

results_lowprev <- run.jags('hw_stronginf.txt', n.chains=2)

# Remember to check convergence!
# plot(results_lowprev)

results_lowprev


#The specificity for test 2 is well estimated but the sensitivity has large confidence intervals.  This is because there are relatively few true positive samples from which sensitivity can be estimated.  The opposite is true with high prevalence i.e. it is harder to estimate specificity.


## Optional solution A

#We need two additional parameters where we fix the value in R and pass that value into JAGS as data:


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
  se[1] ~ dbeta(se_prior[1,1], se_prior[1,2])T(1-sp[1], )
  sp[1] ~ dbeta(sp_prior[1,1], sp_prior[1,2])
  se[2] ~ dbeta(se_prior[2,1], se_prior[2,2])T(1-sp[2], )
  sp[2] ~ dbeta(sp_prior[2,1], sp_prior[2,2])

  #data# Cross_Classified_Data, N, se_prior, sp_prior
  #monitor# prev, prob, se, sp, deviance
  #inits# prev, se, sp
}
")
cat(hw_definition, sep='', file='hw_hyperprior.txt')
cleanup <- c(cleanup, 'hw_hyperprior.txt')
cat(hw_definition)



# Set a random seed so that the data are reproducible:
set.seed(2021-09-20)

se <- c(0.9, 0.6)
sp <- c(0.95, 0.9)
N <- 1000
prevalence <- 0.5

data <- tibble(Status = rbinom(N, 1, prevalence)) %>%
  mutate(Test1 = rbinom(N, 1, se[1]*Status + (1-sp[1])*(1-Status))) %>%
  mutate(Test2 = rbinom(N, 1, se[2]*Status + (1-sp[2])*(1-Status)))

twoXtwo <- with(data, table(Test1, Test2))
Cross_Classified_Data <- as.numeric(twoXtwo)

# Additional data needed for the adapted model:
(se_prior <- matrix(1, ncol=2, nrow=2))
(sp_prior <- matrix(1, ncol=2, nrow=2))

# New priors from the manufacturer:
#  * Sensitivity = 0.95 (95% CI: 0.92 - 0.98)
#  * Specificity = 0.999 (95%CI: 0.99 - 1.00)

PriorGen::findbeta(themedian = 0.95, percentile=0.95, percentile.value = 0.92)
se_prior[1,] <- c(183.59, 9.98)
TeachingDemos::hpd(qbeta, shape1=se_prior[1,1], shape2=se_prior[1,2])

PriorGen::findbeta(themedian = 0.999, percentile=0.95, percentile.value = 0.99)
sp_prior[1,] <- c(199.22, 0.53)
TeachingDemos::hpd(qbeta, shape1=sp_prior[1,1], shape2=sp_prior[1,2])

# Now our hyper-prior parameters look like:
se_prior
sp_prior

# Re-run the model with these prirs:
results_manufacturer <- run.jags('hw_hyperprior.txt', n.chains=2)

# Remember to check convergence!
# plot(results_manufacturer)

results_manufacturer
results_si_1000


# Note: you might see some warnings like `Warning in if (class(temp) == "function") {: the condition has length > 1 and only the first element will be used` - this is due to a bug in runjags that will be fixed soon!

# Compared to the earlier results (results_si_1000), the sensitivity and specificty of test 1 are estimated to be higher, but the sensitivity of test 2 is estimated to be lower. In fact, the 95% CI for specificity of test 2 no longer contains the true simulation value (0.9). Having too much confidence in the performance of one test will always make all of the other tests look worse! So remember to carefully assess the studies on which you base your priors:  if these come from laboratory validation exercises (which typically involve spiked samples with extremely high concentrations of the target vs completely clean/sterile negative controls) then this may not be a realistic estimation of the sensitivity or specificity in the field.
