################### Session 2021 - 09 - 21 ###################


#setwd("/Users/lefmel/Documents/GitHub/harmony_larissa_training/Day_2")

# load necessary libraries
library(tidyverse)
library(runjags)
library(rjags)
library(PriorGen)
runjags.options(silent.jags=TRUE, silent.runjags=TRUE)



# To collect temporary filenames:
cleanup <- character(0)

###   Multi test & Multi-population Hui-Walter models ###

## Why stop at two tests?

#In traditional diagnostic test evaluation, one test is assumed to be a gold standard from which all other tests are evaluated

#So it makes no difference if you assess one test at a time or do multiple tests at the same time

#Using a latent class model each new test adds new information - so we should analyse all available test results in the same model

## Simulating data

#Simulating data using an arbitrary number of independent tests is quite straightforward:
  
# Parameter values to simulate:
N <- 200
sensitivity <- c(0.8, 0.9, 0.95)
specificity <- c(0.95, 0.99, 0.95)

Populations <- 2
prevalence <- c(0.25,0.5)

data <- tibble(Population = sample(seq_len(Populations), N, replace=TRUE)) %>%
  mutate(Status = rbinom(N, 1, prevalence[Population])) %>%
  mutate(Test1 = rbinom(N, 1, sensitivity[1]*Status + (1-specificity[1])*(1-Status))) %>%
  mutate(Test2 = rbinom(N, 1, sensitivity[2]*Status + (1-specificity[2])*(1-Status))) %>%
  mutate(Test3 = rbinom(N, 1, sensitivity[3]*Status + (1-specificity[3])*(1-Status))) %>%
  select(-Status)


## Model specification

#Like for two tests, except it is now a 2x2x2 table

#The model code and data format for an arbitrary number of populations (and tests) can be determined automatically 
#using the template_huiwalter function from the runjas package:

template_huiwalter(data, covariance = FALSE)
cleanup <- c(cleanup, "huiwalter_model.txt")


#Tally[1:8,p] ~ dmulti(prob[1:8,p], TotalTests[p])

# Probability of observing Test1- Test2- Test3-
#prob[1,p] <-  prev[p] * ((1-se[1])*(1-se[2])*(1-se[3]) +  (1-prev[p]) * (sp[1]*sp[2]*sp[3])
                         
# Probability of observing Test1+ Test2- Test3-
#prob[2,p] <-  prev[p] * (se[1]*(1-se[2])*(1-se[3])) + (1-prev[p]) * ((1-sp[1])*sp[2]*sp[3])
                         
# Probability of observing Test1+ Test2+ Test3+
#prob[3,p] <-  prev[p] * (se[1]*se[2]*se[3]) + (1-prev[p]) * ((1-sp[1])*(1-sp[2])*(1-sp[3]))

#- We need to take **extreme** care with these equations, and the multinomial tabulation!!!
                           
                           

## Are the tests conditionally independent?
                           
#- Example:  we have one blood (test 1), one milk ((test 2)), and one faecal test ((test 3))
                         
#* But the blood and milk test are basically the same test
#* Therefore they are more likely to give the same result
                         
                        
#- Example:  we test people for COVID using an antigen test on a nasal swab, a PCR test on a throat swab, 
#and the same antigen test on the same throat swab
                         
#* The virus may be present in the throat, nose, neither, or both
#* But we use the same antigen test twice
#* Might it cross-react with the same non-target virus?

# - In both situations we have pairwise correlation between some of the tests
                         
                         
## Dealing with correlation
                         
#It helps to consider the data simulation as a (simplified) biological process (where my parameters are not representative of real life!):
                      
# The probability of infection with COVID in two populations:
prevalence <- c(0.01,0.05)
# The probability of shedding COVID in the nose conditional on infection:
nose_shedding <- 0.8

# The probability of shedding COVID in the throat conditional on infection:
throat_shedding <- 0.8
                         
# The probability of detecting virus with the antigen test: 
antigen_detection <- 0.75
                         
# The probability of detecting virus with the PCR test: 
pcr_detection <- 0.999
                         
# The probability of random cross-reaction with the antigen test (1-sp): 
antigen_crossreact <- 0.05
                         
# The probability of random cross-reaction with the PCR test (1-sp): 
pcr_crossreact <- 0.01
  
                         
# Note:  cross-reactions are assumed to be independent!
                          
#Simulating latent states:
                           
N <- 20000
Populations <- length(prevalence)
                         
covid_data <- tibble(Population = sample(seq_len(Populations), N, replace=TRUE)) %>%
## True infection status:
mutate(Status = rbinom(N, 1, prevalence[Population])) %>%
## Nose shedding status:
mutate(Nose = Status * rbinom(N, 1, nose_shedding)) %>%
## Throat shedding status:
mutate(Throat = Status * rbinom(N, 1, throat_shedding))

#Simulating test results:
covid_data <- covid_data %>%
## The nose swab antigen test may be false or true positive:
mutate(NoseAG = case_when(
Nose == 1 ~ rbinom(N, 1, antigen_detection),
Nose == 0 ~ rbinom(N, 1, antigen_crossreact)
)) %>%
## The throat swab antigen test may be false or true positive:
mutate(ThroatAG = case_when(
Throat == 1 ~ rbinom(N, 1, antigen_detection),
Throat == 0 ~ rbinom(N, 1, antigen_crossreact) 
)) %>%
## The PCR test may be false or true positive:
mutate(ThroatPCR = case_when(
Throat == 1 ~ rbinom(N, 1, pcr_detection),
Throat == 0 ~ rbinom(N, 1, pcr_crossreact)
))

#The overall sensitivity of the tests can be calculated as follows:
covid_sensitivity <- c(
# Nose antigen:
nose_shedding*antigen_detection + (1-nose_shedding)*antigen_crossreact,
# Throat antigen:
throat_shedding*antigen_detection + (1-throat_shedding)*antigen_crossreact,
# Throat PCR:
throat_shedding*pcr_detection + (1-throat_shedding)*pcr_crossreact
)
covid_sensitivity

#The overall specificity of the tests is more straightforward:
covid_specificity <- c(
# Nose antigen:
1 - antigen_crossreact,
# Throat antigen:
1 - antigen_crossreact,
# Throat PCR:
1 - pcr_crossreact
)
covid_specificity
                         
#However:  this assumes that cross-reactions are independent!
                           
## Model specification 
#prob[1,p] <-  prev[p] * ((1-se[1])*(1-se[2])*(1-se[3]) +covse12 +covse13 +covse23) +
#                           (1-prev[p]) * (sp[1]*sp[2]*sp[3] +covsp12 +covsp13 +covsp23)
                         
                         
#prob[2,p] <- prev[p] * (se[1]*(1-se[2])*(1-se[3]) -covse12 -covse13 +covse23) +
#                           (1-prev[p]) * ((1-sp[1])*sp[2]*sp[3] -covsp12 -covsp13 +covsp23)
                         
                         
# Covariance in sensitivity between tests 1 and 2:
#covse12 ~ dunif( (se[1]-1)*(1-se[2]), min(se[1],se[2]) - se[1]*se[2] )
                         
# Covariance in specificity between tests 1 and 2:
#covsp12 ~ dunif( (sp[1]-1)*(1-sp[2]) , min(sp[1],sp[2]) - sp[1]*sp[2] )
           
                          
#It is quite easy to get the terms slightly wrong!
                           
## Template Hui-Walter

#The model code and data format for an arbitrary number of populations (and tests) can be determined automatically using the template_huiwalter function from the runjas package:
                           
template_huiwalter(covid_data %>% select(Population, NoseAG, ThroatAG, ThroatPCR), outfile = 'covidmodel.txt')
                         
#This generates self-contained model/data/initial values etc 
cleanup <- c(cleanup, 'covidmodel.txt')
cat(readLines('covidmodel.txt')[3:111], sep='\n')

                        
cat(readLines('covidmodel.txt')[-(1:111)], sep='\n')
                     
#And can be run directly from R:
results <- run.jags('covidmodel.txt') 

results
                         
res <- summary(results)[,c(1:3,9,11)]
                         
res[] <- round(res, 3)
res
                         
## Template Hui-Walter
                         
                         
#- Modifying priors must still be done directly in the model file
#- The model needs to be re-generated if the data changes
#* But remember that your modified priors will be reset
#- There must be a single column for the population (as a factor), and all of the other columns (either factor, logical or numeric) are interpreted as being test results
             
# - Covariance terms are all deactivated by default
                         
## Activating covariance terms
                         
#Find the lines for the covariances that we want to activate (i.e. the two Throat tests):
indexes <- 101:106
cleanup <- c(cleanup, "covidmodel.txt")
ml <- readLines('covidmodel.txt')
cat(gsub('\t','',ml[indexes]), sep='\n')

#And edit so it looks like:

ml[indexes] <- c('	# Covariance in sensitivity between ThroatAG and ThroatPCR tests:', '	covse23 ~ dunif( (se[1]-1)*(1-se[2]) , min(se[1],se[2]) - se[1]*se[2] )  ## if the sensitivity of these tests may be correlated', '	 # covse23 <- 0  ## if the sensitivity of these tests can be assumed to be independent','	# Covariance in specificity between ThroatAG and ThroatPCR tests:', '	covsp23 ~ dunif( (sp[1]-1)*(1-sp[2]) , min(sp[1],sp[2]) - sp[1]*sp[2] )  ## if the specificity of these tests may be correlated', '	 # covsp23 <- 0  ## if the specificity of these tests can be assumed to be independent')
cat(ml, file='covidmodel.txt', sep='\n')
ml <- readLines('covidmodel.txt')
cat(gsub('\t','',ml[indexes]), sep='\n')

                           
#You will also need to uncomment out the relevant initial values for BOTH chains (on lines 117-122 and 128-133):
                           
                     
ml <- readLines('covidmodel.txt')
cat(gsub('\t','',ml[128:133]), sep='\n')

                         
#So that they look like:

ml[c(119,122)] <- c('"covse23" <- 0', '"covsp23" <- 0')
ml[c(130,133)] <- c('"covse23" <- 0', '"covsp23" <- 0')
cat(ml, file='covidmodel.txt', sep='\n')
ml <- readLines('covidmodel.txt')
cat(gsub('\t','',ml[128:133]), sep='\n')
                         
results <- run.jags('covidmodel.txt', sample=50000)
                         
results
                         
## Practical considerations
          
                         
#- Correlation terms add complexity to the model in terms of:
#* Opportunity to make a coding mistake
#* Reduced identifiability
                         
                        
#- The template_huiwalter function helps us with coding mistakes

#- Only careful consideration of covariance terms can help us with identifiability
                           
# Practical session 4
                           
                           
## Points to consider
                        
#How does including a third test impact the inference for the first two tests?
# What happens if we include correlation between tests?
                           
# Can we include correlation if we only have 2 tests?
                           
  ## Exercise 1 {.fragile}
                         
  #Use the template_huiwalter function to look at the simple 2-test 5-population example from earlier.  Use this data simulation code:
                           
# Set a random seed so that the data are reproducible:
set.seed(2021-09-21)
                         
sensitivity <- c(0.9, 0.6)
specificity <- c(0.95, 0.9)
N <- 1000
                         
# Change the number of populations here:
Populations <- 5
# Change the variation in prevalence here:
(prevalence <- runif(Populations, min=0.1, max=0.9))

data <- tibble(Population = sample(seq_len(Populations), N, replace=TRUE)) %>%
     mutate(Status = rbinom(N, 1, prevalence[Population])) %>%
     mutate(Test1 = rbinom(N, 1, sensitivity[1]*Status + (1-specificity[1])*(1-Status))) %>%
     mutate(Test2 = rbinom(N, 1, sensitivity[2]*Status + (1-specificity[2])*(1-Status))) %>%
     select(-Status)
                         
(twoXtwoXpop <- with(data, table(Test1, Test2, Population)))
(Tally <- matrix(twoXtwoXpop, ncol=Populations))
 (TotalTests <- apply(Tally, 2, sum))
                         
template_huiwalter(data, outfile="template_2test.txt")
                         
#Look at the model code and familiarise yourself with how the model is set out (there are some small differences, but the overall code is equivalent).  Make sure you can add a .RNG.seed, modify the priors, and add a deviance monitor.  Run the model.
                         
#Now activate the correlation terms between tests 1 and 2.  Is anything different about the results?

## Solution 1
                           
#There is no particular solution to the first part of this exercise, but please ask if you have any questions about the model code that template_huiwalter generates.  
#Remember that re-running the template_huiwalter function will over-write your existing model including any changes you made, so be careful!
                           
# We can run the model as follows:
results_nocov <- run.jags("template_2test.txt")
results_nocov

#cleanup <- c(cleanup, "template_2test.txt")
#cleanup <- c(cleanup, "template_2test_cov.txt")
                         
                         
#A shortcut for activating the covariance terms is to re-run template_huiwalter as follows:
                         
template_huiwalter(data, outfile="template_2test_cov.txt", covariance=TRUE)
                         
results_cov <- run.jags("template_2test_cov.txt")
                         
results_cov
                         
                         
#Activating the covariance terms with 2 tests has made the model less identifiable, 
#and has therefore decreased the effective sample size and increased the width of the 95% CI for all of the parameters to the point that the model is no longer very useful.  
#This is not something that we recommend you do in practice, even if the two tests are known to be correlated!
                         
## Exercise 2
                         
#Simulate some Covid data based on the code given above, and analyse the data using the default priors.  
#Interpret the results and compare them to the known values:
                           
covid_sensitivity
covid_specificity
                        
#Now exclude the NoseAG test from the dataset, re-generate the model code (without covariance terms), run the model, and interpret the results.  
#How have the posteriors for the throat tests been affected by excluding the nose swab test?
                           
## Solution 2
                           
                           
#Here are the results for all 3 tests:
                           
results_3t <- run.jags('covidmodel.txt')
                         
results_3t

#You should also notice that the model has detected a positive covariance between tests 2 and 3 (the throat swab tests), 
#although the 95% CI does include zero.  Estimating covariance terms is extremely difficult for the model to do accurately.
                         
#Excluding the nasal swab test gives us these results:
                           
template_huiwalter(covid_data %>% select(Population, ThroatAG, ThroatPCR), outfile = 'covidmodel_2t.txt')
                         
results_2t <- run.jags('covidmodel_2t.txt')
        
results_2t
                        
 cleanup <- c(cleanup, 'covidmodel_2t.txt')

                         
 #The posterior estimates for sensitivity have been affected for both tests, and neither now identify the true simulation parameter:
                           
                         
 test2 <- combine.mcmc(results_2t, vars="se", return.samples = 10000)
                         
 test3 <- combine.mcmc(results_3t, vars="se[2:3]", return.samples = 10000)
                         
                        
  bind_rows(
                          
    tibble(Model = "TwoTest", ThroatAG = test2[,1], ThroatPCR = test2[,2]),                   
    tibble(Model = "ThreeTest", ThroatAG = test3[,1], ThroatPCR = test3[,2])                  
    ) %>%
    pivot_longer(c(ThroatAG, ThroatPCR), names_to = "Test", values_to = "Estimate") %>%
    ggplot() +
    aes(x = Estimate, col = Model) +
    geom_density() +
    facet_wrap( ~ Test)
                         
                         
  #However, they do look more similar to the probabilities of detection conditional on shedding:
                           
                         
  antigen_detection
                         
  pcr_detection
                         
  #Have a think about why this might be the case.  We will spend a lot of time discussing this topic tomorrow!
                           
                           
## Optional Exercise A
#Re-fit a model to this data using all three possible covse and covsp parameters between all 3 tests
                         
#What do you notice about the results?
                           
## Optional Solution A
                           
#You can either manually change all 3 covse/covsp from before, or regenerate the model using the covariance=TRUE option:
 
                                                 
template_huiwalter(covid_data %>% select(Population, NoseAG, ThroatAG, ThroatPCR), outfile='covidmodel_allcov.txt', covariance=TRUE)
cleanup <- c(cleanup, 'covidmodel_allcov.txt')
                         
results_allcov <- run.jags('covidmodel_allcov.txt')
                         
results_allcov
                         
#The effective sample size is much lower, because the model is less identifiable.  
#But otherwise the model does a reasonable job of estimating the parameters due to the large sample size, albeit with wider 95% CI then with only covariance between throat swab tests included.  
#This might not be the case with a smaller sample size!
                         
## Summary
#- Including multiple tests is technically easy
#- But philosophically more difficult!!!

#- Complexity of adding correlation terms increases non-linearly with more tests
#- Probably best to stick to correlations with biological justification?

#- Adding/removing test results may change the posterior for
                         #- Other test Se / Sp
                         #- Prevalence
                         