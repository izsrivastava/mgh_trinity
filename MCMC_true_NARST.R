###############################################################################################
## Attempting MCMC on Booton x Srivastava 
###############################################################################################


###############################################################################################
# PREP WORK
###############################################################################################
#Setting my working directory 
#----------------------------------------------------------------------------------------------
setwd("C:/Users/isrivastava/Desktop/Project/My_Models")

#Calling up all the libraries that seem important
#----------------------------------------------------------------------------------------------
#From Booton et al. 
library(deSolve) 
library(ggplot2) 
library(tidyr)


#From Stan tutorial 
library(gridExtra)
install.packages("rstudioapi")
library(rstudioapi)
install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(StanHeaders)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() )


###############################################################################################
#Let's get Bayesian...theoretically 
###############################################################################################
#Probability of parameters = sampling distribution/likelihood * prior 

#Sampling distribution: link unique solution of ODEs with observed data = noisy estimate of reality
## Hence model noise estimate with some kind of distribution -> our sampling distribution 

#Specify a prior for each of the parameters... 
## Can check if priors are consistent with domain knowledge using prior predictive checks

#Here's what we have from Booton...
#----------------------------------------------------------------------------------------------
## dH <- gamma*LAMBDA_H*(1-H) + LAMBDA_H*beta_HH*H*(1-H) + LAMBDA_H*beta_AH*(1-H)*A + LAMBDA_H*beta_EH*(1-H)*E - mu_H*H
## dA <-  gamma*LAMBDA_A*(1-A) + LAMBDA_A*beta_AA*A*(1-A) + LAMBDA_A*beta_HA*(1-A)*H + LAMBDA_A*beta_EA*(1-A)*E - mu_A*A
## dE <-  LAMBDA_E*beta_EE*E*(1-E) + LAMBDA_E*beta_HE*(1-E)*H + LAMBDA_E*beta_AE*(1-E)*A - mu_E*E

#What's happening in Stan: 
#----------------------------------------------------------------------------------------------
#We code f(t,y), where y=(H,A,E) in the functions block -> function called sir
#We declare fixed data in the data block, then transform the data "to match our signature of sir"
#We declare the parameters, then transform the parameters 
#We evaluate the solution numerically 
#We code the prior and the sampling distributions 
#We make a block of generated quantities (predicted human prevalence, for example)


###############################################################################################
#FITTING THE MODEL IN R
###############################################################################################

#Time series of "cases" 
#----------------------------------------------------------------------------------------------
#In the tutorial, this has to do with compartment I...
#My compartment of greatest interest is compartment H, so let's...think about that
#I do have a time series for each bin, technically! 
### Should I just supply the data I have for H (my main interest), and let A and E be point estimates in the ICs? 
### Or...should I run the Markov model 3 times, one for each compartment? 
### Include a line of animal_prev, env_prev? And run generated quantities for those, too? 
### I should start with just human_prev to confirm I can run it...
human_prev <- c(rep(20923,365), rep(28574, 365)) 
#Rounded from 20922.768 -> 20923 and 28573.62 -> 28572.62 (need integer values)


#Total count
#----------------------------------------------------------------------------------------------
#In the tutorial, this is N
#But I...don't have N! 
#Just leave this blank and hope? 

#Times 
#----------------------------------------------------------------------------------------------
#n_years <- length(human_prev) // based on the tutorial -- but that's only going to give 2, which seems like...not enough? 
n_days <- length(human_prev)
t <- seq(0, n_days, by=1)
t0=0
t <- t[-1]

#Initial conditions
#----------------------------------------------------------------------------------------------
#I believe Booton et al.'s ICs were 1/1000 for each compartment 
### h0 <- 1/1000
### a0 <- 1/1000
### e0 <- 1/1000
#Or...I could use point estimates of the prevalence from NARST & Booton?
h0 <- 28574/51484 #NARST - 2021 
a0 <- 12/54 #rectal swabs in animals that Booton et al. used 
e0 <- 3/25 #environmental info that Booton et al. used 
y0 <- c(H=h0, A=a0, E=e0)
#data for Stan 
data_sir <- list(n_days = n_days, y0=y0, t0=t0, ts=t, human_prev=human_prev)
#number of MCMC steps 
niter <- 2000

#Next, we compile the model and save it in a file
#----------------------------------------------------------------------------------------------
model <- stan_model("C:/Users/isrivastava/Desktop/Project/My_Models/MCMC/MCMC_true_NARST.stan")
#This is causing an error :( 
#I'm supposed to restart R, then run this - install.packages(c("StanHeaders","rstan"),type="source")
#THE ERROR!!! IT'S GONE!!!! 
#This was the page I used to address the issue - https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Windows

#We use Stan defaults, including 4 Markov chains
fit_sir_negbin <- sampling(model, 
                           data=data_sir, 
                           iter=niter, 
                           chains=4, 
                           seed=0)

#We specify the parameters of interest
#----------------------------------------------------------------------------------------------
pars=c("gamma", 
       "LAMBDA_H", "beta_HH", "beta_AH", "beta_EH", 
       "LAMBDA_A", "beta_HA", "beta_AA", "beta_EA", 
       "LAMBDA_E", "beta_HE", "beta_AE", "beta_EE", 
       "mu")

#And start with a summary table of the results, which will give us...
#----------------------------------------------------------------------------------------------
### posterior mean, 
### standard error, 
### quantiles, 
### some useful diagnostics 
print(fit_sir_negbin, pars=pars)
data <- as.data.frame(fit_sir_negbin, pars=pars)
data$gamma

#Next, we plot the marginal posterior densities to confirm the Markov chains are in agreement with one another
#----------------------------------------------------------------------------------------------
stan_dens(fit_sir_negbin, pars=pars, separate_chains=TRUE)

#Now we plot the results of the MCMC in a box-and-whiskerly fashion
#----------------------------------------------------------------------------------------------
par(mfrow=c(1,1))
#Plots of lambdas
boxplot(data$LAMBDA_H, data$LAMBDA_A,data$LAMBDA_E,
        col="cadetblue3", 
        names=c(expression(Lambda[H]),expression(Lambda[A]),expression(Lambda[E])),show.names=TRUE,
        medcol="black",whiskcol="black",staplecol="black",boxcol="black",outcol="firebrick",outbg="red",
        boxwex=0.5,las=2, 
        xlab="Type of antibiotic use", ylab="Prevalence of antibiotic use", 
        main = "Box-and-Whisker Plot of Lambda Values")

#Plots of betas
boxplot(data$beta_HH, data$beta_AA,data$beta_EE,data$beta_HE,
        data$beta_AH, data$beta_EH,data$beta_HA,
        data$beta_EA, data$beta_AE,
        col="cadetblue3",
        names=c(expression(beta[HH]),expression(beta[AA]),expression(beta[EE]),expression(beta[HE]),expression(beta[AH]),expression(beta[EH]),expression(beta[HA]), expression(beta[EA]),expression(beta[AE])),show.names=TRUE,
        medcol="black",whiskcol="black",staplecol="black",boxcol="black",outcol="firebrick",outbg="red",
        boxwex=0.5,las=2,cex.lab=1, 
        xlab="Transmission path", ylab="Magnitude of transmission", 
        main = "Box-and-Whisker Plot of Beta Values")

#Plot of gamma 
boxplot(data$gamma,
        col="cadetblue3", 
        names=c(expression(gamma)),show.names=TRUE,
        medcol="black",whiskcol="black",staplecol="black",boxcol="black",outcol="firebrick",outbg="red",
        boxwex=0.5,las=2,cex.lab=1,
        xlab="Gamma (rate of colonization)", ylab="Prevalence of colonization", 
        main = "Box-and-Whisker Plot of Gamma Values")

#Plot of mu
boxplot(data$mu,
        col="cadetblue3", 
        names=c(expression(gamma)),show.names=TRUE,
        medcol="black",whiskcol="black",staplecol="black",boxcol="black",outcol="firebrick",outbg="red",
        boxwex=0.5,las=2,cex.lab=1,
        xlab="Mu (rate of loss of resistance)", ylab="Magnitude of resistance loss", 
        main = "Box-and-Whisker Plot of Gamma Values")


#Now, we check the utility of our model (problem-specific, includes precise estimation of a quantity or prediction of future behaviors)
#----------------------------------------------------------------------------------------------
### In general, we should check if our fitted model produces simulations that are consistent with the observed data (the idea behind posterior predictive checks)
### We sample predictions (Y_pred) from p(Y_pred | Y) and use these samples to construct a fitted curve for human infections, together with the uncertainty interval
### See if the model gives a satisfying fit to the data and the model uncertainty captures the variation of the data 
smr_pred <- cbind(as.data.frame(summary(
  fit_sir_negbin, pars="pred_human_prev", probs=c(0.05, 0.05, 0.95))$summary), t, human_prev)
colnames(smr_pred) <- make.names(colnames(smr_pred))
smr_pred

ggplot(smr_pred, mapping = aes(x = t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "pink", alpha = 0.35) +
  geom_line(mapping = aes(x = t, y = X50.), color = "pink") + 
  geom_point(mapping = aes(y = human_prev)) +
  labs(x = "Year", y = "Human prevalence of resistant bacteria")
#This seems problematic, because it's returning mean=0 for both years...
#But it's not my main focus, so I'm going to let it go for now! 


###############################################################################################
#PRIOR PREDICTIVE CHECKS 
###############################################################################################
#We could conduct a prior predictive check -- put parms of interested in generated quantities and remove sampling distribution from the model
### Without sampling distribution, parms are not fitted to the data, they are sampled from their prior distribution
### Can add a compute_likelihood to the data to make prior predictive check easy

if(computer_likelihood == 1)
  human_prev ~ neg_binomial_2(col(to_matrix(y),2), phi)
#Then compile model without likelihood and sample from it
#This gives us samples from the a priori distribution of parameters, which we can visualize