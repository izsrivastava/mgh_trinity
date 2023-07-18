###############################################################################################
## Attempting MCMC on Booton x Srivastava 
###############################################################################################


###############################################################################################
# PREP WORK
###############################################################################################
#Setting my working directory 
#----------------------------------------------------------------------------------------------
setwd("C:/Users/isrivastava/Desktop/Project/My_Models/MCMC")

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
install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos"))) #Click yes on the pop-up, Iz! 
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
#BELOW THE THRESHOLD - 20-35% prevalence
###############################################################################################

#Time series of "cases" 
#----------------------------------------------------------------------------------------------
human_prev <- c(rep(204, 36), rep(208, 36), rep(213, 36),
                rep(219, 36), rep(225, 36), rep(233, 36),
                rep(242, 36),rep(245, 36), rep(253, 36),
                rep(259, 36), rep(261, 36), rep(270, 36),
                rep(279, 36), rep(284, 36), rep(298, 36),
                rep(309, 36), rep (317, 36), rep(324, 36), 
                rep(338, 36), rep(350, 36))
animal_prev <- c(rep(304, 36), rep(308, 36), rep(313, 36),
                 rep(319, 36), rep(325, 36), rep(333, 36),
                 rep(342, 36), rep(345, 36), rep(353, 36),
                 rep(361, 36), rep(370, 36),
                 rep(379, 36), rep(384, 36), rep(398, 36),
                 rep(409, 36), rep(417, 36), rep(424, 36), 
                 rep(438, 36), rep(442, 36), rep(450, 36)) #This prev is factually wrong, need to change
enviro_prev <- c(rep(201, 36), rep(202, 36), rep(203, 36),
                 rep(205, 36), rep(205, 36), rep(207, 36),
                 rep(208, 36),rep(210, 36), rep(211, 36),
                 rep(212, 36), rep(213, 36),
                 rep(215, 36), rep(216, 36), rep(217, 36),
                 rep(217, 36), rep (219, 36), rep(221, 36), 
                 rep(222, 36), rep(223, 36), 
                 rep(225, 36)) #This prev is factually wrong, need to change

#Times 
#----------------------------------------------------------------------------------------------
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
#Or...I could use point estimates, given the data that I do have... 
h0 <- .25 #IC < threshold 
a0 <- 0.375 #splitting the difference
e0 <- 0.375 #splitting the difference
y0 <- c(H=h0, A=a0, E=e0)
#data for Stan 
data_sir <- list(n_days = n_days, y0=y0, t0=t0, ts=t, human_prev=human_prev, animal_prev=animal_prev, enviro_prev=enviro_prev)
#number of MCMC steps 
niter <- 2000

#Next, we compile the model and save it in a file
#----------------------------------------------------------------------------------------------
model <- stan_model("C:/Users/isrivastava/Desktop/Project/My_Models/MCMC/MCMC_test_threshold.stan")
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


#Next, we plot the marginal posterior densities to confirm the Markov chains are in agreement with one another
#----------------------------------------------------------------------------------------------
stan_dens(fit_sir_negbin, pars=pars, separate_chains=TRUE)

traceplot(fit_sir_negbin, pars, include = TRUE, unconstrain = FALSE, 
          inc_warmup = FALSE) 

#Now we plot the results of the MCMC in a box-and-whiskerly fashion
#----------------------------------------------------------------------------------------------
#First, we need to be able to work with the data in a data frame 
data <- as.data.frame(fit_sir_negbin, pars=pars)
data$gamma #Checking it looks as expected...

#Set some margins, bestie...
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
        names=c(expression(mu)),show.names=TRUE,
        medcol="black",whiskcol="black",staplecol="black",boxcol="black",outcol="firebrick",outbg="red",
        boxwex=0.5,las=2,cex.lab=1,
        xlab="Mus (rate of loss of resistance)", ylab="Magnitude of resistance loss", 
        main = "Box-and-Whisker Plot of Mu Values")


###############################################################################################
#AT THE THRESHOLD - 50-66% PREVALENCE
###############################################################################################

#Time series of "cases" 
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
human_prev <- c(rep(500, 36), rep(514, 36), rep(530, 36), rep(526, 36), rep(539, 36),
                rep(555, 36), rep(546, 36), rep(564, 36), rep (578, 36), rep(572, 36),
                rep(586, 36), rep(588, 36), rep(596, 36), rep(591, 36), rep(598, 36),
                rep(604, 36), rep(622, 36), rep(641, 36), rep(657, 36), rep(666, 36))
animal_prev <- c(rep(400, 36), rep(414, 36), rep(430, 36), rep(426, 36), rep(439, 36),
                rep(555, 36), rep(546, 36), rep(564, 36), rep (578, 36), rep(572, 36),
                rep(486, 36), rep(488, 36), rep(496, 36), rep(491, 36), rep(498, 36),
                rep(404, 36), rep(422, 36), rep(441, 36), rep(457, 36), rep(466, 36))#This prev is factually wrong, need to change
enviro_prev <- c(rep(251, 36), rep(252, 36), rep(253, 36),
                 rep(255, 36), rep(255, 36), rep(257, 36),
                 rep(258, 36),rep(260, 36), rep(261, 36),
                 rep(262, 36), rep(263, 36),
                 rep(265, 36), rep(266, 36), rep(267, 36),
                 rep(267, 36), rep(269, 36), rep(271, 36), 
                 rep(272, 36), rep(273, 36), 
                 rep(275, 36))#This prev is factually wrong, need to change

#Times 
#----------------------------------------------------------------------------------------------
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
#Or...I could use point estimates, given the data that I do have... 
h0 <- .50 #IC = threshold 
a0 <- 0.25 #splitting the difference 
e0 <- 0.25 #splitting the difference
y0 <- c(H=h0, A=a0, E=e0)
#data for Stan 
data_sir <- list(n_days = n_days, y0=y0, t0=t0, ts=t, human_prev=human_prev, animal_prev=animal_prev, enviro_prev=enviro_prev)
#number of MCMC steps 
niter <- 2000

#Next, we compile the model and save it in a file
#----------------------------------------------------------------------------------------------
model <- stan_model("C:/Users/isrivastava/Desktop/Project/My_Models/MCMC/MCMC_thresholds.stan")
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


#Next, we plot the marginal posterior densities to confirm the Markov chains are in agreement with one another
#----------------------------------------------------------------------------------------------
stan_dens(fit_sir_negbin, pars=pars, separate_chains=TRUE)
traceplot(fit_sir_negbin, pars, include = TRUE, unconstrain = FALSE, 
          inc_warmup = FALSE) 

#Now we plot the results of the MCMC in a box-and-whiskerly fashion
#----------------------------------------------------------------------------------------------
#First, we need to be able to work with the data in a data frame 
data <- as.data.frame(fit_sir_negbin, pars=pars)
data$gamma #Checking it looks as expected...
#Set some margins, bestie...
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
        main = "Box-and-Whisker Plot of Mu Values")


###############################################################################################
#ABOVE THE THRESHOLD - 70-85% PREVALENCE
###############################################################################################

#Time series of "cases" 
#----------------------------------------------------------------------------------------------
#Here, I'm going to assume a theoretical population of 1000 humans and use that assumption to convert
## the NARST percentages into integers 
human_prev <- c(rep(700, 36), rep(714, 36), rep(730, 36), rep(726, 36), rep(739, 36),
                rep(755, 36), rep(746, 36), rep(764, 36), rep (778, 36), rep(772, 36),
                rep(786, 36), rep(788, 36), rep(796, 36), rep(791, 36), rep(798, 36),
                rep(804, 36), rep(822, 36), rep(841, 36), rep(857, 36), rep(866, 36))
animal_prev <- c(rep(600, 36), rep(614, 36), rep(630, 36), rep(626, 36), rep(639, 36),
                 rep(655, 36), rep(646, 36), rep(664, 36), rep (678, 36), rep(672, 36),
                 rep(686, 36), rep(688, 36), rep(696, 36), rep(691, 36), rep(698, 36),
                 rep(704, 36), rep(722, 36), rep(741, 36), rep(757, 36), rep(766, 36))#This prev is factually wrong, need to change
enviro_prev <- c(rep(281, 36), rep(282, 36), rep(283, 36),
                 rep(285, 36), rep(285, 36), rep(287, 36),
                 rep(288, 36),rep(280, 36), rep(291, 36),
                 rep(292, 36), rep(293, 36),
                 rep(295, 36), rep(296, 36), rep(297, 36),
                 rep(297, 36), rep(299, 36), rep(301, 36), 
                 rep(302, 36), rep(303, 36), 
                 rep(305, 36))#This prev is factually wrong, need to change

#Total count
#----------------------------------------------------------------------------------------------
#In the tutorial, this is N
#But I...don't have an N for my bacterial population
#So I simply leave this blank :) 

#Times 
#----------------------------------------------------------------------------------------------
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
#Or...I could use point estimates, given the data that I do have... 
h0 <- .75 #IC < threshold 
a0 <- 0.125 #rectal swabs in animals that Booton et al. used 
e0 <- 0.125 #environmental info that Booton et al. used 
y0 <- c(H=h0, A=a0, E=e0)
#data for Stan 
data_sir <- list(n_days = n_days, y0=y0, t0=t0, ts=t, human_prev=human_prev, animal_prev=animal_prev, enviro_prev=enviro_prev)
#number of MCMC steps 
niter <- 2000

#Next, we compile the model and save it in a file
#----------------------------------------------------------------------------------------------
model <- stan_model("C:/Users/isrivastava/Desktop/Project/My_Models/MCMC/MCMC_thresholds.stan")
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


#Next, we plot the marginal posterior densities to confirm the Markov chains are in agreement with one another
#----------------------------------------------------------------------------------------------
stan_dens(fit_sir_negbin, pars=pars, separate_chains=TRUE)
traceplot(fit_sir_negbin, pars, include = TRUE, unconstrain = FALSE, 
          inc_warmup = FALSE) 

#Now we plot the results of the MCMC in a box-and-whiskerly fashion
#----------------------------------------------------------------------------------------------
#First, we need to be able to work with the data in a data frame 
data <- as.data.frame(fit_sir_negbin, pars=pars)
data$gamma #Checking it looks as expected...
#Set some margins, bestie...
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
        main = "Box-and-Whisker Plot of Mu Values")
