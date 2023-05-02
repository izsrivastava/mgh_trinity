################################################################################################
# LHS Practice
################################################################################################
#Set-up
setwd("C:/Users/isrivastava/Desktop/Project")
library(deSolve) 
library(ggplot2) 
library(tidyr)
library(data.table)
install.packages("lhs")
library(lhs)
install.packages("dplyr")
library(dplyr)
require(gridExtra)

#########################################################################################################################################################
# Now we try the parallel package...
#########################################################################################################################################################
library(parallel)
detectCores() #8

#Sampling function
f <- function(aaa){
  Samples_SIR <- randomLHS(aaa, 17) #from https://www.sciencedirect.com/science/article/pii/S2542519619301305
  
  Samples_SIR[,2] <- 0.1 + (1-0.1)*Samples_SIR[,2] #LAMBDA_A #vary to max
  Samples_SIR[,1] <- 0 + (Samples_SIR[,2]-0)*Samples_SIR[,1]  #LAMBDA_H #vary to max
  Samples_SIR[,3] <- 0 + (Samples_SIR[,1]-0)*Samples_SIR[,3] #LAMBDA_E #vary to max
  
  Samples_SIR[,4] <- 0.3 + (1-0.3)*Samples_SIR[,4]  #beta_HH
  Samples_SIR[,5] <- 0 + (1-0)*Samples_SIR[,5]  #beta_AA
  
  
  Samples_SIR[,7] <- 0 + (Samples_SIR[,4]-0)*Samples_SIR[,7] #beta_AH
  Samples_SIR[,9] <- 0 + (Samples_SIR[,7]-0)*Samples_SIR[,9] #beta_HA #beta_AA should always be bigger than beta_HA
  
  Samples_SIR[,10] <- 0 + (Samples_SIR[,5]-0)*Samples_SIR[,10] #beta_EA
  Samples_SIR[,8] <- 0 + (Samples_SIR[,4]-0)*Samples_SIR[,8]  #beta_EH
  
  Samples_SIR[,6] <- 0 + (1-0)*Samples_SIR[,6]  #beta_HE
  Samples_SIR[,11] <- 0 + (1-0)*Samples_SIR[,11]  #beta_AE
  
  Samples_SIR[,12] <- 0 + (0.5-0)*Samples_SIR[,12]  #mu_H #Carriage of ESBL or CRE at 12 months  community 25.4% patients 35.2%
  Samples_SIR[,13] <-  Samples_SIR[,12]#0 + (1-0)*Samples_SIR[,13]  #mu_A
  Samples_SIR[,14] <-  Samples_SIR[,12]  #0 + (1-0)*Samples_SIR[,14]    #mu_E
  Samples_SIR[,15] <- 0 + (1-0)*Samples_SIR[,15]    #gamma
  Samples_SIR[,16] <-  0 + (1-0)*Samples_SIR[,16] #beta_EE
  Samples_SIR[,17] <- 0 + (1-0)*Samples_SIR[,17] #epsilon
  
  
  paramsMat_SIR <- data.frame(
    LAMBDA_H=Samples_SIR[,1],LAMBDA_A=Samples_SIR[,2],LAMBDA_E=Samples_SIR[,3],
    beta_HH=Samples_SIR[,4],  beta_AA=Samples_SIR[,5],   beta_HA=Samples_SIR[,9], beta_AH=Samples_SIR[,7],
    beta_HE=Samples_SIR[,6],  beta_EH=Samples_SIR[,8],    beta_AE=Samples_SIR[,11],   beta_EA=Samples_SIR[,10],
    mu_H = Samples_SIR[,12], mu_A =Samples_SIR[,13], mu_E = Samples_SIR[,14],gamma = Samples_SIR[,15],beta_EE = Samples_SIR[,16],epsilon=Samples_SIR[,17]
  )
  return(paramsMat_SIR)
}

#AMR function
state <- c(H=(1/1000),A=(1/1000),E=(1/1000)) 
epid.start<- 2004
epid.duration <- 50
vectTime <- c(0,1:epid.duration)

AMRmodel <- function(time,state,parameters){ #using package deSolve
  with(as.list(c(state,parameters)),{
    dH <- gamma*LAMBDA_H*(1-H) + LAMBDA_H*beta_HH*H*(1-H) + LAMBDA_H*beta_AH*(1-H)*A + LAMBDA_H*beta_EH*(1-H)*E - mu_H*H
    dA <-  gamma*LAMBDA_A*(1-A) + LAMBDA_A*beta_AA*A*(1-A) + LAMBDA_A*beta_HA*(1-A)*H + LAMBDA_A*beta_EA*(1-A)*E - mu_A*A
    dE <-  LAMBDA_E*beta_EE*E*(1-E) + LAMBDA_E*beta_HE*(1-E)*H + LAMBDA_E*beta_AE*(1-E)*A - mu_E*E
    return(  list(c(dH,dA,dE)))
  })}

#Epid function
epid <- function(LAMBDA_H, LAMBDA_A, LAMBDA_E,
                 beta_HH, beta_AA, 
                 beta_HA, beta_AH, 
                 beta_HE, beta_EH, 
                 beta_AE, beta_EA,
                 mu_H, mu_A, mu_E,gamma,returnout,beta_EE,epsilon
                 
){params <- c(LAMBDA_H=LAMBDA_H,LAMBDA_A=LAMBDA_A,LAMBDA_E=LAMBDA_E,
              beta_HH=beta_HH,  beta_AA=beta_AA,  beta_HE=beta_HE,  beta_AH=beta_AH,
              beta_EH=beta_EH,  beta_HA=beta_HA,  beta_EA=beta_EA,  beta_AE=beta_AE,  
              mu_H = mu_H, mu_A = mu_A, mu_E = mu_E,gamma=gamma,beta_EE=beta_EE,epsilon=epsilon)

#run the model using desolve
out <- as.data.frame(ode(y=state,time=vectTime,func=AMRmodel,parms=params))
out$time = out$time+epid.start #rescale the time so that it runs from 2005 onwards 
#Average compartment contents between 2012-2013
model2012.H <- out$H[out$time==2012] #(out$H[out$time==2013] +out$H[out$time==2012])/2
model2012.A <- out$A[out$time==2012] #(out$A[out$time==2013] +out$A[out$time==2012])/2
model2012.E <- out$E[out$time==2012] #(out$E[out$time==2013] +out$E[out$time==2012])/2

model2005.H <- out$H[out$time==2005] 
model2005.A <- out$A[out$time==2005] 
model2005.E <- out$E[out$time==2005] 
#model2004.H <- out$H[out$time==2004] 
#model2004.A <- out$A[out$time==2004] 
#model2004.E <- out$E[out$time==2004] 

model2008.H <- out$H[out$time==2008]
model2009.H <- out$H[out$time==2009] 
model2010.H <- out$H[out$time==2010]
model2011.H <- out$H[out$time==2011]

if (returnout ==1){
  return(out)} else{
    return(c(model2012.H= model2012.H,
             model2012.A= model2012.A,
             model2012.E= model2012.E,
             model2005.H=model2005.H,
             model2005.A=model2005.A,
             model2005.E=model2005.E,
             # model2004.H=model2004.H,
             # model2004.A=model2004.A,
             # model2004.E=model2004.E,
             model2008.H=model2008.H,
             model2009.H=model2009.H,
             model2010.H=model2010.H,
             model2011.H=model2011.H))}
}

cl <- makeCluster(3)
clusterEvalQ(cl, {
  library(deSolve) 
  library(ggplot2) 
  library(tidyr)
  library(data.table)
  library(lhs)
  library(dplyr)
  
  f <- function(aaa){
    Samples_SIR <- randomLHS(aaa, 17) #from https://www.sciencedirect.com/science/article/pii/S2542519619301305
    
    Samples_SIR[,2] <- 0.1 + (1-0.1)*Samples_SIR[,2] #LAMBDA_A #vary to max
    Samples_SIR[,1] <- 0 + (Samples_SIR[,2]-0)*Samples_SIR[,1]  #LAMBDA_H #vary to max
    Samples_SIR[,3] <- 0 + (Samples_SIR[,1]-0)*Samples_SIR[,3] #LAMBDA_E #vary to max
    
    Samples_SIR[,4] <- 0.3 + (1-0.3)*Samples_SIR[,4]  #beta_HH
    Samples_SIR[,5] <- 0 + (1-0)*Samples_SIR[,5]  #beta_AA
    
    
    Samples_SIR[,7] <- 0 + (Samples_SIR[,4]-0)*Samples_SIR[,7] #beta_AH
    Samples_SIR[,9] <- 0 + (Samples_SIR[,7]-0)*Samples_SIR[,9] #beta_HA #beta_AA should always be bigger than beta_HA
    
    Samples_SIR[,10] <- 0 + (Samples_SIR[,5]-0)*Samples_SIR[,10] #beta_EA
    Samples_SIR[,8] <- 0 + (Samples_SIR[,4]-0)*Samples_SIR[,8]  #beta_EH
    
    Samples_SIR[,6] <- 0 + (1-0)*Samples_SIR[,6]  #beta_HE
    Samples_SIR[,11] <- 0 + (1-0)*Samples_SIR[,11]  #beta_AE
    
    Samples_SIR[,12] <- 0 + (0.5-0)*Samples_SIR[,12]  #mu_H #Carriage of ESBL or CRE at 12 months  community 25.4% patients 35.2%
    Samples_SIR[,13] <-  Samples_SIR[,12]#0 + (1-0)*Samples_SIR[,13]  #mu_A
    Samples_SIR[,14] <-  Samples_SIR[,12]  #0 + (1-0)*Samples_SIR[,14]    #mu_E
    Samples_SIR[,15] <- 0 + (1-0)*Samples_SIR[,15]    #gamma
    Samples_SIR[,16] <-  0 + (1-0)*Samples_SIR[,16] #beta_EE
    Samples_SIR[,17] <- 0 + (1-0)*Samples_SIR[,17] #epsilon
    
    
    paramsMat_SIR <- data.frame(
      LAMBDA_H=Samples_SIR[,1],LAMBDA_A=Samples_SIR[,2],LAMBDA_E=Samples_SIR[,3],
      beta_HH=Samples_SIR[,4],  beta_AA=Samples_SIR[,5],   beta_HA=Samples_SIR[,9], beta_AH=Samples_SIR[,7],
      beta_HE=Samples_SIR[,6],  beta_EH=Samples_SIR[,8],    beta_AE=Samples_SIR[,11],   beta_EA=Samples_SIR[,10],
      mu_H = Samples_SIR[,12], mu_A =Samples_SIR[,13], mu_E = Samples_SIR[,14],gamma = Samples_SIR[,15],beta_EE = Samples_SIR[,16],epsilon=Samples_SIR[,17]
    )
    return(paramsMat_SIR)
  }
  
  #AMR function
  state <- c(H=(1/1000),A=(1/1000),E=(1/1000)) 
  epid.start<- 2004
  epid.duration <- 50
  vectTime <- c(0,1:epid.duration)
  
  AMRmodel <- function(time,state,parameters){ #using package deSolve
    with(as.list(c(state,parameters)),{
      dH <- gamma*LAMBDA_H*(1-H) + LAMBDA_H*beta_HH*H*(1-H) + LAMBDA_H*beta_AH*(1-H)*A + LAMBDA_H*beta_EH*(1-H)*E - mu_H*H
      dA <-  gamma*LAMBDA_A*(1-A) + LAMBDA_A*beta_AA*A*(1-A) + LAMBDA_A*beta_HA*(1-A)*H + LAMBDA_A*beta_EA*(1-A)*E - mu_A*A
      dE <-  LAMBDA_E*beta_EE*E*(1-E) + LAMBDA_E*beta_HE*(1-E)*H + LAMBDA_E*beta_AE*(1-E)*A - mu_E*E
      return(  list(c(dH,dA,dE)))
    })}
  
  #Epid function
  epid <- function(LAMBDA_H, LAMBDA_A, LAMBDA_E,
                   beta_HH, beta_AA, 
                   beta_HA, beta_AH, 
                   beta_HE, beta_EH, 
                   beta_AE, beta_EA,
                   mu_H, mu_A, mu_E,gamma,returnout,beta_EE,epsilon
                   
  ){params <- c(LAMBDA_H=LAMBDA_H,LAMBDA_A=LAMBDA_A,LAMBDA_E=LAMBDA_E,
                beta_HH=beta_HH,  beta_AA=beta_AA,  beta_HE=beta_HE,  beta_AH=beta_AH,
                beta_EH=beta_EH,  beta_HA=beta_HA,  beta_EA=beta_EA,  beta_AE=beta_AE,  
                mu_H = mu_H, mu_A = mu_A, mu_E = mu_E,gamma=gamma,beta_EE=beta_EE,epsilon=epsilon)
  
  #run the model using desolve
  out <- as.data.frame(ode(y=state,time=vectTime,func=AMRmodel,parms=params))
  out$time = out$time+epid.start #rescale the time so that it runs from 2005 onwards 
  #Average compartment contents between 2012-2013
  model2012.H <- out$H[out$time==2012] #(out$H[out$time==2013] +out$H[out$time==2012])/2
  model2012.A <- out$A[out$time==2012] #(out$A[out$time==2013] +out$A[out$time==2012])/2
  model2012.E <- out$E[out$time==2012] #(out$E[out$time==2013] +out$E[out$time==2012])/2
  
  model2005.H <- out$H[out$time==2005] 
  model2005.A <- out$A[out$time==2005] 
  model2005.E <- out$E[out$time==2005] 
  #model2004.H <- out$H[out$time==2004] 
  #model2004.A <- out$A[out$time==2004] 
  #model2004.E <- out$E[out$time==2004] 
  
  model2008.H <- out$H[out$time==2008]
  model2009.H <- out$H[out$time==2009] 
  model2010.H <- out$H[out$time==2010]
  model2011.H <- out$H[out$time==2011]
  
  if (returnout ==1){
    return(out)} else{
      return(c(model2012.H= model2012.H,
               model2012.A= model2012.A,
               model2012.E= model2012.E,
               model2005.H=model2005.H,
               model2005.A=model2005.A,
               model2005.E=model2005.E,
               # model2004.H=model2004.H,
               # model2004.A=model2004.A,
               # model2004.E=model2004.E,
               model2008.H=model2008.H,
               model2009.H=model2009.H,
               model2010.H=model2010.H,
               model2011.H=model2011.H))}
  }
})
  
  save3 <- parallel::parLapply(cl,1:100, f)
  stopCluster(cl)

#Trying this with n=100 & variable # of clusters... 
#---------------------------------------------------------------------------------------------------------------------------------------------------------
system.time({
  cl <- makeCluster(3)
  clusterEvalQ(cl, {
    library(deSolve) 
    library(ggplot2) 
    library(tidyr)
    library(data.table)
    library(lhs)
    library(dplyr)
    
    f <- function(aaa){
      Samples_SIR <- randomLHS(aaa, 17) #from https://www.sciencedirect.com/science/article/pii/S2542519619301305
      
      Samples_SIR[,2] <- 0.1 + (1-0.1)*Samples_SIR[,2] #LAMBDA_A #vary to max
      Samples_SIR[,1] <- 0 + (Samples_SIR[,2]-0)*Samples_SIR[,1]  #LAMBDA_H #vary to max
      Samples_SIR[,3] <- 0 + (Samples_SIR[,1]-0)*Samples_SIR[,3] #LAMBDA_E #vary to max
      
      Samples_SIR[,4] <- 0.3 + (1-0.3)*Samples_SIR[,4]  #beta_HH
      Samples_SIR[,5] <- 0 + (1-0)*Samples_SIR[,5]  #beta_AA
      
      
      Samples_SIR[,7] <- 0 + (Samples_SIR[,4]-0)*Samples_SIR[,7] #beta_AH
      Samples_SIR[,9] <- 0 + (Samples_SIR[,7]-0)*Samples_SIR[,9] #beta_HA #beta_AA should always be bigger than beta_HA
      
      Samples_SIR[,10] <- 0 + (Samples_SIR[,5]-0)*Samples_SIR[,10] #beta_EA
      Samples_SIR[,8] <- 0 + (Samples_SIR[,4]-0)*Samples_SIR[,8]  #beta_EH
      
      Samples_SIR[,6] <- 0 + (1-0)*Samples_SIR[,6]  #beta_HE
      Samples_SIR[,11] <- 0 + (1-0)*Samples_SIR[,11]  #beta_AE
      
      Samples_SIR[,12] <- 0 + (0.5-0)*Samples_SIR[,12]  #mu_H #Carriage of ESBL or CRE at 12 months  community 25.4% patients 35.2%
      Samples_SIR[,13] <-  Samples_SIR[,12]#0 + (1-0)*Samples_SIR[,13]  #mu_A
      Samples_SIR[,14] <-  Samples_SIR[,12]  #0 + (1-0)*Samples_SIR[,14]    #mu_E
      Samples_SIR[,15] <- 0 + (1-0)*Samples_SIR[,15]    #gamma
      Samples_SIR[,16] <-  0 + (1-0)*Samples_SIR[,16] #beta_EE
      Samples_SIR[,17] <- 0 + (1-0)*Samples_SIR[,17] #epsilon
      
      
      paramsMat_SIR <- data.frame(
        LAMBDA_H=Samples_SIR[,1],LAMBDA_A=Samples_SIR[,2],LAMBDA_E=Samples_SIR[,3],
        beta_HH=Samples_SIR[,4],  beta_AA=Samples_SIR[,5],   beta_HA=Samples_SIR[,9], beta_AH=Samples_SIR[,7],
        beta_HE=Samples_SIR[,6],  beta_EH=Samples_SIR[,8],    beta_AE=Samples_SIR[,11],   beta_EA=Samples_SIR[,10],
        mu_H = Samples_SIR[,12], mu_A =Samples_SIR[,13], mu_E = Samples_SIR[,14],gamma = Samples_SIR[,15],beta_EE = Samples_SIR[,16],epsilon=Samples_SIR[,17]
      )
      return(paramsMat_SIR)
    }
    
    #AMR function
    state <- c(H=(1/1000),A=(1/1000),E=(1/1000)) 
    epid.start<- 2004
    epid.duration <- 50
    vectTime <- c(0,1:epid.duration)
    
    AMRmodel <- function(time,state,parameters){ #using package deSolve
      with(as.list(c(state,parameters)),{
        dH <- gamma*LAMBDA_H*(1-H) + LAMBDA_H*beta_HH*H*(1-H) + LAMBDA_H*beta_AH*(1-H)*A + LAMBDA_H*beta_EH*(1-H)*E - mu_H*H
        dA <-  gamma*LAMBDA_A*(1-A) + LAMBDA_A*beta_AA*A*(1-A) + LAMBDA_A*beta_HA*(1-A)*H + LAMBDA_A*beta_EA*(1-A)*E - mu_A*A
        dE <-  LAMBDA_E*beta_EE*E*(1-E) + LAMBDA_E*beta_HE*(1-E)*H + LAMBDA_E*beta_AE*(1-E)*A - mu_E*E
        return(  list(c(dH,dA,dE)))
      })}
    
    #Epid function
    epid <- function(LAMBDA_H, LAMBDA_A, LAMBDA_E,
                     beta_HH, beta_AA, 
                     beta_HA, beta_AH, 
                     beta_HE, beta_EH, 
                     beta_AE, beta_EA,
                     mu_H, mu_A, mu_E,gamma,returnout,beta_EE,epsilon
                     
    ){params <- c(LAMBDA_H=LAMBDA_H,LAMBDA_A=LAMBDA_A,LAMBDA_E=LAMBDA_E,
                  beta_HH=beta_HH,  beta_AA=beta_AA,  beta_HE=beta_HE,  beta_AH=beta_AH,
                  beta_EH=beta_EH,  beta_HA=beta_HA,  beta_EA=beta_EA,  beta_AE=beta_AE,  
                  mu_H = mu_H, mu_A = mu_A, mu_E = mu_E,gamma=gamma,beta_EE=beta_EE,epsilon=epsilon)
    
    #run the model using desolve
    out <- as.data.frame(ode(y=state,time=vectTime,func=AMRmodel,parms=params))
    out$time = out$time+epid.start #rescale the time so that it runs from 2005 onwards 
    #Average compartment contents between 2012-2013
    model2012.H <- out$H[out$time==2012] #(out$H[out$time==2013] +out$H[out$time==2012])/2
    model2012.A <- out$A[out$time==2012] #(out$A[out$time==2013] +out$A[out$time==2012])/2
    model2012.E <- out$E[out$time==2012] #(out$E[out$time==2013] +out$E[out$time==2012])/2
    
    model2005.H <- out$H[out$time==2005] 
    model2005.A <- out$A[out$time==2005] 
    model2005.E <- out$E[out$time==2005] 
    #model2004.H <- out$H[out$time==2004] 
    #model2004.A <- out$A[out$time==2004] 
    #model2004.E <- out$E[out$time==2004] 
    
    model2008.H <- out$H[out$time==2008]
    model2009.H <- out$H[out$time==2009] 
    model2010.H <- out$H[out$time==2010]
    model2011.H <- out$H[out$time==2011]
    
    if (returnout ==1){
      return(out)} else{
        return(c(model2012.H= model2012.H,
                 model2012.A= model2012.A,
                 model2012.E= model2012.E,
                 model2005.H=model2005.H,
                 model2005.A=model2005.A,
                 model2005.E=model2005.E,
                 # model2004.H=model2004.H,
                 # model2004.A=model2004.A,
                 # model2004.E=model2004.E,
                 model2008.H=model2008.H,
                 model2009.H=model2009.H,
                 model2010.H=model2010.H,
                 model2011.H=model2011.H))}
    }  
    
})
  save3 <- parallel::parLapply(cl,1:100, f)
  stopCluster(cl)
})
#With 1 cluster...
#user = 0.07
#system = 0.01
#elapsed = 1.36

#With 2 clusters...
#user = 0.03
#system = 0.04
#elapsed = 1.20

#With 3 clusters...
#user = 0.11
#system = 0.02
#elapsed = 1.99

#With 4 clusters... 
#user = 0.03
#system = 0.02
#elapsed = 1.67

#With 5 clusters... 
#user = 0.03
#system = 0.04
#elapsed = 1.97

save3 #Okay, this returns 100 rows of LAMBDA_H, LAMBDA_A, LAMBDA_E, beta_HH...

#I'm now writing the csv outside the parallel loop, but maybe...that's okay? 
write.csv(save3[[100]], paste0("OUTIE_", 100,".csv"))


#Trying with n=1000 & variable # of clusters... 
#---------------------------------------------------------------------------------------------------------------------------------------------------------
system.time({
  cl <- makeCluster(3)
  clusterEvalQ(cl, {
    library(deSolve) 
    library(ggplot2) 
    library(tidyr)
    library(data.table)
    library(lhs)
    library(dplyr)
    
    f <- function(aaa){
      Samples_SIR <- randomLHS(aaa, 17) #from https://www.sciencedirect.com/science/article/pii/S2542519619301305
      
      Samples_SIR[,2] <- 0.1 + (1-0.1)*Samples_SIR[,2] #LAMBDA_A #vary to max
      Samples_SIR[,1] <- 0 + (Samples_SIR[,2]-0)*Samples_SIR[,1]  #LAMBDA_H #vary to max
      Samples_SIR[,3] <- 0 + (Samples_SIR[,1]-0)*Samples_SIR[,3] #LAMBDA_E #vary to max
      
      Samples_SIR[,4] <- 0.3 + (1-0.3)*Samples_SIR[,4]  #beta_HH
      Samples_SIR[,5] <- 0 + (1-0)*Samples_SIR[,5]  #beta_AA
      
      
      Samples_SIR[,7] <- 0 + (Samples_SIR[,4]-0)*Samples_SIR[,7] #beta_AH
      Samples_SIR[,9] <- 0 + (Samples_SIR[,7]-0)*Samples_SIR[,9] #beta_HA #beta_AA should always be bigger than beta_HA
      
      Samples_SIR[,10] <- 0 + (Samples_SIR[,5]-0)*Samples_SIR[,10] #beta_EA
      Samples_SIR[,8] <- 0 + (Samples_SIR[,4]-0)*Samples_SIR[,8]  #beta_EH
      
      Samples_SIR[,6] <- 0 + (1-0)*Samples_SIR[,6]  #beta_HE
      Samples_SIR[,11] <- 0 + (1-0)*Samples_SIR[,11]  #beta_AE
      
      Samples_SIR[,12] <- 0 + (0.5-0)*Samples_SIR[,12]  #mu_H #Carriage of ESBL or CRE at 12 months  community 25.4% patients 35.2%
      Samples_SIR[,13] <-  Samples_SIR[,12]#0 + (1-0)*Samples_SIR[,13]  #mu_A
      Samples_SIR[,14] <-  Samples_SIR[,12]  #0 + (1-0)*Samples_SIR[,14]    #mu_E
      Samples_SIR[,15] <- 0 + (1-0)*Samples_SIR[,15]    #gamma
      Samples_SIR[,16] <-  0 + (1-0)*Samples_SIR[,16] #beta_EE
      Samples_SIR[,17] <- 0 + (1-0)*Samples_SIR[,17] #epsilon
      
      
      paramsMat_SIR <- data.frame(
        LAMBDA_H=Samples_SIR[,1],LAMBDA_A=Samples_SIR[,2],LAMBDA_E=Samples_SIR[,3],
        beta_HH=Samples_SIR[,4],  beta_AA=Samples_SIR[,5],   beta_HA=Samples_SIR[,9], beta_AH=Samples_SIR[,7],
        beta_HE=Samples_SIR[,6],  beta_EH=Samples_SIR[,8],    beta_AE=Samples_SIR[,11],   beta_EA=Samples_SIR[,10],
        mu_H = Samples_SIR[,12], mu_A =Samples_SIR[,13], mu_E = Samples_SIR[,14],gamma = Samples_SIR[,15],beta_EE = Samples_SIR[,16],epsilon=Samples_SIR[,17]
      )
      return(paramsMat_SIR)
    }
    
    #AMR function
    state <- c(H=(1/1000),A=(1/1000),E=(1/1000)) 
    epid.start<- 2004
    epid.duration <- 50
    vectTime <- c(0,1:epid.duration)
    
    AMRmodel <- function(time,state,parameters){ #using package deSolve
      with(as.list(c(state,parameters)),{
        dH <- gamma*LAMBDA_H*(1-H) + LAMBDA_H*beta_HH*H*(1-H) + LAMBDA_H*beta_AH*(1-H)*A + LAMBDA_H*beta_EH*(1-H)*E - mu_H*H
        dA <-  gamma*LAMBDA_A*(1-A) + LAMBDA_A*beta_AA*A*(1-A) + LAMBDA_A*beta_HA*(1-A)*H + LAMBDA_A*beta_EA*(1-A)*E - mu_A*A
        dE <-  LAMBDA_E*beta_EE*E*(1-E) + LAMBDA_E*beta_HE*(1-E)*H + LAMBDA_E*beta_AE*(1-E)*A - mu_E*E
        return(  list(c(dH,dA,dE)))
      })}
    
    #Epid function
    epid <- function(LAMBDA_H, LAMBDA_A, LAMBDA_E,
                     beta_HH, beta_AA, 
                     beta_HA, beta_AH, 
                     beta_HE, beta_EH, 
                     beta_AE, beta_EA,
                     mu_H, mu_A, mu_E,gamma,returnout,beta_EE,epsilon
                     
    ){params <- c(LAMBDA_H=LAMBDA_H,LAMBDA_A=LAMBDA_A,LAMBDA_E=LAMBDA_E,
                  beta_HH=beta_HH,  beta_AA=beta_AA,  beta_HE=beta_HE,  beta_AH=beta_AH,
                  beta_EH=beta_EH,  beta_HA=beta_HA,  beta_EA=beta_EA,  beta_AE=beta_AE,  
                  mu_H = mu_H, mu_A = mu_A, mu_E = mu_E,gamma=gamma,beta_EE=beta_EE,epsilon=epsilon)
    
    #run the model using desolve
    out <- as.data.frame(ode(y=state,time=vectTime,func=AMRmodel,parms=params))
    out$time = out$time+epid.start #rescale the time so that it runs from 2005 onwards 
    #Average compartment contents between 2012-2013
    model2012.H <- out$H[out$time==2012] #(out$H[out$time==2013] +out$H[out$time==2012])/2
    model2012.A <- out$A[out$time==2012] #(out$A[out$time==2013] +out$A[out$time==2012])/2
    model2012.E <- out$E[out$time==2012] #(out$E[out$time==2013] +out$E[out$time==2012])/2
    
    model2005.H <- out$H[out$time==2005] 
    model2005.A <- out$A[out$time==2005] 
    model2005.E <- out$E[out$time==2005] 
    #model2004.H <- out$H[out$time==2004] 
    #model2004.A <- out$A[out$time==2004] 
    #model2004.E <- out$E[out$time==2004] 
    
    model2008.H <- out$H[out$time==2008]
    model2009.H <- out$H[out$time==2009] 
    model2010.H <- out$H[out$time==2010]
    model2011.H <- out$H[out$time==2011]
    
    if (returnout ==1){
      return(out)} else{
        return(c(model2012.H= model2012.H,
                 model2012.A= model2012.A,
                 model2012.E= model2012.E,
                 model2005.H=model2005.H,
                 model2005.A=model2005.A,
                 model2005.E=model2005.E,
                 # model2004.H=model2004.H,
                 # model2004.A=model2004.A,
                 # model2004.E=model2004.E,
                 model2008.H=model2008.H,
                 model2009.H=model2009.H,
                 model2010.H=model2010.H,
                 model2011.H=model2011.H))}
    }  
    
  })
  save_t <- parallel::parLapply(cl,1:1000, f)
  stopCluster(cl)
})
#With 1 cluster...
#user = 0.17
#system = 0.11
#elapsed = 3.86

#With 2 clusters...
#user = 0.17
#system = 0.15
#elapsed = 4.30

#With 3 clusters...
#user = 0.22
#system = 0.08
#elapsed = 3.64

#With 4 clusters... 
#user = 0.14
#system = 0.16
#elapsed = 2.78

#With 5 clusters... 
#user = 0.11
#system = 0.11
#elapsed = 2.95

write.csv(save_t[[1000]], paste0("OUTIE_", 1000,".csv"))


#Trying with n=10k & 3 clusters... 
#---------------------------------------------------------------------------------------------------------------------------------------------------------
system.time({
  cl <- makeCluster(3)
  clusterEvalQ(cl, {
    library(deSolve) 
    library(ggplot2) 
    library(tidyr)
    library(data.table)
    library(lhs)
    library(dplyr)
    
    #Sampling function
    f <- function(aaa){
      Samples_SIR <- randomLHS(aaa, 17) #from https://www.sciencedirect.com/science/article/pii/S2542519619301305
      
      Samples_SIR[,2] <- 0.1 + (1-0.1)*Samples_SIR[,2] #LAMBDA_A #vary to max
      Samples_SIR[,1] <- 0 + (Samples_SIR[,2]-0)*Samples_SIR[,1]  #LAMBDA_H #vary to max
      Samples_SIR[,3] <- 0 + (Samples_SIR[,1]-0)*Samples_SIR[,3] #LAMBDA_E #vary to max
      
      Samples_SIR[,4] <- 0.3 + (1-0.3)*Samples_SIR[,4]  #beta_HH
      Samples_SIR[,5] <- 0 + (1-0)*Samples_SIR[,5]  #beta_AA
      
      
      Samples_SIR[,7] <- 0 + (Samples_SIR[,4]-0)*Samples_SIR[,7] #beta_AH
      Samples_SIR[,9] <- 0 + (Samples_SIR[,7]-0)*Samples_SIR[,9] #beta_HA #beta_AA should always be bigger than beta_HA
      
      Samples_SIR[,10] <- 0 + (Samples_SIR[,5]-0)*Samples_SIR[,10] #beta_EA
      Samples_SIR[,8] <- 0 + (Samples_SIR[,4]-0)*Samples_SIR[,8]  #beta_EH
      
      Samples_SIR[,6] <- 0 + (1-0)*Samples_SIR[,6]  #beta_HE
      Samples_SIR[,11] <- 0 + (1-0)*Samples_SIR[,11]  #beta_AE
      
      Samples_SIR[,12] <- 0 + (0.5-0)*Samples_SIR[,12]  #mu_H #Carriage of ESBL or CRE at 12 months  community 25.4% patients 35.2%
      Samples_SIR[,13] <-  Samples_SIR[,12]#0 + (1-0)*Samples_SIR[,13]  #mu_A
      Samples_SIR[,14] <-  Samples_SIR[,12]  #0 + (1-0)*Samples_SIR[,14]    #mu_E
      Samples_SIR[,15] <- 0 + (1-0)*Samples_SIR[,15]    #gamma
      Samples_SIR[,16] <-  0 + (1-0)*Samples_SIR[,16] #beta_EE
      Samples_SIR[,17] <- 0 + (1-0)*Samples_SIR[,17] #epsilon
      
      
      paramsMat_SIR <- data.frame(
        LAMBDA_H=Samples_SIR[,1],LAMBDA_A=Samples_SIR[,2],LAMBDA_E=Samples_SIR[,3],
        beta_HH=Samples_SIR[,4],  beta_AA=Samples_SIR[,5],   beta_HA=Samples_SIR[,9], beta_AH=Samples_SIR[,7],
        beta_HE=Samples_SIR[,6],  beta_EH=Samples_SIR[,8],    beta_AE=Samples_SIR[,11],   beta_EA=Samples_SIR[,10],
        mu_H = Samples_SIR[,12], mu_A =Samples_SIR[,13], mu_E = Samples_SIR[,14],gamma = Samples_SIR[,15],beta_EE = Samples_SIR[,16],epsilon=Samples_SIR[,17]
      )
      return(paramsMat_SIR)
    }
    
    #AMR function
    state <- c(H=(1/1000),A=(1/1000),E=(1/1000)) 
    epid.start<- 2004
    epid.duration <- 50
    vectTime <- c(0,1:epid.duration)
    
    AMRmodel <- function(time,state,parameters){ #using package deSolve
      with(as.list(c(state,parameters)),{
        dH <- gamma*LAMBDA_H*(1-H) + LAMBDA_H*beta_HH*H*(1-H) + LAMBDA_H*beta_AH*(1-H)*A + LAMBDA_H*beta_EH*(1-H)*E - mu_H*H
        dA <-  gamma*LAMBDA_A*(1-A) + LAMBDA_A*beta_AA*A*(1-A) + LAMBDA_A*beta_HA*(1-A)*H + LAMBDA_A*beta_EA*(1-A)*E - mu_A*A
        dE <-  LAMBDA_E*beta_EE*E*(1-E) + LAMBDA_E*beta_HE*(1-E)*H + LAMBDA_E*beta_AE*(1-E)*A - mu_E*E
        return(  list(c(dH,dA,dE)))
      })}
    
    #Epid function
    epid <- function(LAMBDA_H, LAMBDA_A, LAMBDA_E,
                     beta_HH, beta_AA, 
                     beta_HA, beta_AH, 
                     beta_HE, beta_EH, 
                     beta_AE, beta_EA,
                     mu_H, mu_A, mu_E,gamma,returnout,beta_EE,epsilon
                     
    ){params <- c(LAMBDA_H=LAMBDA_H,LAMBDA_A=LAMBDA_A,LAMBDA_E=LAMBDA_E,
                  beta_HH=beta_HH,  beta_AA=beta_AA,  beta_HE=beta_HE,  beta_AH=beta_AH,
                  beta_EH=beta_EH,  beta_HA=beta_HA,  beta_EA=beta_EA,  beta_AE=beta_AE,  
                  mu_H = mu_H, mu_A = mu_A, mu_E = mu_E,gamma=gamma,beta_EE=beta_EE,epsilon=epsilon)
    
    
    #run the model using desolve
    out <- as.data.frame(ode(y=state,time=vectTime,func=AMRmodel,parms=params))
    out$time = out$time+epid.start #rescale the time so that it runs from 2005 onwards 
    #Average compartment contents between 2012-2013
    model2012.H <- out$H[out$time==2012] #(out$H[out$time==2013] +out$H[out$time==2012])/2
    model2012.A <- out$A[out$time==2012] #(out$A[out$time==2013] +out$A[out$time==2012])/2
    model2012.E <- out$E[out$time==2012] #(out$E[out$time==2013] +out$E[out$time==2012])/2
    
    model2005.H <- out$H[out$time==2005] 
    model2005.A <- out$A[out$time==2005] 
    model2005.E <- out$E[out$time==2005] 
    #model2004.H <- out$H[out$time==2004] 
    #model2004.A <- out$A[out$time==2004] 
    #model2004.E <- out$E[out$time==2004] 
    
    model2008.H <- out$H[out$time==2008]
    model2009.H <- out$H[out$time==2009] 
    model2010.H <- out$H[out$time==2010]
    model2011.H <- out$H[out$time==2011]
    
    
    if (returnout ==1){
      return(out)} else{
        
        return(c(model2012.H= model2012.H,
                 model2012.A= model2012.A,
                 model2012.E= model2012.E,
                 model2005.H=model2005.H,
                 model2005.A=model2005.A,
                 model2005.E=model2005.E,
                 # model2004.H=model2004.H,
                 # model2004.A=model2004.A,
                 # model2004.E=model2004.E,
                 model2008.H=model2008.H,
                 model2009.H=model2009.H,
                 model2010.H=model2010.H,
                 model2011.H=model2011.H))}
    }
  })
  save_tt <- parLapply(cl,1:10000, f)
  stopCluster(cl)
})
#With 3 clusters...
#user = 6.96
#system = 4.97
#elapsed = 161.72

#write.csv(save_tt[[10000]], paste0("OUTIE_", 10000,".csv"))


#Trying with n=100k & 3 clusters...
#---------------------------------------------------------------------------------------------------------------------------------------------------------
system.time({
  cl <- makeCluster(3)
  clusterEvalQ(cl, {
    library(deSolve) 
    library(ggplot2) 
    library(tidyr)
    library(data.table)
    library(lhs)
    library(dplyr)
    
    #Sampling function
    f <- function(aaa){
      Samples_SIR <- randomLHS(aaa, 17) #from https://www.sciencedirect.com/science/article/pii/S2542519619301305
      
      Samples_SIR[,2] <- 0.1 + (1-0.1)*Samples_SIR[,2] #LAMBDA_A #vary to max
      Samples_SIR[,1] <- 0 + (Samples_SIR[,2]-0)*Samples_SIR[,1]  #LAMBDA_H #vary to max
      Samples_SIR[,3] <- 0 + (Samples_SIR[,1]-0)*Samples_SIR[,3] #LAMBDA_E #vary to max
      
      Samples_SIR[,4] <- 0.3 + (1-0.3)*Samples_SIR[,4]  #beta_HH
      Samples_SIR[,5] <- 0 + (1-0)*Samples_SIR[,5]  #beta_AA
      
      
      Samples_SIR[,7] <- 0 + (Samples_SIR[,4]-0)*Samples_SIR[,7] #beta_AH
      Samples_SIR[,9] <- 0 + (Samples_SIR[,7]-0)*Samples_SIR[,9] #beta_HA #beta_AA should always be bigger than beta_HA
      
      Samples_SIR[,10] <- 0 + (Samples_SIR[,5]-0)*Samples_SIR[,10] #beta_EA
      Samples_SIR[,8] <- 0 + (Samples_SIR[,4]-0)*Samples_SIR[,8]  #beta_EH
      
      Samples_SIR[,6] <- 0 + (1-0)*Samples_SIR[,6]  #beta_HE
      Samples_SIR[,11] <- 0 + (1-0)*Samples_SIR[,11]  #beta_AE
      
      Samples_SIR[,12] <- 0 + (0.5-0)*Samples_SIR[,12]  #mu_H #Carriage of ESBL or CRE at 12 months  community 25.4% patients 35.2%
      Samples_SIR[,13] <-  Samples_SIR[,12]#0 + (1-0)*Samples_SIR[,13]  #mu_A
      Samples_SIR[,14] <-  Samples_SIR[,12]  #0 + (1-0)*Samples_SIR[,14]    #mu_E
      Samples_SIR[,15] <- 0 + (1-0)*Samples_SIR[,15]    #gamma
      Samples_SIR[,16] <-  0 + (1-0)*Samples_SIR[,16] #beta_EE
      Samples_SIR[,17] <- 0 + (1-0)*Samples_SIR[,17] #epsilon
      
      
      paramsMat_SIR <- data.frame(
        LAMBDA_H=Samples_SIR[,1],LAMBDA_A=Samples_SIR[,2],LAMBDA_E=Samples_SIR[,3],
        beta_HH=Samples_SIR[,4],  beta_AA=Samples_SIR[,5],   beta_HA=Samples_SIR[,9], beta_AH=Samples_SIR[,7],
        beta_HE=Samples_SIR[,6],  beta_EH=Samples_SIR[,8],    beta_AE=Samples_SIR[,11],   beta_EA=Samples_SIR[,10],
        mu_H = Samples_SIR[,12], mu_A =Samples_SIR[,13], mu_E = Samples_SIR[,14],gamma = Samples_SIR[,15],beta_EE = Samples_SIR[,16],epsilon=Samples_SIR[,17]
      )
      return(paramsMat_SIR)
    }
    
    #AMR function
    state <- c(H=(1/1000),A=(1/1000),E=(1/1000)) 
    epid.start<- 2004
    epid.duration <- 50
    vectTime <- c(0,1:epid.duration)
    
    AMRmodel <- function(time,state,parameters){ #using package deSolve
      with(as.list(c(state,parameters)),{
        dH <- gamma*LAMBDA_H*(1-H) + LAMBDA_H*beta_HH*H*(1-H) + LAMBDA_H*beta_AH*(1-H)*A + LAMBDA_H*beta_EH*(1-H)*E - mu_H*H
        dA <-  gamma*LAMBDA_A*(1-A) + LAMBDA_A*beta_AA*A*(1-A) + LAMBDA_A*beta_HA*(1-A)*H + LAMBDA_A*beta_EA*(1-A)*E - mu_A*A
        dE <-  LAMBDA_E*beta_EE*E*(1-E) + LAMBDA_E*beta_HE*(1-E)*H + LAMBDA_E*beta_AE*(1-E)*A - mu_E*E
        return(  list(c(dH,dA,dE)))
      })}
    
    #Epid function
    epid <- function(LAMBDA_H, LAMBDA_A, LAMBDA_E,
                     beta_HH, beta_AA, 
                     beta_HA, beta_AH, 
                     beta_HE, beta_EH, 
                     beta_AE, beta_EA,
                     mu_H, mu_A, mu_E,gamma,returnout,beta_EE,epsilon
                     
    ){params <- c(LAMBDA_H=LAMBDA_H,LAMBDA_A=LAMBDA_A,LAMBDA_E=LAMBDA_E,
                  beta_HH=beta_HH,  beta_AA=beta_AA,  beta_HE=beta_HE,  beta_AH=beta_AH,
                  beta_EH=beta_EH,  beta_HA=beta_HA,  beta_EA=beta_EA,  beta_AE=beta_AE,  
                  mu_H = mu_H, mu_A = mu_A, mu_E = mu_E,gamma=gamma,beta_EE=beta_EE,epsilon=epsilon)
    
    
    #run the model using desolve
    out <- as.data.frame(ode(y=state,time=vectTime,func=AMRmodel,parms=params))
    out$time = out$time+epid.start #rescale the time so that it runs from 2005 onwards 
    #Average compartment contents between 2012-2013
    model2012.H <- out$H[out$time==2012] #(out$H[out$time==2013] +out$H[out$time==2012])/2
    model2012.A <- out$A[out$time==2012] #(out$A[out$time==2013] +out$A[out$time==2012])/2
    model2012.E <- out$E[out$time==2012] #(out$E[out$time==2013] +out$E[out$time==2012])/2
    
    model2005.H <- out$H[out$time==2005] 
    model2005.A <- out$A[out$time==2005] 
    model2005.E <- out$E[out$time==2005] 
    #model2004.H <- out$H[out$time==2004] 
    #model2004.A <- out$A[out$time==2004] 
    #model2004.E <- out$E[out$time==2004] 
    
    model2008.H <- out$H[out$time==2008]
    model2009.H <- out$H[out$time==2009] 
    model2010.H <- out$H[out$time==2010]
    model2011.H <- out$H[out$time==2011]
    
    
    if (returnout ==1){
      return(out)} else{
        
        return(c(model2012.H= model2012.H,
                 model2012.A= model2012.A,
                 model2012.E= model2012.E,
                 model2005.H=model2005.H,
                 model2005.A=model2005.A,
                 model2005.E=model2005.E,
                 # model2004.H=model2004.H,
                 # model2004.A=model2004.A,
                 # model2004.E=model2004.E,
                 model2008.H=model2008.H,
                 model2009.H=model2009.H,
                 model2010.H=model2010.H,
                 model2011.H=model2011.H))}
    }
  })
  save_ht <- parLapply(cl,1:100000, f)
  stopCluster(cl)
})
##################### ATTEMPT 1 ###############################################################################################################################
#Error in checkForRemoteErros(val): 
#3 nodes produced errors; first error: std::bad_alloc
### Unfortunately, the internet suggests this is an issue with my R...
#Timing stopped at...
#user = 1.22
#system = 0.7
#elapsed = 983

#write.csv(save_ht[[100000]], paste0("OUTIE_", 100000,".csv"))