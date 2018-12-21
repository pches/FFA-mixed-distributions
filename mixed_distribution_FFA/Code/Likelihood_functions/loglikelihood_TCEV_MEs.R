# Author: Kenneth Joel Roop-Eckart
# Date: 8/22/2017
# Purpose: likelihood function for Two Component Extreme Value (TCEV) distribution
# used for fitting TCEV
# Copyright 2017
# uses the TCEV parameterization from: Rossi, Fabio, Mauro Fiorentino, and Pasquale Versace. 1984. "Two-Component Extreme Value Distribution for Flood Frequency Analysis." Water Resources Research 20 (7): 847-56. doi:10.1029/WR020i007p00847.

# TCEV gage
densityTCEV_gage_MEs <- function(parameters,data,gage_uncertainty) {
  
  A1 <- parameters[1]
  theta1 <- parameters[2]
  A2 <- parameters[3]
  theta2 <- parameters[4]
  like <- matrix(data = 0, ncol = length(data), nrow = 1)
  I = length(data)
  if(theta1==0)return(Inf)
  else if(theta2==0)return(Inf)
  else if(A1==0)return(Inf)
  else if(A2==0)return(Inf)
  else {
    for(i in 1:I){
      Xij <- qnorm(seq(from = 0.025, to = 0.975, by = 0.095), mean = data[i], sd = gage_uncertainty[i])
      Fij <- dnorm(Xij, mean = data[i], sd = gage_uncertainty[i])
      Fij <- Fij/sum(Fij)
      if(min(Xij)<0)(return(-Inf))
      else(
        like[1,i] <- sum(
          Fij*(exp(A1*(-exp(-Xij/theta1))-A2*exp(-Xij/theta2))*
                 ((A1*exp(-Xij/theta1))/theta1+(A2*exp(-Xij/theta2))/theta2))
        )
      )
    }
  }
  return(like)
}

####################################### TCEV historical component #############################
densityTCEV_hist_MEs <- function(parameters,hist,hist_uncertainty) {
  
  A1 <- parameters[1]
  theta1 <- parameters[2]
  A2 <- parameters[3]
  theta2 <- parameters[4]
  like <- matrix(data = 0, ncol = length(hist), nrow = 1)
  I = length(hist)
  if(theta1==0)return(Inf)
  else if(theta2==0)return(Inf)
  else if(A1==0)return(Inf)
  else if(A2==0)return(Inf)
  else {
    for(i in 1:I){
      Xij <- qnorm(seq(from = 0.025, to = 0.975, by = 0.095), mean = hist[i], sd = hist_uncertainty[i])
      Fij <- dnorm(Xij, mean = hist[i], sd = hist_uncertainty[i])
      Fij <- Fij/sum(Fij)
      if(min(Xij)<0)(return(-Inf))
      else(
        like[1,i] <- sum(
          Fij*(exp(A1*(-exp(-Xij/theta1))-A2*exp(-Xij/theta2))*
                 ((A1*exp(-Xij/theta1))/theta1+(A2*exp(-Xij/theta2))/theta2))
        )
      )
    }
    out <- like*(dbinom(length(hist), size = hist_record, prob = 1-exp(-A1*exp(-min(hist)/theta1)-A2*exp(-min(hist)/theta2))))/
      (1-exp(-A1*exp(-min(hist)/theta1)-A2*exp(-min(hist)/theta2)))
  }
  return(out)
}

############################################ TCEV paleo component ###################################

densityTCEV_paleo_MEs <- function(parameters, data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty) {
  
  A1 <- parameters[1]
  theta1 <- parameters[2]
  A2 <- parameters[3]
  theta2 <- parameters[4]
  like <- vector(mode = 'numeric', length = length(paleo))
  I = length(paleo_age)
  Xmin <- paleo_flow_lower
  if(theta1==0)return(Inf)
  else if(theta2==0)return(Inf)
  else if(A1==0)return(Inf)
  else if(A2==0)return(Inf)
  else {
    #    for(i in 1:I){
    #      # time component of the likelihood function
    #      Nij <- paleo_age_uncertainty[,1] # dstribution of ages
    #      FUnij <- paleo_age_uncertainty[,2]
    #      # flood size component of the likelihood function
    #      variable <- paleo_flow_uncertainty[,1] # distribution of flows based on a normal distribution around a best estimate
    #      Uik <- (variable)
    #      FUdik <- paleo_flow_uncertainty[,2]
    #      intvect <- exp(-A1*exp(-Uik/theta1)-A2*exp(-Uik/theta2))-exp(-A1*exp(-Xmin/theta1)-A2*exp(-Xmin/theta2)) # the probability space between the upper and lower estimates
    #      like[i] <- sum(FUnij*Nij)*#log(
    #        sum(FUdik*intvect # multiply the discrete probability distribution by the probability space between the paleo upper and lower bounds
    #        )#)
    for(i in 1:I){
      FUdik <- paleo_flow_uncertainty[,2]
      variable <- paleo_flow_uncertainty[,1]
      Uik <- (variable)
      if(parameters[1] < 0)(return(rep(0,length(paleo))))
      else if(parameters[2] < 0)(return(rep(0,length(paleo))))
      else{
        like[i] <- sum(
          FUdik*(exp(-A1*exp(-Uik/theta1)-A2*exp(-Uik/theta2))-exp(-A1*exp(-Xmin/theta1)-A2*exp(-Xmin/theta2)))
        )
      }
    }
    probs <- (paleo_flow_uncertainty[,2]*paleo_age_uncertainty[,2])/(sum(paleo_flow_uncertainty[,2]*paleo_age_uncertainty[,2]))
    
    ############### Addition
    age_like <- vector(mode = 'numeric', length = I)
    for(i in 1:I){
      age_like[i] <- dbinom(length(paleo_age[i]), size = round(paleo_age[i]), prob = (1-pgev1(paleo_flow_lower, xi = parameters[1],  beta = parameters[2], mu = parameters[3])))/
        (1-pgev1(paleo_flow_lower, xi = parameters[1],  beta = parameters[2], mu = parameters[3]))
    }
    if(is.na(min(age_like))==TRUE)(which(is.nan(age_like))==TRUE)
    
    
    ############## Addition
    out <- vector(mode = 'numeric', length = length(age_like))
    for (i in 1:length(age_like)) {
      out[i] <- sum(like*probs*age_like[i])
    }
    out <- sum(out)
    ############## End addition
    
    
#    out <- (like*(dbinom(length(paleo_age), size = paleo_age, prob = 1-exp(-A1*exp(-paleo_flow_lower/theta1)-A2*exp(-paleo_flow_lower/theta2))))/
#              (1-exp(-A1*exp(-paleo_flow_lower/theta1)-A2*exp(-paleo_flow_lower/theta2))))
  }
  return(out)
}
############################### gage hist paleo TCEV with MEs ###########################
TCEV_gage_hist_paleo <- function(parameters, data, gage_uncertainty, hist, hist_record, hist_uncertainty, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty){
  A1 <- parameters[1]
  theta1 <- parameters[2]
  A2 <- parameters[3]
  theta2 <- parameters[4]
  density <- c(densityTCEV_gage_MEs(parameters, data, gage_uncertainty),
               densityTCEV_hist_MEs(parameters, hist, hist_uncertainty),
               densityTCEV_paleo_MEs(parameters, data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty))
  loglik <- sum(log(density))
  if(is.na(loglik)==TRUE)(loglik=-Inf)
  return(loglik)
}


TCEV_gage_hist_paleo_MLE <- function(parameters, data, gage_uncertainty, hist, hist_record, hist_uncertainty, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty){
  A1 <- parameters[1]
  theta1 <- parameters[2]
  A2 <- parameters[3]
  theta2 <- parameters[4]
  density <- c(densityTCEV_gage_MEs(parameters, data, gage_uncertainty),
               densityTCEV_hist_MEs(parameters, hist, hist_uncertainty),
               densityTCEV_paleo_MEs(parameters, data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty))
  loglik <- sum(log(density))
  if(is.na(loglik)==TRUE)(loglik=-Inf)
  return(-loglik)
}


TCEV_gage_hist_paleo_density <- function(parameters, data, gage_uncertainty, hist, hist_record, hist_uncertainty, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty){
  A1 <- parameters[1]
  theta1 <- parameters[2]
  A2 <- parameters[3]
  theta2 <- parameters[4]
  density <- c(densityTCEV_gage_MEs(parameters, data, gage_uncertainty),
               densityTCEV_hist_MEs(parameters, hist, hist_uncertainty),
               densityTCEV_paleo_MEs(parameters, data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty))
  #  loglik <- sum(log(density))
  #  if(is.na(loglik)==TRUE)(loglik=-Inf)
  return(density)
}