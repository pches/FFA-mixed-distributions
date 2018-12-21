

# created by Joel Roop-Eckart and Ben Lee
# papers used: (Stedinger et al., 1986) and (O'Connell, 2002)
# Stedinger, Jery R., and Timothy A. Cohn. 1986. "Flood Frequency Analysis With Historical and Paleoflood Information." Water Resources Research 22 (5): 785-93. doi:10.1029/WR022i005p00785.
# O'Connell, Daniel R. H. 2002. "Bayesian Flood Frequency Analysis with Paleohydrologic Bound Data." Water Resources Research 38 (5). doi:10.1029/2000WR000028.

require(zipfR)

LP3prob <- function(quantiles, parameters){
  tau <- parameters[1]
  alpha <- parameters[2]
  beta <- parameters[3]
  x <- log10(quantiles)
  out <- (beta*Igamma(alpha, (x-tau)/beta))/(gamma(alpha)*abs(beta))
  return(out)
}
################################### Gage Component ###################################
densityMLP3_gage<-function(parameters,data,gage_uncertainty){
  #  if(parameters[1] < log10(min(data)))(return(rep(0,length(data))))
  #  else if(parameters[2] < 0)(return(rep(0,length(data))))
  #  else{
  like <- matrix(data = 0, ncol = length(data), nrow = 1)
  I = length(data)
  for(i in 1:I){
    Xij <- qnorm(seq(from = 0.025, to = 0.975, by = 0.095), mean = data[i], sd = gage_uncertainty[i])
    Fij <- dnorm(Xij, mean = data[i], sd = gage_uncertainty[i])
    Fij <- Fij/sum(Fij)
    if(parameters[1] < 0)(return(rep(0,length(data))))
    else if(parameters[2] < 0)(return(rep(0,length(data))))
    else{
      like[1,i] <- sum(
        Fij*(((((log10(Xij)-parameters[1])/parameters[3])^(parameters[2]-1))*exp(-(log10(Xij)-parameters[1])/parameters[3]))/(abs(parameters[3])*gamma(parameters[2])))*(1/(Xij*log(10)))
      )
    }
  }
  return(like)
  #  }
}
###################################### Historical Component ####################################
densityMLP3_hist<-function(parameters,hist,hist_uncertainty){
  #  if(parameters[1] < log10(min(hist)))(return(rep(0,length(hist))))
  #  else if(parameters[2] < 0)(return(rep(0,length(hist))))
  #  else{
  like <- matrix(data = 0, ncol = length(hist), nrow = 1)
  I = length(hist)
  for(i in 1:I){
    Xij <- qnorm(seq(from = 0.025, to = 0.975, by = 0.095), mean = hist[i], sd = hist_uncertainty[i])
    Fij <- dnorm(Xij, mean = hist[i], sd = hist_uncertainty[i])
    Fij <- Fij/sum(Fij)
    if(parameters[1] < 0)(return(rep(0,length(hist))))
    else if(parameters[2] < 0)(return(rep(0,length(hist))))
    else{
      like[1,i] <- sum(
        Fij*(((((log10(Xij)-parameters[1])/parameters[3])^(parameters[2]-1))*exp(-(log10(Xij)-parameters[1])/parameters[3]))/(abs(parameters[3])*gamma(parameters[2])))*(1/(Xij*log(10)))
      )
    }
  }
  out <- (like*dbinom(length(hist), size = hist_record, prob = 1-LP3prob(min(hist),parameters))/
            (1-LP3prob(min(hist), parameters)))
  #  }
  return(out)
}
########################################### Paleo Component ################################
densityMLP3_paleo <- function(parameters, data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty){
  I <- length(paleo_age)
  paleo_age <- paleo_age
  variable <- paleo_flow_uncertainty[,1]
  Uik <- (variable)
  Xmin <- paleo_flow_lower
  parameters <- parameters
  like <- vector(mode = 'numeric', length = length(paleo_age))
  #  for(i in 1:I){
  #    # time component of the likelihood function
  #    Nij <- paleo_age_uncertainty[,1] # dstribution of ages
  #    FUnij <- paleo_age_uncertainty[,2]
  #    # flood size component of the likelihood function
  #    variable <- paleo_flow_uncertainty[,1] # distribution of flows based on a normal distribution around a best estimate
  #    Uik <- (variable)
  #    FUdik <- paleo_flow_uncertainty[,2]
  #    intvect <- pgengamma(log(Uik), mu = parameters[1], sigma = parameters[2], Q = parameters[3])-
  #      pgengamma(log(Xmin), mu = parameters[1], sigma = parameters[2], Q = parameters[3]) # the probability space between the upper and lower estimates
  #    
  #    like[i] <- sum(FUnij*Nij)*#log(
  #      sum(FUdik*intvect # multiply the discrete probability distribution by the probability space between the paleo upper and lower bounds
  #      )#)
  #  }
  
  Uik <- paleo_flow_uncertainty[,1]
  paleo_age <- paleo_age_uncertainty[,1]
  I <- length(paleo_age)
  probs <- (paleo_flow_uncertainty[,2]*paleo_age_uncertainty[,2])/(sum(paleo_flow_uncertainty[,2]*paleo_age_uncertainty[,2]))
  
  like <- (((((log10(Uik)-parameters[1])/parameters[3])^(parameters[2]-1))*exp(-(log10(Uik)-parameters[1])/parameters[3]))/(abs(parameters[3])*gamma(parameters[2])))*(1/(Uik*log(10)))
  
  age_like <- vector(mode = 'numeric', length = I)
  for(i in 1:I){
    age_like[i] <- dbinom(length(paleo_age[i]), size = round(paleo_age[i]), prob = (1-LP3prob(paleo_flow_lower, parameters)))/
      (1-LP3prob(paleo_flow_lower, parameters))
  }
  if(is.na(min(age_like))==TRUE)(which(is.nan(age_like))==TRUE)
  
#  out <- sum(like*age_like*probs)
  
  ############## Addition
  out <- vector(mode = 'numeric', length = length(age_like))
  for (i in 1:length(age_like)) {
    out[i] <- sum(like*probs*age_like[i])
  }
  out <- sum(out)
  ############## End addition
  
  
  #  for(i in 1:I){
  #    FUdik <- paleo_flow_uncertainty[,2]
  #    variable <- paleo_flow_uncertainty[,1]
  #    Uik <- (variable)
  #    if(parameters[1] < log10(min(Uik)))(return(rep(0,length(paleo))))
  #    else if(parameters[2] < 0)(return(rep(0,length(paleo))))
  #    else{
  #      like <- sum(
  #        FUdik*(((((log10(Uik)-parameters[1])/parameters[3])^(parameters[2]-1))*exp(-(log10(Uik)-parameters[1])/parameters[3]))/(abs(parameters[3])*gamma(parameters[2])))*(1/(Uik*log(10)))
  #      )
  #    }
  #  }
  
  
  #  out <- (like*dbinom(length(800), size = 800, prob = (1-LP3prob(paleo_flow_lower, parameters)))/
  #            (1-LP3prob(paleo_flow_lower, parameters)))
  return(out)
}
####################################### Complete Likelihood Function ##############################
LP3_gage_hist_paleo_mixed_test <- function(parameters, data, gage_uncertainty, hist, hist_record, hist_uncertainty, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty){
  density <- parameters[7]*c(densityMLP3_gage(parameters[1:3], data, gage_uncertainty),
                             densityMLP3_hist(parameters[1:3], hist, hist_uncertainty),
                             densityMLP3_paleo(parameters[1:3], data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty))+
    (1-parameters[7])*c(densityMLP3_gage(parameters[4:6], data, gage_uncertainty),
                        densityMLP3_hist(parameters[4:6], hist, hist_uncertainty),
                        densityMLP3_paleo(parameters[4:6], data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty))
  loglik <- sum(log(density))
  if(is.na(loglik)==TRUE)(loglik=-Inf)
  return(loglik)
}

# For MLE (positive log like values)

LP3_gage_hist_paleo_mixed_MLE_test <- function(parameters, data, gage_uncertainty, hist, hist_record, hist_uncertainty, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty){
  density <- parameters[7]*c(densityMLP3_gage(parameters[1:3], data, gage_uncertainty),
                             densityMLP3_hist(parameters[1:3], hist, hist_uncertainty),
                             densityMLP3_paleo(parameters[1:3], data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty))+
    (1-parameters[7])*c(densityMLP3_gage(parameters[4:6], data, gage_uncertainty),
                        densityMLP3_hist(parameters[4:6], hist, hist_uncertainty),
                        densityMLP3_paleo(parameters[4:6], data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty))
  loglik <- sum(log(density))
  if(is.na(loglik)==TRUE)(loglik=-Inf)
  #  if(parameters[7] > 0.999)(loglik=-Inf)
  #  if(parameters[7] < 0.001)(loglik=-Inf)
  #  loglik<-(log(dnorm(parameters[5], mean = (1-0.245749),sd = 0.047778)))+loglik
  return(-loglik)
}

# density functions

densityMLP3 <- function(parameters, data, gage_uncertainty, hist, hist_record, hist_uncertainty, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty){
  density <- parameters[7]*c(densityMLP3_gage(parameters[1:3], data, gage_uncertainty),
                             densityMLP3_hist(parameters[1:3], hist, hist_uncertainty),
                             densityMLP3_paleo(parameters[1:3], data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty))+
    (1-parameters[7])*c(densityMLP3_gage(parameters[4:6], data, gage_uncertainty),
                        densityMLP3_hist(parameters[4:6], hist, hist_uncertainty),
                        densityMLP3_paleo(parameters[4:6], data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty))
  return(density)
}