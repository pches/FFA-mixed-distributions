



# created by Joel Roop-Eckart and Ben Lee
# papers used: (Stedinger et al., 1986) and (O'Connell, 2002)
# Stedinger, Jery R., and Timothy A. Cohn. 1986. "Flood Frequency Analysis With Historical and Paleoflood Information." Water Resources Research 22 (5): 785-93. doi:10.1029/WR022i005p00785.
# O'Connell, Daniel R. H. 2002. "Bayesian Flood Frequency Analysis with Paleohydrologic Bound Data." Water Resources Research 38 (5). doi:10.1029/2000WR000028.

################################### Gage Component ###################################
densityLN2_gage<-function(parameters,data,gage_uncertainty){
  #m<-min((1+(parameters[1]*(data-parameters[3])/parameters[2])))
  #if(m<0.00001)(return(Inf))
  #else if(any(parameters[2]<0.00001))(return(Inf))
  #if(parameters[1] < 0)(return(Inf))
  if(parameters[2] < 0)(return(rep(0,length(data))))
  else{
    like <- vector(mode = 'numeric', length = length(data))
    I = length(data)
    for(i in 1:I){
      Xij <- qnorm(seq(from = 0.025, to = 0.975, by = 0.095), mean = data[i], sd = gage_uncertainty[i])
      Fij <- dnorm(Xij, mean = data[i], sd = gage_uncertainty[i])
      Fij <- Fij/sum(Fij)
      #m<-min((1+(parameters[1]*(data-parameters[3])/parameters[2])))
      #if(m<0.00001)(return(Inf))
      #if(parameters[1] < 0)(return(Inf))
      if(parameters[2] < 0)(return(rep(0,length(data))))
      else{
        like[i] <- sum(
          Fij*dlnorm(Xij, parameters[1], parameters[2])
        )
      }
    }
    
    #    loglik_MEs <-  sum(like)
    #    if(is.nan(loglik_MEs)==TRUE)(return(Inf))
    return(like)
  }
}
###################################### Historical Component ####################################
densityLN2_hist<-function(parameters,hist,hist_uncertainty){
  #m<-min((1+(parameters[1]*(data-parameters[3])/parameters[2])))
  #if(m<0.00001)(return(Inf))
  #else if(any(parameters[2]<0.00001))(return(Inf))
  #if(parameters[1] < 0)(return(Inf))
  if(parameters[2] < 0)(return(rep(0,length(hist))))
  else{
    like <- matrix(data = 0, ncol = length(hist), nrow = 1)
    I = length(hist)
    for(i in 1:I){
      Xij <- qnorm(seq(from = 0.025, to = 0.975, by = 0.095), mean = hist[i], sd = hist_uncertainty[i])
      Fij <- dnorm(Xij, mean = hist[i], sd = hist_uncertainty[i])
      Fij <- Fij/sum(Fij)
      #m<-min((1+(parameters[1]*(Xij-parameters[3])/parameters[2])))
      #if(m<0.00001)(return(Inf))
      #if(parameters[1] < 0)(return(Inf))
      if(parameters[2] < 0)(return(rep(0,length(hist))))
      else{
        like[1,i] <- sum(
          Fij*dlnorm(Xij, parameters[1], parameters[2])
        )
      }
    }
    
    #    loglik_MEs <-  sum(like)
    out <- (like*dbinom(length(hist), size = hist_record, prob = 1-plnorm(min(hist), parameters[1], parameters[2]))/
              (1-plnorm(min(hist), parameters[1], parameters[2])))
  }
  return(out)
}
########################################### Paleo Component ################################
densityLN2_paleo <- function(parameters, data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty){
  #Nij <- paleo_age
  #  I <- length(paleo_age)
  #  paleo_age <- paleo_age
  #  Uik <- paleo_flow_upper
  #  Xmin <- paleo_flow_lower
  #  parameters <- parameters
  #  like <- vector(mode = 'numeric', length = length(paleo_age))
  #  for(i in 1:I){
  #    # time component of the likelihood function
  #    Nij <- paleo_age_uncertainty[,1] # distribution of ages
  #    FUnij <- paleo_age_uncertainty[,2]
  #    # flood size component of the likelihood function
  #    variable <- paleo_flow_uncertainty[,1] # distribution of flows based on a normal distribution around a best estimate
  #    Uik <- (variable)
  #    FUdik <- paleo_flow_uncertainty[,2]
  #    intvect <- plnorm(Uik, parameters[1], parameters[2])-plnorm(Xmin, parameters[1], parameters[2]) # the probability space between the upper and lower estimates
  #    
  #    like[i] <- sum(FUnij*Nij)*#log(
  #      sum(FUdik*intvect # multiply the discrete probability distribution by the probability space between the paleo upper and lower bounds
  #      )#)
  #  }
  #  loglik_MEs <- sum(like)
  
  I <- length(paleo_age)
  paleo_age <- paleo_age
  variable <- paleo_flow_uncertainty[,1]
  Uik <- (variable)
  Xmin <- paleo_flow_lower
  parameters <- parameters
  like <- vector(mode = 'numeric', length = length(paleo_age))
  
  Uik <- paleo_flow_uncertainty[,1]
  paleo_age <- paleo_age_uncertainty[,1]
  I <- length(paleo_age)
  probs <- (paleo_flow_uncertainty[,2]*paleo_age_uncertainty[,2])/(sum(paleo_flow_uncertainty[,2]*paleo_age_uncertainty[,2]))
  
  like <- dlnorm(Uik, parameters[1], parameters[2])
  
  age_like <- vector(mode = 'numeric', length = I)
  for(i in 1:I){
    age_like[i] <- dbinom(length(paleo_age[i]), size = round(paleo_age[i]), prob = (1-plnorm(paleo_flow_lower, parameters[1], parameters[2])))/
      (1-plnorm(paleo_flow_lower, parameters[1], parameters[2]))
  }
  if(is.na(min(age_like))==TRUE)(which(is.nan(age_like))==TRUE)
  
 # out <- sum(like*age_like*probs)
  ############## Addition
  out <- vector(mode = 'numeric', length = length(age_like))
  for (i in 1:length(age_like)) {
    out[i] <- sum(like*probs*age_like[i])
  }
  out <- sum(out)
  ############## End addition
  
  
  ############ testing paleo-hist combo
  #  for(i in 1:I){
  #    FUdik <- paleo_flow_uncertainty[,2]
  #    variable <- paleo_flow_uncertainty[,1]
  #    Uik <- (variable)
  #    if(parameters[1] < 0)(return(rep(0,length(paleo))))
  #    else if(parameters[2] < 0)(return(rep(0,length(paleo))))
  #    else{
  #      like[i] <- sum(
  #        FUdik*dlnorm(Uik, parameters[1], parameters[2])
  #      )
  #    }
  #  }
  #  
  #  
  #  out <-(like*dbinom(length(paleo_age), size = paleo_age, prob = 1-plnorm(paleo_flow_lower, parameters[1], parameters[2]))/
  #          (1-plnorm(paleo_flow_lower, parameters[1], parameters[2])))
  return(out)
}
####################################### Complete Likelihood Function ##############################
# gage historical and paleo likelihood function
LN2_gage_hist_paleo <- function(parameters, data, gage_uncertainty, hist, hist_record, hist_uncertainty, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty){
  density <- c(densityLN2_gage(parameters, data, gage_uncertainty),
               densityLN2_hist(parameters, hist, hist_uncertainty),
               densityLN2_paleo(parameters, data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty))
  loglik <- sum(log(density))
  if(is.na(loglik)==TRUE)(loglik=-Inf)
  return(loglik)
}
# for MLE
LN2_gage_hist_paleo_MLE <- function(parameters, data, gage_uncertainty, hist, hist_record, hist_uncertainty, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty){
  density <- c(densityLN2_gage(parameters, data, gage_uncertainty),
               densityLN2_hist(parameters, hist, hist_uncertainty),
               densityLN2_paleo(parameters, data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty))
  loglik <- sum(log(density))
  if(is.na(loglik)==TRUE)(loglik=-Inf)
  return(-loglik)
}
######################################### density functions #########################################
# gage historical and paleo likelihood function
densityLN2 <- function(parameters, data, gage_uncertainty, hist, hist_record, hist_uncertainty, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty){
  density <- c(densityLN2_gage(parameters, data, gage_uncertainty),
               densityLN2_hist(parameters, hist, hist_uncertainty),
               densityLN2_paleo(parameters, data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty))
  return(density)
}