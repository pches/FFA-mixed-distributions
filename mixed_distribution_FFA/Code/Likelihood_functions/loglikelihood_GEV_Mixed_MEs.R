



# created by Joel Roop-Eckart with significant input from Ben Lee
# papers used: (Stedinger et al., 1986) and (O'Connell, 2002)
# Stedinger, Jery R., and Timothy A. Cohn. 1986. "Flood Frequency Analysis With Historical and Paleoflood Information." Water Resources Research 22 (5): 785-93. doi:10.1029/WR022i005p00785.
# O'Connell, Daniel R. H. 2002. "Bayesian Flood Frequency Analysis with Paleohydrologic Bound Data." Water Resources Research 38 (5). doi:10.1029/2000WR000028.

# dgev1 and pgev1 are edited from the dgev and pgev1 fExtremes functions to prevent the creation of a data frame every iteration.
# Wuertz, Diethelm, Tobias Setz, and Yohan Chalabi. 2013. fExtremes: Rmetrics - Extreme Financial Market Data. R package version  3010.81. https://CRAN.R-project.org/package=fExtremes
dgev1 <- function(q, xi, beta, mu, log = FALSE){
  shape <- xi
  scale <- beta
  location <- mu
  x <- q
  
  stopifnot(min(scale) > 0)
  stopifnot(length(shape) == 1)
  x = (x - location)/scale
  if (shape == 0) {
    d = log(1/scale) - x - exp(-x)
  }
  else {
    nn = length(x)
    xx = 1 + shape * x
    xxpos = xx[xx > 0 | is.na(xx)]
    scale = rep(scale, length.out = nn)[xx > 0 | is.na(xx)]
    d = numeric(nn)
    d[xx > 0 | is.na(xx)] = log(1/scale) - xxpos^(-1/shape) - 
      (1/shape + 1) * log(xxpos)
    d[xx <= 0 & !is.na(xx)] = -Inf
  }
  if (!log) {
    d = exp(d)
  }
  #  attr(d, "control") = data.frame(location = location[1], scale = scale[1], 
  #                                  shape = shape[1], log = log, row.names = "")
  return(d)
}

pgev1 <- function(q, xi, beta, mu, lower.tail = TRUE) {
  shape <- xi
  scale <- beta
  location <- mu
  
  stopifnot(min(scale) > 0)
  stopifnot(length(shape) == 1)
  q = (q - location)/scale
  if (shape == 0) {
    p = exp(-exp(-q))
  }
  else {
    p = exp(-pmax(1 + shape * q, 0)^(-1/shape))
  }
  if (!lower.tail) {
    p = 1 - p
  }
  #  attr(p, "control") = data.frame(location = location[1], scale = scale[1], 
  #                                  shape = shape[1], lower.tail = lower.tail, row.names = "")
  return(p)
}
################################### Gage Component ###################################
densityMGEV_gage<-function(parameters,data,gage_uncertainty){
  like <- matrix(data = 0, ncol = length(data), nrow = 1)
  I = length(data)
  for(i in 1:I){
    Xij <- qnorm(seq(from = 0.1, to = 0.9, by = 0.1), mean = data[i], sd = gage_uncertainty[i])
    Fij <- dnorm(Xij, mean = data[i], sd = gage_uncertainty[i])
    Fij <- Fij/sum(Fij)
    like[1,i] <- (
      sum(
        Fij*dgev1(Xij, xi = parameters[1],  beta = parameters[2], mu = parameters[3])
      )
    )
  }
  if(is.na(min(like))==TRUE)(like[which(is.na(like)==TRUE)]<-0)
  return(like)
}
###################################### Historical Component ####################################
densityMGEV_hist<-function(parameters,hist,hist_uncertainty){
  like <- matrix(data = 0, ncol = length(hist), nrow = 1)
  I = length(hist)
  for(i in 1:I){
    Xij <- qnorm(seq(from = 0.1, to = 0.9, by = 0.1), mean = hist[i], sd = hist_uncertainty[i])
    Fij <- dnorm(Xij, mean = hist[i], sd = hist_uncertainty[i])
    Fij <- Fij/sum(Fij)
    like[1,i] <- (
      sum(
        Fij*dgev1(Xij, xi = parameters[1], beta = parameters[2], mu = parameters[3])
      )
    )
  }
  out<-(like*(dbinom(length(hist), size = hist_record, prob = (1-pgev1(min(hist), xi = parameters[1], beta = parameters[2], mu = parameters[3]))))/
          (1-pgev1(min(hist), xi = parameters[1], beta = parameters[2], mu = parameters[3])))
  if(is.na(min(out))==TRUE)(out[which(is.na(out)==TRUE)]<-0)
  return(out)
}
########################################### Paleo Component ################################
densityMGEV_paleo <- function(parameters, data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty){
  #  I <- length(paleo_age)
  #  paleo_age <- paleo_age
  #  Uik <- paleo_flow_upper
  #  Xmin <- paleo_flow_lower
  #  like <- vector(mode = 'numeric', length = length(paleo_age))
  #  for(i in 1:I){
  #    # time component of the likelihood function
  #    Nij <- paleo_age_uncertainty[,1] # dstribution of ages
  #    FUnij <- paleo_age_uncertainty[,2]
  #    # flood size component of the likelihood function
  #    variable <- paleo_flow_uncertainty[,1] # distribution of flows based on a normal distribution around a best estimate
  #    Uik <- (variable)
  #    FUdik <- paleo_flow_uncertainty[,2]
  #    intvect <- pgev1(Uik, xi = parameters[1], beta = parameters[2], mu = parameters[3])-pgev1(Xmin, xi = parameters[1], beta = parameters[2], mu = parameters[3]) # the probability space between the upper and lower estimates
  #    
  #    # NOTE: was sum(FUnij*Nij)*log(sum(FUdik*intvect)), used exponential to remove log
  #    like[i] <- sum(FUnij*Nij)*#log(
  #      sum(FUdik*intvect # multiply the discrete probability distribution by the probability space between the paleo upper and lower bounds
  #      )#)
  #  }
  
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
  
  like <- dgev1(Uik, xi = parameters[1], beta = parameters[2], mu = parameters[3])
  
  age_like <- vector(mode = 'numeric', length = I)
  for(i in 1:I){
    age_like[i] <- dbinom(length(paleo_age[i]), size = round(paleo_age[i]), prob = (1-pgev1(paleo_flow_lower, xi = parameters[1], beta = parameters[2], mu = parameters[3])))/
      (1-pgev1(paleo_flow_lower, xi = parameters[1], beta = parameters[2], mu = parameters[3]))
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
  
  # thesis way
  #    for(i in 1:I){
  #    FUdik <- paleo_flow_uncertainty[,2]
  #    variable <- paleo_flow_uncertainty[,1]
  #    Uik <- (variable)
  #    if(parameters[1] < 0)(return(rep(0,length(paleo))))
  #    else if(parameters[2] < 0)(return(rep(0,length(paleo))))
  #    else{
  #      like[i] <- sum(
  #        FUdik*dgev1(Uik, xi = parameters[1], beta = parameters[2], mu = parameters[3])
  #      )
  #    }
  #  }
  
  #  loglik_MEs <- sum(like)
  #  out<-(like*((dbinom(length(paleo_age), size = paleo_age, prob = 1-pgev1(paleo_flow_lower, xi = parameters[1], beta = parameters[2], mu = parameters[3])))/
  #      (1-pgev1(paleo_flow_lower, xi = parameters[1], beta = parameters[2], mu = parameters[3]))))
  if(is.na(min(out))==TRUE)(out[which(is.na(out)==TRUE)]<-0)
  return(out)
}
####################################### Complete Likelihood Function ##############################
GEV_gage_hist_paleo_mixed_test <- function(parameters, data, gage_uncertainty, hist, hist_record, hist_uncertainty, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty){
  density <-parameters[7]*c(densityMGEV_gage(parameters[1:3], data, gage_uncertainty),
                            densityMGEV_hist(parameters[1:3], hist, hist_uncertainty),
                            densityMGEV_paleo(parameters[1:3], data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty))+
    (1-parameters[7])*c(densityMGEV_gage(parameters[4:6], data, gage_uncertainty),
                        densityMGEV_hist(parameters[4:6], hist, hist_uncertainty),
                        densityMGEV_paleo(parameters[4:6], data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty))
  loglik <- sum(log(density))
  if(is.na(loglik)==TRUE)(loglik=-Inf)
  return(loglik)
}


# for MLE fitting, returns positive loglikelihood value for input into DEoptim
GEV_gage_hist_paleo_mixed_MLE_test <- function(parameters, data, gage_uncertainty, hist, hist_record, hist_uncertainty, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty){
  density <- parameters[7]*c(densityMGEV_gage(parameters[1:3], data, gage_uncertainty),
                             densityMGEV_hist(parameters[1:3], hist, hist_uncertainty),
                             densityMGEV_paleo(parameters[1:3], data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty))+
    (1-parameters[7])*c(densityMGEV_gage(parameters[4:6], data, gage_uncertainty),
                        densityMGEV_hist(parameters[4:6], hist, hist_uncertainty),
                        densityMGEV_paleo(parameters[4:6], data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty))
  loglik <- sum(log(density))
  if(is.na(loglik)==TRUE)(loglik=-Inf)
  #  loglik<-(log(dnorm(parameters[7], mean = (1-0.245749),sd = 0.047778)))+loglik
  return(-loglik)
}

# function
densityMGEV <- function(parameters, data, gage_uncertainty, hist, hist_record, hist_uncertainty, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty){
  density <- parameters[7]*c(densityMGEV_gage(parameters[1:3], data, gage_uncertainty),
                             densityMGEV_hist(parameters[1:3], hist, hist_uncertainty),
                             densityMGEV_paleo(parameters[1:3], data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty))+
    (1-parameters[7])*c(densityMGEV_gage(parameters[4:6], data, gage_uncertainty),
                        densityMGEV_hist(parameters[4:6], hist, hist_uncertainty),
                        densityMGEV_paleo(parameters[4:6], data, paleo_flow_upper, paleo_flow_lower, paleo_age, paleo_age_uncertainty, paleo_flow_uncertainty))
  return(density)
}