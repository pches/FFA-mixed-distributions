#################################################################################
#
#  -file = "MLE_Fitting_singleandtcev.R"   Code written June 2018
#  - Author: Kenneth Joel Roop-Eckart (kjr30@psu.edu)
#
#  - This code fits the LN2, LP3, GEV, and TCEV models to the gage, historical,
#    and paleo data using maximum likelihoods implemented with DEoptim and 
#    calculates corresponding AIC and BIC values.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# To use this function, simply source this file:
#   source("MLE_Fitting_singleandtcev.R")
#
###################################################################################
# load relevant libraries
library(DEoptim) # for DEoptim function
library(fExtremes) # for GEV functions
library(zipfR) # for the Igamma function 
# set the seed
set.seed(1)
####################################### LN2 MLE FIT #########################################

MLE_LN2_optim<-DEoptim(LN2_gage_hist_paleo_MLE,upper = c(10000, 10000), lower = c(1, 0), 
                       data = data,gage_uncertainty=gage_uncertainty, hist = hist, 
                       hist_record = hist_record, hist_uncertainty=hist_uncertainty,
                       paleo_flow_upper = paleo_flow_upper, paleo_flow_lower = paleo_flow_lower, 
                       paleo_age = paleo_age,paleo_age_uncertainty=paleo_age_uncertainty, 
                       paleo_flow_uncertainty=paleo_flow_uncertainty, control = DEoptim.control(itermax = 500))

###################################### LP3 MLE FIT #########################################

MLE_LP3_optim<-DEoptim(LP3_gage_hist_paleo_MLE,upper = c(100, 100,10), lower = c(0, 0, -4), 
                       data = data,gage_uncertainty=gage_uncertainty, hist = hist, 
                       hist_record = hist_record, hist_uncertainty=hist_uncertainty,
                       paleo_flow_upper = paleo_flow_upper, paleo_flow_lower = paleo_flow_lower, 
                       paleo_age = paleo_age,paleo_age_uncertainty=paleo_age_uncertainty, 
                       paleo_flow_uncertainty=paleo_flow_uncertainty, control = DEoptim.control(itermax = 1000))

##################################### GEV MLE FIT ###########################################

MLE_GEV_optim<-DEoptim(GEV_gage_hist_paleo_MLE,upper = c(3, 1000, 1000), lower = c(-3, 1, 1), 
                       data = data,gage_uncertainty=gage_uncertainty, hist = hist, 
                       hist_record = hist_record, hist_uncertainty=hist_uncertainty,
                       paleo_flow_upper = paleo_flow_upper, paleo_flow_lower = paleo_flow_lower, 
                       paleo_age = paleo_age,paleo_age_uncertainty=paleo_age_uncertainty, 
                       paleo_flow_uncertainty=paleo_flow_uncertainty, control = DEoptim.control(itermax = 500))

####################################### TCEV MLE FIT #############################################

MLE_TCEV_optim<-DEoptim(TCEV_gage_hist_paleo_MLE,upper = c(5000, 5000, 5000, 5000), lower = c(0, 0, 0, 0), 
                        data = data,gage_uncertainty=gage_uncertainty, hist = hist, 
                        hist_record = hist_record, hist_uncertainty=hist_uncertainty,
                        paleo_flow_upper = paleo_flow_upper, paleo_flow_lower = paleo_flow_lower, 
                        paleo_age = paleo_age,paleo_age_uncertainty=paleo_age_uncertainty, 
                        paleo_flow_uncertainty=paleo_flow_uncertainty, control = DEoptim.control(itermax = 2000, trace = 20))

###################### Akaike Information Criteria #############################
LN2_AIC <- 2*MLE_LN2_optim$optim$bestval+2*length(MLE_LN2_optim$optim$bestmem) # LN2 AIC
LP3_AIC <- 2*MLE_LP3_optim$optim$bestval+2*length(MLE_LP3_optim$optim$bestmem) # LP3 AIC
GEV_AIC <- 2*MLE_GEV_optim$optim$bestval+2*length(MLE_GEV_optim$optim$bestmem) # GEV AIC
TCEV_AIC <- 2*MLE_TCEV_optim$optim$bestval+2*length(MLE_TCEV_optim$optim$bestmem) # TCEV AIC


####################### Bayesian Information Criteria ###########################
LN2_BIC <- 2*MLE_LN2_optim$optim$bestval+length(MLE_LN2_optim$optim$bestmem)*log(length(c(data, hist, paleo_flow_lower))) # LN2 BIC
LP3_BIC <- 2*MLE_LP3_optim$optim$bestval+length(MLE_LP3_optim$optim$bestmem)*log(length(c(data, hist, paleo_flow_lower))) # LP3 BIC
GEV_BIC <- 2*MLE_GEV_optim$optim$bestval+length(MLE_GEV_optim$optim$bestmem)*log(length(c(data, hist, paleo_flow_lower))) # GEV BIC
TCEV_BIC <- 2*MLE_TCEV_optim$optim$bestval+length(MLE_TCEV_optim$optim$bestmem)*log(length(c(data, hist, paleo_flow_lower))) # TCEV BIC


