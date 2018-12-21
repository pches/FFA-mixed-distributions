#################################################################################
#
#  -file = "MLE_mixedLP3.R"   Code written June 2018
#  - Author: Kenneth Joel Roop-Eckart (kjr30@psu.edu)
#
#  - This code fits the Mixed LP3 to the gage, historical, and paleo data
#    by using two indpendently initiated runs of DEoptim to search for
#    the maximum likelihood parameter values. The independent runs converge 
#    to a single MLE.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# To use this function, simply source this file:
#   source("MLE_mixedLP3.R")
#
###################################################################################
################################### Run two DEoptim runs from different seeds ###########################
# load relevant libraries
library(DEoptim) # for DEoptim function
library(zipfR) # for the Igamma function 
# start
iter <- 2000
set.seed(1)
MLP3run1<-DEoptim(LP3_gage_hist_paleo_mixed_MLE_test,upper = c(10, 100, 5, 10, 100, 5, 0.9), lower = c(0, 0, 0, 0, 0, 0, 0.5), 
                  data = data,gage_uncertainty=gage_uncertainty, hist = hist, 
                  hist_record = hist_record, hist_uncertainty=hist_uncertainty,
                  paleo_flow_upper = paleo_flow_upper, paleo_flow_lower = paleo_flow_lower, 
                  paleo_age = paleo_age,paleo_age_uncertainty=paleo_age_uncertainty, 
                  paleo_flow_uncertainty=paleo_flow_uncertainty, control = DEoptim.control(itermax = iter, NP = 500, F = 0.7, CR = 0.5))

set.seed(2)
MLP3run2<-DEoptim(LP3_gage_hist_paleo_mixed_MLE_test,upper = c(10, 100, 5, 10, 100, 5, 0.9), lower = c(0, 0, 0, 0, 0, 0, 0.5), 
                  data = data,gage_uncertainty=gage_uncertainty, hist = hist, 
                  hist_record = hist_record, hist_uncertainty=hist_uncertainty,
                  paleo_flow_upper = paleo_flow_upper, paleo_flow_lower = paleo_flow_lower, 
                  paleo_age = paleo_age,paleo_age_uncertainty=paleo_age_uncertainty, 
                  paleo_flow_uncertainty=paleo_flow_uncertainty, control = DEoptim.control(itermax = iter, NP = 500, F = 0.7, CR = 0.5))

############################## AIC ##################################
MLP3run1_AIC <- 2*MLP3run1$optim$bestval+2*length(MLP3run1$optim$bestmem) # MLP3 AIC
MLP3run2_AIC <- 2*MLP3run2$optim$bestval+2*length(MLP3run2$optim$bestmem) # MLP3 AIC
############################## BIC ###################################
MLP3run1_BIC <- 2*MLP3run1$optim$bestval+length(MLP3run1$optim$bestmem)*log(length(data)+length(hist)+length(paleo)) # MLP3 BIC
MLP3run2_BIC <- 2*MLP3run2$optim$bestval+length(MLP3run2$optim$bestmem)*log(length(data)+length(hist)+length(paleo)) # MLP3 BIC

############################# set run 1 (arbitrarily chosen) to MLE value
MLE_MLP3_optim <- MLP3run1
