#################################################################################
#
#  -file = "MLE_mixedGEV.R"   Code written June 2018
#  - Author: Kenneth Joel Roop-Eckart (kjr30@psu.edu)
#
#  - This code fits the Mixed GEV to the gage, historical, and paleo data
#    by using four indpendently initiated runs of DEoptim to search for
#    the maximum likelihood parameter values. The parameter estimate with
#    the best likelihood value is accepted as the MLE.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# To use this function, simply source this file:
#   source("MLE_mixedGEV.R")
#
###################################################################################
################################### Run four DEoptim runs from different seeds ###########################
# load relevant libraries
library(DEoptim) # for DEoptim function
library(fExtremes) # for GEV functions
# start
iter <- 2500
set.seed(1)
MGEVrun1<-DEoptim(GEV_gage_hist_paleo_mixed_MLE_test,upper = c(3, 250, 500,3, 250, 500,0.9), lower = c(-4, 0, 0,-4, 0, 0, 0.6), 
                  data = data,gage_uncertainty=gage_uncertainty, hist = hist, 
                  hist_record = hist_record, hist_uncertainty=hist_uncertainty,
                  paleo_flow_upper = paleo_flow_upper, paleo_flow_lower = paleo_flow_lower, 
                  paleo_age = paleo_age,paleo_age_uncertainty=paleo_age_uncertainty, 
                  paleo_flow_uncertainty=paleo_flow_uncertainty, control = DEoptim.control(itermax = iter, NP = 500, F = 0.7, CR = 0.5))
set.seed(2)
MGEVrun2<-DEoptim(GEV_gage_hist_paleo_mixed_MLE_test,upper = c(3, 250, 500,3, 250, 500,0.9), lower = c(-4, 0, 0,-4, 0, 0, 0.6), 
                  data = data,gage_uncertainty=gage_uncertainty, hist = hist, 
                  hist_record = hist_record, hist_uncertainty=hist_uncertainty,
                  paleo_flow_upper = paleo_flow_upper, paleo_flow_lower = paleo_flow_lower, 
                  paleo_age = paleo_age,paleo_age_uncertainty=paleo_age_uncertainty, 
                  paleo_flow_uncertainty=paleo_flow_uncertainty, control = DEoptim.control(itermax = iter, NP = 500, F = 0.7, CR = 0.5))
set.seed(4)
MGEVrun4<-DEoptim(GEV_gage_hist_paleo_mixed_MLE_test,upper = c(3, 250, 500,3, 250, 500,0.9), lower = c(-4, 0, 0,-4, 0, 0, 0.6), 
                  data = data,gage_uncertainty=gage_uncertainty, hist = hist, 
                  hist_record = hist_record, hist_uncertainty=hist_uncertainty,
                  paleo_flow_upper = paleo_flow_upper, paleo_flow_lower = paleo_flow_lower, 
                  paleo_age = paleo_age,paleo_age_uncertainty=paleo_age_uncertainty, 
                  paleo_flow_uncertainty=paleo_flow_uncertainty, control = DEoptim.control(itermax = iter, NP = 500, F = 0.7, CR = 0.5))
set.seed(5)
MGEVrun5<-DEoptim(GEV_gage_hist_paleo_mixed_MLE_test,upper = c(3, 250, 500,3, 250, 500,0.9), lower = c(-4, 0, 0,-4, 0, 0, 0.6), 
                  data = data,gage_uncertainty=gage_uncertainty, hist = hist, 
                  hist_record = hist_record, hist_uncertainty=hist_uncertainty,
                  paleo_flow_upper = paleo_flow_upper, paleo_flow_lower = paleo_flow_lower, 
                  paleo_age = paleo_age,paleo_age_uncertainty=paleo_age_uncertainty, 
                  paleo_flow_uncertainty=paleo_flow_uncertainty, control = DEoptim.control(itermax = iter, NP = 500, F = 0.7, CR = 0.5))

############################## AIC ##################################
MGEVrun1_AIC <- 2*MGEVrun1$optim$bestval+2*length(MGEVrun1$optim$bestmem) # MGEV AIC
MGEVrun2_AIC <- 2*MGEVrun2$optim$bestval+2*length(MGEVrun1$optim$bestmem) # MGEV AIC
MGEVrun4_AIC <- 2*MGEVrun4$optim$bestval+2*length(MGEVrun1$optim$bestmem) # MGEV AIC
MGEVrun5_AIC <- 2*MGEVrun5$optim$bestval+2*length(MGEVrun1$optim$bestmem) # MGEV AIC
############################## BIC ###################################
MGEVrun1_BIC <- 2*MGEVrun1$optim$bestval+length(MGEVrun1$optim$bestmem)*
  log(length(data)+length(hist)+length(paleo)) # MGEV BIC
MGEVrun2_BIC <- 2*MGEVrun2$optim$bestval+length(MGEVrun1$optim$bestmem)*
  log(length(data)+length(hist)+length(paleo)) # MGEV BIC
MGEVrun4_BIC <- 2*MGEVrun4$optim$bestval+length(MGEVrun1$optim$bestmem)*
  log(length(data)+length(hist)+length(paleo)) # MGEV BIC
MGEVrun5_BIC <- 2*MGEVrun5$optim$bestval+length(MGEVrun1$optim$bestmem)*
  log(length(data)+length(hist)+length(paleo)) # MGEV BIC
############################# set run 1 chosen as MLE value due to having the best likelihood
MLE_MGEV_optim <- MGEVrun5
