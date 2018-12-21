#################################################################################
#
#  -file = "Master_script_thesis.R"   Code written June 2018
#  - Author: Kenneth Joel Roop-Eckart (kjr30@psu.edu)
#
#  - This code is the master script that runs each script for the analysis in
#    the proper order.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# To use this code, simply source this file:
#   source("Master_script_thesis.R")
#
###################################################################################

setwd("C:/Users/kjr30.SCRiM-Win-DT01/Desktop/Thesis_code_forBOX") # set working directory to the "Thesis" folder
# IMPORTANT NOTE: the working directory must be changed in three places.
# Line 20 of this code, line 70 of this code, and line 27 of the hydrograph_overtopping_script.R 

######################## Install packages
install.packages("DEoptim")
install.packages("fExtremes")
install.packages("triangle")
install.packages("RColorBrewer")
install.packages("zipfR")
######################## Load Packages
library(DEoptim) # for DEoptim function
library(fExtremes) # for GEV functions
library(triangle) # for the dtriangle function
library(RColorBrewer) # for colorblind friendly colors
library(zipfR) # for the Igamma function 

####################### Load Pueblo Data ################################
# script loads the data for Pueblo, calculates uncertainties

# NOTE: when using Mac, change "method" from "wininet" to "curl on line 5
source('Code/LoadingData_Uncertainty.R')

#################### MLE: LN2, LP3, GEV, TCEV #################################
# load likelihood functions
source("Code/Likelihood_functions/loglikelihood_LN2_MEs.R") # likelihood function for the LN2
source("Code/Likelihood_functions/loglikelihood_LP3_MEs.R") # likelihood function for the LP3
source("Code/Likelihood_functions/loglikelihood_GEV_MEs.R") # likelihood function for the GEV
source('Code/Likelihood_functions/loglikelihood_TCEV_MEs.R') # likelihood function for the TCEV
# script fits the single distirubitions and the TCEV by maximum likelihood
source('Code/MLE_Fitting_singleandtcev.R')

######################## MLE: MLP3, MGEV #################################
# load likelihood functions
source('Code/Likelihood_functions/loglikelihood_mixedLP3_MEs_serioustesting.R')
source('Code/Likelihood_functions/loglikelihood_mixedGEV_MEs_serioustesting.R')
# script fits the Mixed Log Pearson III and the Mixed Generalized Extreme Value distributions by maximum likelihood
source('Code/MLE_mixedLP3.R')
source('Code/MLE_mixedGEV.R')

#################### scale hydrographs, find overtopping magnitudes ################################
# script scales the three hydrographs to find overtopping peak flows for the Pueblo Dam
source('Code/TREX_hydrographs_analysis.R') # load the TREX hydrographs
source('Code/hydro_scale.R') # load the hydrograph scaling function
source('Code/hydrograph_overtopping_script.R') # determine peak flows for overtopping for each hydrograph

# NOTE: hydrograph_overtopping_script.R changes the working directory to the "Data" folder
# so it must be changed back to "Thesis" for figure plotting
############################# create figures ######################################################
# script creates output figures from the analysis
setwd("C:/Users/kjr30.SCRiM-Win-DT01/Desktop/Thesis_code_forBOX") # set working directory to the "Thesis" folder
source('Code/PlottingPositions.R') # calculate plotting position return periods of peak flow data
source('Code/Thesis Figure Script Final.R')

print("FINISHED")
