#################################################################################
#
#  -file = "LoadingData_Uncertainty.R"   Code written June 2018
#  - Author: Kenneth Joel Roop-Eckart (kjr30@psu.edu)
#
#  - This code loads the annual peak flow gage data and contains the explicitly
#    stated historical flood and paleoflood values. It also calculates the 
#    uncertianty based on general gage uncertainty estimates 
#    (O'Connell et al., 2002).
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# To use this function, simply source this file:
#   source("LoadingData_Uncertainty.R")
#
###################################################################################
# load relevant libararies
library(triangle) # for the dtriangle function
library(zipfR) # for the Igamma function 
################################### Load Data ######################################

################### Gage data
download.file('https://nwis.waterdata.usgs.gov/nwis/peak?site_no=07099500&agency_cd=USGS&format=rdb', method = 'wininet', 'Data/pueblo_yrpkflow.txt')
hist.flow1 = read.table('Data/pueblo_yrpkflow.txt', header = FALSE, sep = '\t', skip = 66)
# transform data into a dataframe
puebloflow_emperial <- data.frame(hist.flow1$V5)
# convert to cubic meters per second
puebloflow <- puebloflow_emperial$hist.flow1.V5*rep(0.0283168, by = length(puebloflow_emperial))
data = puebloflow

gage_error <- vector(mode = 'numeric', length = length(data))
gage_error[which(data<283)]<-0.1
gage_error[which(data>=283)]<-0.25
gage_uncertainty <- (gage_error*data/1.96)

#################### Historical data
hist_record <- 36 # 36 years (1859-1894) of record prior to the gage (England, 2010)
# alternative interpretation used in (England, 2010) has gage record at 64 yrs (1859-1921)
flood1864 <- 650 # estimated graphically from England et al. (2010)
flood1893 <- 1000 # estimated graphically from England et al. (2010)
flood1894 <- 39100*0.0283168 # (Follansbee and Jones 1922) converted from cubic feet to cubic meters
flood1921 <- 103000*0.0283168 # (Follansbee and Jones 1922) converted from cubic feet to cubic meters
# the 1921 flood was ultimately considered a gage data point for purposes of data processing and fitting
hist <- c(flood1864,flood1893,flood1894)
hist_error <- 0.25 # assume 25% error on all historical flood estimates (O'Connell et al., 2002)
hist_prob <- 0.975 # assume the 25% corresponds to the 95% confidence interval (O'Connell et al., 2002)
hist_uncertainty <- matrix(data = NA, nrow = 1, ncol = length(hist))
uncert<-function(param, flood, error, prob){
  q<-qnorm(p=prob, mean = flood, sd = param)
  return(abs(q-(error*flood+flood)))
}
hist_uncertainty <- (hist_error*hist/1.96)

################### Paleo Data
paleo_age <- (730+840)/2 # preferred age (England et al., 2010) was...
paleo_age_upper <- 840 # (England et al., 2010)
paleo_age_lower <- 730 # (England et al., 2010)

paleo_flow_lower <- 3680 # lower estimate (England et al., 2010)
paleo_flow_upper <- 4530 # upper estimate (England et al., 2010)
paleo <- 4250 # preferred estimate (England et al., 2010) is 4250

# create 13 possible paleoflood magnitudes and assign them probabilities
vals <- seq(from=paleo_flow_lower,to=paleo_flow_upper,by=(paleo_flow_upper-paleo_flow_lower)/12)
probs <- dtriangle(vals,paleo_flow_lower,paleo_flow_upper,paleo)+mean(dtriangle(seq(from=paleo_flow_lower,to=paleo_flow_upper,by=(paleo_flow_upper-paleo_flow_lower)/12),paleo_flow_lower,paleo_flow_upper,paleo))

# create uncertainty matrix
paleo_flow_uncertainty <- matrix(data = NA,nrow = length(vals),ncol = 2*length(paleo_age))
paleo_flow_uncertainty[,1] <- vals
paleo_flow_uncertainty[,2] <- probs/sum(probs)

# create 13 possible paleoflood ages and assign them probabilities
vals <- seq(from=paleo_age_lower,to=paleo_age_upper,by=(paleo_age_upper-paleo_age_lower)/12)
probs<- dtriangle(vals,paleo_age_lower,paleo_age_upper,mean(c(paleo_age_lower,paleo_age_upper)))+mean(dtriangle(vals,paleo_age_lower,paleo_age_upper,mean(c(paleo_age_lower,paleo_age_upper))))

paleo_age_uncertainty <- matrix(data = NA,nrow = length(vals),ncol = 2*length(paleo_age))
paleo_age_uncertainty[,1] <- vals
paleo_age_uncertainty[,2] <- probs/sum(probs)

