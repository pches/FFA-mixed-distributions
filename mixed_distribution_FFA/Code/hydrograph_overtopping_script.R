#################################################################################
#
#  -file = "hydrograph_overtopping_script.R"   Code written June 2018
#  - Author: Kenneth Joel Roop-Eckart (kjr30@psu.edu)
#
#  - This code calculates peak flows for each hydrograph that cause a number
#    of reservoir elevations, including overtopping (4931.27 ft) using the
#    FLROUT program and the hydrograph scaling function ("hydro_scale.R").
#    It also calculates a timeseries of reservoir elevations for each
#    hydrograph scaled to the 1921 peak flow.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# To use this code, simply source this file:
#   source("hydrograph_overtopping_script.R")
#
###################################################################################
setwd("C:/Users/kjr30.SCRiM-Win-DT01/Desktop/Thesis_code_forBOX/Data") # set working directory to the "Data" folder in the "mixed_distribution_FFA" folder
set.seed(1)
############################# Functions to wrap Fortran FLROUT #######################################
# background variables
CV_elev <- 4931.27
overtop_elev <- 4931.27
params = c(10000)
hydrograph <- read.csv(file = 'PMFnosnowrout.csv', sep = ',', header = FALSE, skip = 10)

# run fortran command function for Windows
cmd_exec = function(command){ # string must be quoted
  shell(command, intern = TRUE, mustWork = TRUE)
}

# run fortran command function for Mac/Linux (Remove # from lines 37-39 and 42, comment out 32-34 and 41 to run on Max/Linux)
#cmd_exec = function(command){ # string must be quoted
#  system(command, intern = TRUE)
#}

flrout <- paste("flrout-kjr-edited.exe") # for Windows
#flrout <- paste("flrout-kjr-edited.out") # for Mac/Linux

# function for overtopping without snowmelt included
overtopping <- function(params, hydrograph, overtop_elev){
  hydrograph_nosnow(params, hydrograph)
  # call edited, compiled FLROUT model, input hydrograph data file
  cmd_exec("echo data.in | flrout-kjr-edited.exe") # get(flrout) == flrout-kjr-edited.exe, can change out if script does not work
  # load the output data file
  output <- read.fwf(file = "data.out", 
                     widths = c(10,10,11,11,11,11,11,9), skip = 53, n = 500,
                     col.names = c('TIME', 'INFLOW', 'OUTFLOW', 'STORAGE', 'ELEVATION', 
                                   'WATERWAY 1', 'WATERWAY 2', 'RATING CURVE'))
  return(abs(overtop_elev - max(output[,5])))
}
####################### Load Data for peakflow/elevation calculations (no snowmelt) ################
a <- PMFhydrograph_nosnow <- read.csv(file = 'PMFnosnowrout.csv', sep = ',', header = FALSE, skip = 10)
b <- floodofrecord <- read.csv(file = '1921flood.csv', sep = ',', header = FALSE, skip = 1)
c <- PMFsnowonly <- read.csv(file = 'pmfsnowonly.csv', sep = ',', header = FALSE, skip = 10)
d <- TREX_list # from TREX_hydrographs_analysis.R
################################ Dam Elevation 4931.27 #####################################
# based on the simplicity of the problem, we use the following solution method
# it has numerous pitfalls and is not recommended for any other uses, check the third column
# if it has converged on the correct value, the third column will equal 0.00

overtop_elev <- 4931.27
peakflows_4931.27 <- data.frame(0)
peakflows_4931.27[1,1] <- "PMF_nosnow_4931.27"
peakflows_4931.27[2,1] <- "floodofrecord_4931.27"
peakflows_4931.27[3,1] <- "TREX_nosnow_4931.27"
hydrograph <- a
elev_old = 10
elev_new = 10
peakflow_old = 10000
peakflow = 10000
iterations = 1000
for(i in 1:iterations){
  peakflow = peakflow_old + sample(rnorm(100,mean = 0, sd = (500*elev_old)), 1)
  elev_new <- overtopping(peakflow, hydrograph, overtop_elev)
  if(elev_new < elev_old){peakflow_old <- peakflow 
  elev_old <- elev_new}
  print(c((100/iterations*i),peakflow_old,elev_old))
  if(elev_old == 0){break}
}
peakflows_4931.27[1,2] <- round(peakflow_old)

hydrograph <- b
elev_old = 10
elev_new = 10
peakflow_old = 10000
peakflow = 10000
iterations = 1000
for(i in 1:iterations){
  peakflow = peakflow_old + sample(rnorm(100,mean = 0, sd = (500*elev_old)), 1)
  elev_new <- overtopping(peakflow, hydrograph, overtop_elev)
  if(elev_new < elev_old){peakflow_old <- peakflow 
  elev_old <- elev_new}
  print(c((100/iterations*i),peakflow_old,elev_old))
  if(elev_old == 0){break}
}
peakflows_4931.27[2,2] <- round(peakflow_old)

#hydrograph <- d
elev_old = 10
elev_new = 10
peakflow_old = 10000
peakflow = 10000
iterations = 1000
for(i in 1:iterations){
  peakflow = peakflow_old + sample(rnorm(100,mean = 0, sd = (300*elev_old)), 1)
  hydrograph[,2] <- abs(as.numeric(unlist(d[which(abs(TREX_list_max[1:15]-peakflow)==min(abs(TREX_list_max[1:15]-peakflow)))])))
  elev_new <- overtopping(peakflow, hydrograph, overtop_elev)
  if(elev_new < elev_old){peakflow_old <- peakflow 
  elev_old <- elev_new}
  print(c((100/iterations*i),peakflow_old,elev_old))
  if(elev_old == 0){break}
}
peakflows_4931.27[3,2] <- round(peakflow_old)

# all data reported in m^3/s
# scaled hydrographs, none include snowmelt
# dam crest with freeboard: 4931.27
# overtopping flow, PMF, nosnow: 17261.0
print(peakflows_4931.27[1,])
# overtopping flow, 1921 flood: 22293.4
print(peakflows_4931.27[2,])
# overtopping flow, TREX hydrograph: 14221.02
print(peakflows_4931.27[3,])

################################### Dam Elevation 4925 ##################################
overtop_elev <- 4925 # Target reservoir elevation, not actual overtopping elevation
peakflows_4925 <- data.frame(0)
peakflows_4925[1,1] <- "PMF_nosnow_4925"
peakflows_4925[2,1] <- "floodofrecord_4925"
peakflows_4925[3,1] <- "TREX_nosnow_4925"
hydrograph <- a
elev_old = 10
elev_new = 10
peakflow_old = 15000
peakflow = 10000
iterations = 1000
for(i in 1:iterations){
  peakflow = peakflow_old + sample(rnorm(100,mean = 0, sd = (500*elev_old)), 1)
  elev_new <- overtopping(peakflow, hydrograph, overtop_elev)
  if(elev_new < elev_old){peakflow_old <- peakflow 
  elev_old <- elev_new}
  print(c((100/iterations*i),peakflow_old,elev_old))
  if(elev_old == 0){break}
}
peakflows_4925[1,2] <- round(peakflow_old)

hydrograph <- b
elev_old = 10
elev_new = 10
peakflow_old = 15000
peakflow = 10000
iterations = 1000
for(i in 1:iterations){
  peakflow = peakflow_old + sample(rnorm(100,mean = 0, sd = (500*elev_old)), 1)
  elev_new <- overtopping(peakflow, hydrograph, overtop_elev)
  if(elev_new < elev_old){peakflow_old <- peakflow 
  elev_old <- elev_new}
  print(c((100/iterations*i),peakflow_old,elev_old))
  if(elev_old == 0){break}
}
peakflows_4925[2,2] <- round(peakflow_old)

#hydrograph <- d
elev_old = 10
elev_new = 10
peakflow_old = 15000
peakflow = 10000
iterations = 1000
for(i in 1:iterations){
  peakflow = peakflow_old + sample(rnorm(100,mean = 0, sd = (300*elev_old)), 1)
  hydrograph[,2] <- abs(as.numeric(unlist(d[which(abs(TREX_list_max[1:15]-peakflow)==min(abs(TREX_list_max[1:15]-peakflow)))])))  
  elev_new <- overtopping(peakflow, hydrograph, overtop_elev)
  if(elev_new < elev_old){peakflow_old <- peakflow 
  elev_old <- elev_new}
  print(c((100/iterations*i),peakflow_old,elev_old))
  if(elev_old == 0){break}
}
peakflows_4925[3,2] <- round(peakflow_old)

# all data reported in m^3/s
# scaled hydrographs, none include snowmelt

################################### Dam Elevation 4919 ##################################
overtop_elev <- 4919
peakflows_4919 <- data.frame(0)
peakflows_4919[1,1] <- "PMF_nosnow_4919"
peakflows_4919[2,1] <- "floodofrecord_4919"
peakflows_4919[3,1] <- "TREX_nosnow_4919"
hydrograph <- a
elev_old = 10
elev_new = 10
peakflow_old = 10000
peakflow = 10000
iterations = 1000
for(i in 1:iterations){
  peakflow = peakflow_old + sample(rnorm(100,mean = 0, sd = (500*elev_old)), 1)
  elev_new <- overtopping(peakflow, hydrograph, overtop_elev)
  if(elev_new < elev_old){peakflow_old <- peakflow 
  elev_old <- elev_new}
  print(c((100/iterations*i),peakflow_old,elev_old))
  if(elev_old == 0){break}
}
peakflows_4919[1,2] <- round(peakflow_old)

hydrograph <- b
elev_old = 10
elev_new = 10
peakflow_old = 10000
peakflow = 10000
iterations = 1000
for(i in 1:iterations){
  peakflow = peakflow_old + sample(rnorm(100,mean = 0, sd = (500*elev_old)), 1)
  elev_new <- overtopping(peakflow, hydrograph, overtop_elev)
  if(elev_new < elev_old){peakflow_old <- peakflow 
  elev_old <- elev_new}
  print(c((100/iterations*i),peakflow_old,elev_old))
  if(elev_old == 0){break}
}
peakflows_4919[2,2] <- round(peakflow_old)

#hydrograph <- d
elev_old = 10
elev_new = 10
peakflow_old = 10000
peakflow = 10000
iterations = 1000
for(i in 1:iterations){
  peakflow = peakflow_old + sample(rnorm(100,mean = 0, sd = (300*elev_old)), 1)  
  hydrograph[,2] <- abs(as.numeric(unlist(d[which(abs(TREX_list_max[1:15]-peakflow)==min(abs(TREX_list_max[1:15]-peakflow)))])))  
  elev_new <- overtopping(peakflow, hydrograph, overtop_elev)
  if(elev_new < elev_old){peakflow_old <- peakflow 
  elev_old <- elev_new}
  print(c((100/iterations*i),peakflow_old,elev_old))
  if(elev_old == 0){break}
}
peakflows_4919[3,2] <- round(peakflow_old)

# all data reported in m^3/s
# scaled hydrographs, none include snowmelt

################################### Dam Elevation 4898.7 ##################################
overtop_elev <- 4898.7
peakflows_4898.7 <- data.frame(0)
peakflows_4898.7[1,1] <- "PMF_nosnow_4898.7"
peakflows_4898.7[2,1] <- "floodofrecord_4898.7"
peakflows_4898.7[3,1] <- "TREX_nosnow_4898.7"
hydrograph <- a
elev_old = 10
elev_new = 10
peakflow_old = 10000
peakflow = 10000
iterations = 1000
for(i in 1:iterations){
  peakflow = peakflow_old + sample(rnorm(100,mean = 0, sd = (500*elev_old)), 1)
  elev_new <- overtopping(peakflow, hydrograph, overtop_elev)
  if(elev_new < elev_old){peakflow_old <- peakflow 
  elev_old <- elev_new}
  print(c((100/iterations*i),peakflow_old,elev_old))
  if(elev_old == 0){break}
}
peakflows_4898.7[1,2] <- round(peakflow_old)

hydrograph <- b
elev_old = 10
elev_new = 10
peakflow_old = 10000
peakflow = 10000
iterations = 1000
for(i in 1:iterations){
  peakflow = peakflow_old + sample(rnorm(100,mean = 0, sd = (500*elev_old)), 1)
  elev_new <- overtopping(peakflow, hydrograph, overtop_elev)
  if(elev_new < elev_old){peakflow_old <- peakflow 
  elev_old <- elev_new}
  print(c((100/iterations*i),peakflow_old,elev_old))
  if(elev_old == 0){break}
}
peakflows_4898.7[2,2] <- round(peakflow_old)

#hydrograph <- d
elev_old = 10
elev_new = 10
peakflow_old = 10000
peakflow = 10000
iterations = 1000
for(i in 1:iterations){
  peakflow = peakflow_old + sample(rnorm(100,mean = 0, sd = (300*elev_old)), 1)
  hydrograph[,2] <- abs(as.numeric(unlist(d[which(abs(TREX_list_max[1:15]-peakflow)==min(abs(TREX_list_max[1:15]-peakflow)))])))  
  elev_new <- overtopping(peakflow, hydrograph, overtop_elev)
  if(elev_new < elev_old){peakflow_old <- peakflow 
  elev_old <- elev_new}
  print(c((100/iterations*i),peakflow_old,elev_old))
  if(elev_old == 0){break}
}
peakflows_4898.7[3,2] <- round(peakflow_old)

# all data reported in m^3/s
# scaled hydrographs, none include snowmelt

################################# All Data ###########################
allpeaks <- data.frame(0,check.names = FALSE, check.rows = FALSE)
allpeaks <- cbind(allpeaks, peakflows_4931.27)
allpeaks <- cbind(allpeaks, peakflows_4925)
allpeaks <- cbind(allpeaks, peakflows_4919)
allpeaks <- cbind(allpeaks, peakflows_4898.7)
############ Reservoir surface elevation for each hydrograph for the 1921 peak flow ###################

# set peak flow to peak flow of 1921 flood (flood of record)
peakflow <- max(floodofrecord[,2])*0.0283168 # convert to cubic meters
########### Flood of record
hydrograph_nosnow(peakflow, floodofrecord)
cmd_exec("echo data.in | flrout-kjr-edited.exe")
# load the output data file
output <- read.fwf(file = "data.out", 
                   widths = c(10,10,11,11,11,11,11,9), skip = 53, n = 500,
                   col.names = c('TIME', 'INFLOW', 'OUTFLOW', 'STORAGE', 'ELEVATION', 
                                 'WATERWAY 1', 'WATERWAY 2', 'RATING CURVE'))
floodofrecord_ts  <- output[,5]

########### PMF scaled to flood of record
hydrograph_nosnow(peakflow, PMFhydrograph_nosnow[1:720,1:2])
cmd_exec("echo data.in | flrout-kjr-edited.exe")
# load the output data file
output <- read.fwf(file = "data.out", 
                   widths = c(10,10,11,11,11,11,11,9), skip = 53, n = 500,
                   col.names = c('TIME', 'INFLOW', 'OUTFLOW', 'STORAGE', 'ELEVATION', 
                                 'WATERWAY 1', 'WATERWAY 2', 'RATING CURVE'))
PMFscaled_ts  <- output[,5]

########### TREX scaled to flood of record
duration <- PMFhydrograph_nosnow[1:720,1]
TREXtimeseries_nosnow <- abs(as.numeric(unlist(d[which(abs(TREX_list_max[1:15]-peakflow)==min(abs(TREX_list_max[1:15]-peakflow)))])))
TREXtimeseries_nosnow <- round(TREXtimeseries_nosnow)
TREXhydrograph <- matrix(data = NA, nrow = length(TREXtimeseries_nosnow), ncol = 2)
TREXhydrograph[,1] <- duration
TREXhydrograph[,2] <- TREXtimeseries_nosnow
hydrograph_nosnow(peakflow, TREXhydrograph)
cmd_exec("echo data.in | flrout-kjr-edited.exe")
# load the output data file
output <- read.fwf(file = "data.out", 
                   widths = c(10,10,11,11,11,11,11,9), skip = 53, n = 500,
                   col.names = c('TIME', 'INFLOW', 'OUTFLOW', 'STORAGE', 'ELEVATION', 
                                 'WATERWAY 1', 'WATERWAY 2', 'RATING CURVE'))
TREXscaled_ts  <- output[,5]