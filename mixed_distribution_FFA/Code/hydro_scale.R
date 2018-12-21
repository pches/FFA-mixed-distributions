#################################################################################
#
#  -file = "hydro_scale.R"   Code written June 2018
#  - Author: Kenneth Joel Roop-Eckart (kjr30@psu.edu)
#
#  - This function scales input hydrographs to an input peak flow. For more
#    information see Appendix B of the thesis.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# To use this function, simply source this file:
#   source("hydro_scale.R")
#
###################################################################################

########################################## Hydrograph scaling function #########################################
# scale hydrographs without a snowmelt contribution
hydrograph_nosnow <- function(peakflow, hydrograph){
  peakflow <- peakflow
  hydrograph <- hydrograph
  dt = 30*60
  
  PMFtimeseries_nosnow <- hydrograph[1:720,2]
  time <- hydrograph[1:720,1]
  vol <- dt*PMFtimeseries_nosnow
  
  peakflow_PMF_nosnow <- max(PMFtimeseries_nosnow)
  
  Scaled_Hydrograph <- PMFtimeseries_nosnow*(peakflow/peakflow_PMF_nosnow)/0.0283168 # converted to m^3/s
  scaledpmfnosnow <- PMFtimeseries_nosnow*(peakflow/peakflow_PMF_nosnow)/0.0283168 # converted to m^3/s
  # save scaled hydrograph
  #Scaled_Hydrograph[,2] <- timeseries*35.3147 # convert back to ft^3/s for FLROUT (which uses ft^3/s)
  x <- NA
  for(i in 1:length(time)){
    piece <- paste(time[i], round(Scaled_Hydrograph[i], digits = 0), sep = '\t')
    x <- append(x, piece, after = i-1)
  }
  Scaled_Hydrograph <- x[1:length(time)]
  Scaled_Hydrograph <- as.vector(Scaled_Hydrograph)
  Scaled_Hydrograph <- append(Scaled_Hydrograph, c('Pueblo Dam
Rainfall only PMF Run
Full PMF Routing
3100100021
2,1,.01,0
4880.49,4960,4864
3,Top of Dam
3,Embankment "wings"
2.67,0,4931.27,1200.5,0,0,0,0
2.62,0,4931.27,8480.00,0,0,0,0'), after = 0) # elevation where overtopping begins changed from 4925.25 to 4931.27, and from 4925 to 4931.27
  Scaled_Hydrograph <- append(Scaled_Hydrograph, c('9999 9999
4753.20	0
4754.00	5
4760.00	853
4770.00	5056
4780.00	11677
4790.00	20701
4800.00	32214
4810.00	46676
4820.00	64861
4830.00	86629
4840.00	112286
4850.00	141401
4860.00	174439
4870.00	212059
4880.00	254696
4890.00	303066
4900.00	357369
4910.00	418864
4920.00	488871
4924.00	519654
4924.50	523636
4925.00	527626
4925.50	531626
4926.00	535626
4930.00	569747
4935.00	604527
9999 9999
4880.49	5000
4890.00	5000
4898.70	5000
4900.00	5000
4902.00	10000
4904.00	22500
4906.00	37500
4908.00	56000
4910.00	77000
4912.00	100000
4914.00	125000
4916.00	153000
4920.00	212000
4924.00	270000
4928.00	328000
4930.00	357000
4935.00	429500
9999 9999'), after = length(Scaled_Hydrograph))
# save as .in file
write(Scaled_Hydrograph, file = 'data.in', ncolumns = 1, sep = '\t')
}
