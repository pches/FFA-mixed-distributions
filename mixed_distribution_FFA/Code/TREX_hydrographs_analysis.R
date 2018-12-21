#################################################################################
#
#  -file = "TREX_hydrographs_analysis.R"   Code written June 2018
#  - Author: Kenneth Joel Roop-Eckart (kjr30@psu.edu)
#
#  - This code loads the TREX hydrographs from England et al. (2014)
#    for hydrograph scaling.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# To use this code, simply source this file:
#   source("TREX_hydrographs_analysis.R")
#
###################################################################################

# load PMF with snowmelt only hydrograph to subtract from the TREX hydrographs 
# to get TREX hydrographs with only rainfall inputs for scaling
PMFsnowonly <- read.csv(file = 'Data/pmfsnowonly.csv', sep = ',', header = FALSE, skip = 10)
PMFsnowonly<- PMFsnowonly$V2[1:720]

# load TREX hydrographs
TREX5_2 <- read.csv(file = 'Data/fs5-2sn.in', sep = '', header = FALSE, skip = 10)
TREX5_2 <- as.numeric((TREX5_2$V2[1:720]-PMFsnowonly)*0.0283168)
TREX5_3 <- read.csv(file = 'Data/fs5-3sn.in', sep = '', header = FALSE, skip = 10)
TREX5_3 <- as.numeric((TREX5_3$V2[1:720]-PMFsnowonly)*0.0283168)
TREX5_4 <- read.csv(file = 'Data/fs5-4sn.in', sep = '', header = FALSE, skip = 10)
TREX5_4 <- as.numeric((TREX5_4$V2[1:720]-PMFsnowonly)*0.0283168)
TREX5_5 <- read.csv(file = 'Data/fs5-5sn.in', sep = '', header = FALSE, skip = 10)
TREX5_5 <- as.numeric((TREX5_5$V2[1:720]-PMFsnowonly)*0.0283168)
TREX5_6 <- read.csv(file = 'Data/fs5-6sn.in', sep = '', header = FALSE, skip = 10)
TREX5_6 <- as.numeric((TREX5_6$V2[1:720]-PMFsnowonly)*0.0283168)
TREX5_7 <- read.csv(file = 'Data/fs5-7sn.in', sep = '', header = FALSE, skip = 10)
TREX5_7 <- as.numeric((TREX5_7$V2[1:720]-PMFsnowonly)*0.0283168)
TREX5_8 <- read.csv(file = 'Data/fs5-8sn.in', sep = '', header = FALSE, skip = 10)
TREX5_8 <- as.numeric((TREX5_8$V2[1:720]-PMFsnowonly)*0.0283168)
TREX5_9 <- read.csv(file = 'Data/fs5-9sn.in', sep = '', header = FALSE, skip = 10)
TREX5_9 <- as.numeric((TREX5_9$V2[1:720]-PMFsnowonly)*0.0283168)
TREX5_10 <- read.csv(file = 'Data/fs5-10sn.in', sep = '', header = FALSE, skip = 10)
TREX5_10 <- as.numeric((TREX5_10$V2[1:720]-PMFsnowonly)*0.0283168)
TREX5_11 <- read.csv(file = 'Data/fs5-11sn.in', sep = '', header = FALSE, skip = 10)
TREX5_11 <- as.numeric((TREX5_11$V2[1:720]-PMFsnowonly)*0.0283168)
TREX5_12 <- read.csv(file = 'Data/fs5-12sn.in', sep = '', header = FALSE, skip = 10)
TREX5_12 <- as.numeric((TREX5_12$V2[1:720]-PMFsnowonly)*0.0283168)
TREX5_13 <- read.csv(file = 'Data/fs5-13sn.in', sep = '', header = FALSE, skip = 10)
TREX5_13 <- as.numeric((TREX5_13$V2[1:720]-PMFsnowonly)*0.0283168)
TREX5_14 <- read.csv(file = 'Data/fs5-14sn.in', sep = '', header = FALSE, skip = 10)
TREX5_14 <- as.numeric((TREX5_14$V2[1:720]-PMFsnowonly)*0.0283168)
TREX5_15 <- read.csv(file = 'Data/fs5-15sn.in', sep = '', header = FALSE, skip = 10)
TREX5_15 <- as.numeric((TREX5_15$V2[1:720]-PMFsnowonly)*0.0283168)
TREX5_16 <- read.csv(file = 'Data/fs5-16snow.in', sep = '', header = FALSE, skip = 10)
TREX5_16 <- as.numeric((TREX5_16$V2[1:720]-PMFsnowonly)*0.0283168)

TREX12_1 <- read.csv(file = 'Data/fs12-1sn.in', sep = '', header = FALSE, skip = 10)
TREX12_1 <- as.numeric((TREX12_1$V2[1:720]-PMFsnowonly)*0.0283168)
TREX12_2 <- read.csv(file = 'Data/fs12-2sn.in', sep = '', header = FALSE, skip = 10)
TREX12_2 <- as.numeric((TREX12_2$V2[1:720]-PMFsnowonly)*0.0283168)
TREX12_3 <- read.csv(file = 'Data/fs12-3sn.in', sep = '', header = FALSE, skip = 10)
TREX12_3 <- as.numeric((TREX12_3$V2[1:720]-PMFsnowonly)*0.0283168)
TREX12_4 <- read.csv(file = 'Data/fs12-4sn.in', sep = '', header = FALSE, skip = 10)
TREX12_4 <- as.numeric((TREX12_4$V2[1:720]-PMFsnowonly)*0.0283168)
TREX12_5 <- read.csv(file = 'Data/fs12-5sn.in', sep = '', header = FALSE, skip = 10)
TREX12_5 <- as.numeric((TREX12_5$V2[1:720]-PMFsnowonly)*0.0283168)
TREX12_6 <- read.csv(file = 'Data/fs12-6sn.in', sep = '', header = FALSE, skip = 10)
TREX12_6 <- as.numeric((TREX12_6$V2[1:720]-PMFsnowonly)*0.0283168)
TREX12_7 <- read.csv(file = 'Data/fs12-7sn.in', sep = '', header = FALSE, skip = 10)
TREX12_7 <- as.numeric((TREX12_7$V2[1:720]-PMFsnowonly)*0.0283168)
TREX12_8 <- read.csv(file = 'Data/fs12-8sn.in', sep = '', header = FALSE, skip = 10)
TREX12_8 <- as.numeric((TREX12_8$V2[1:720]-PMFsnowonly)*0.0283168)
TREX12_9 <- read.csv(file = 'Data/fs12-9sn.in', sep = '', header = FALSE, skip = 10)
TREX12_9 <- as.numeric((TREX12_9$V2[1:720]-PMFsnowonly)*0.0283168)
TREX12_10 <- read.csv(file = 'Data/fs12-10sn.in', sep = '', header = FALSE, skip = 10)
TREX12_10 <- as.numeric((TREX12_10$V2[1:720]-PMFsnowonly)*0.0283168)
TREX12_11 <- read.csv(file = 'Data/fs12-11sn.in', sep = '', header = FALSE, skip = 10)
TREX12_11 <- as.numeric((TREX12_11$V2[1:720]-PMFsnowonly)*0.0283168)
TREX12_12 <- read.csv(file = 'Data/fs12-12sn.in', sep = '', header = FALSE, skip = 10)
TREX12_12 <- as.numeric((TREX12_12$V2[1:720]-PMFsnowonly)*0.0283168)
TREX12_13 <- read.csv(file = 'Data/fs12-13sn.in', sep = '', header = FALSE, skip = 10)
TREX12_13 <- as.numeric((TREX12_13$V2[1:720]-PMFsnowonly)*0.0283168)
TREX12_14 <- read.csv(file = 'Data/fs12-14sn.in', sep = '', header = FALSE, skip = 10)
TREX12_14 <- as.numeric((TREX12_14$V2[1:720]-PMFsnowonly)*0.0283168)

# create list of all TREX hydrographs
TREX_list<-list(TREX5_2,TREX5_3,TREX5_4,TREX5_5,TREX5_6,TREX5_7,TREX5_8,TREX5_9,TREX5_10,TREX5_11,TREX5_12,
                TREX5_13,TREX5_14,TREX5_15,TREX5_16,
                TREX12_1,TREX12_2,TREX12_3,TREX12_4,TREX12_5,TREX12_6,TREX12_7,TREX12_8,TREX12_9,TREX12_10,TREX12_11,
                TREX12_12,TREX12_13,TREX12_14)

# create list of peak flows for each TREX hydrograph
TREX_list_max <- c(max(TREX5_2),max(TREX5_3),max(TREX5_4),max(TREX5_5),max(TREX5_6),max(TREX5_7),
                   max(TREX5_8),max(TREX5_9),max(TREX5_10),max(TREX5_11),max(TREX5_12),
                   max(TREX5_13),max(TREX5_14),max(TREX5_15),max(TREX5_16),
                   max(TREX12_1),max(TREX12_2),max(TREX12_3),max(TREX12_4),max(TREX12_5),max(TREX12_6),
                   max(TREX12_7),max(TREX12_8),max(TREX12_9),max(TREX12_10),max(TREX12_11),max(TREX12_12),
                   max(TREX12_13),max(TREX12_14))