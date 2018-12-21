#################################################################################
#
#  -file = "Thesis Figure Script Final.R"   Code written June 2018
#  - Author: Kenneth Joel Roop-Eckart (kjr30@psu.edu)
#
#  - This code creates pdf files of every plot in the thesis and saves
#    them to the "Figures" folder.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# To use this function, simply source this file:
#   source("Thesis Figure Script Final.R")
#
###################################################################################

elev_Q_PMF_nosnow <- c(allpeaks[1,9],allpeaks[1,7],allpeaks[1,5],allpeaks[1,3])
elev_Q_1921 <- c(allpeaks[2,9],allpeaks[2,7],allpeaks[2,5],allpeaks[2,3])
elev_Q_TREX_hydrograph <- c(allpeaks[3,9],allpeaks[3,7],allpeaks[3,5],allpeaks[3,3])

######################## Functions ##########################################
library(fExtremes)
pmgev = function(q,parameters) 
{
  location1 <- parameters[3]
  scale1 <- parameters[2]
  shape1 <- parameters[1]
  location2 <- parameters[6]
  scale2 <- parameters[5]
  shape2 <- parameters[4]
  alpha <- parameters[7]
  if(max(fExtremes::pgev(q, xi = shape2, beta = scale2, mu = location2))>1){
    out <- ((alpha)*fExtremes::pgev(q, xi = shape1, beta = scale1, mu = location1)+(1-alpha)*1)}
  else if(max(fExtremes::pgev(q, xi = shape1, beta = scale1, mu = location1))>1){
    out <- ((alpha)*1+(1-alpha)*fExtremes::pgev(q, xi = shape2, beta = scale2, mu = location2))}
  else(out <- ((1-alpha)*fExtremes::pgev(q, xi = shape1, beta = scale1, mu = location1)+
                 (alpha)*fExtremes::pgev(q, xi = shape2, beta = scale2, mu = location2)))
  return(out)
}

ptcev <- function(q, parameters){
  A1<-parameters[1]
  theta1<-parameters[2]
  A2<-parameters[3]
  theta2<-parameters[4]
  out<- (exp(-A1*exp(-(q)/theta1)-A2*exp(-(q)/theta2)))
  return(out)
}
LP3prob <- function(quantiles, parameters){
  tau <- parameters[1]
  alpha <- parameters[2]
  beta <- parameters[3]
  x <- log10(quantiles)
  
  out <- (beta*Igamma(alpha, (x-tau)/beta))/(gamma(alpha)*abs(beta))
  return(out)
}
pmlp3 = function(q, parameters) 
{
  mu1 <- parameters[1]
  sigma1 <- parameters[2]
  Q1 <- parameters[3]
  mu2 <- parameters[4]
  sigma2 <- parameters[5]
  Q2 <- parameters[6]
  alpha <- parameters[7]
  x <- q
  out <- (alpha*LP3prob(x, parameters[1:3]) + 
            (1-alpha)*LP3prob(x, parameters[4:6]))
  return(out)
}
######################## (1) return period plot with data ############################
max_return_period <- 1000
q <- seq(min(data),2*paleo_flow_upper, 100)
par(mfrow=c(1,1))
cols <- brewer.pal(6,"Dark2")
pdf("Figures/Thesis_Figure_returnperiodsData.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
par(cex = 1.3)
par(mar=c(5,5,4,2))
plot((length(data)+1)/(1:length(data)), sort(data, decreasing = TRUE), log = "x", xlim = c(1.1, max_return_period), 
     ylim = c(min(data), paleo_flow_upper), pch = 16, cex = 0, xlab = 'Return period (years)', ylab = expression('Annual peak discharge   '*'(m'^3/s*')'),
     xaxt = 'n')
axis(1, lwd = 1.5, at=10^(seq(-1,log10(max_return_period), by = 1)), label=parse(text=paste("10^", seq(-1,log10(max_return_period), by = 1), sep="")))


# add labels
abline(v = 11, lty = 2, lwd = 2)
#text(4,2000, labels = paste('Snowmelt?'))
#text(25,2000, labels = paste('Rainfall?'))

# paleo data
arrows((mean(paleo_age)+1)/(1:length(paleo)), sort(paleo_flow_lower, decreasing = TRUE), (mean(paleo_age)+1)/(1:length(paleo)), sort(paleo_flow_upper, decreasing = TRUE), length=0.05, angle=90, code=3)
arrows((paleo_age_lower+1)/(1:length(paleo)), sort(paleo, decreasing = TRUE), (paleo_age_upper+1)/(1:length(paleo)), sort(paleo, decreasing = TRUE), length=0.05, angle=90, code=3)
points((mean(paleo_age)+1)/(1:length(paleo)), sort(paleo, decreasing = TRUE), pch = 16, cex = 1, col = 'brown')

# historical floods
nums <- c(which(floods == hist[1]),which(floods == hist[2]),which(floods == hist[3]))

# uncertainty
arrows(1/(Phat), sort(floods, decreasing = TRUE)-1.96*sort(total_uncertainty, decreasing = TRUE), 1/(Phat), sort(floods, decreasing = TRUE)+1.96*sort(total_uncertainty, decreasing = TRUE), length=0.05, angle=90, code=3)

points(sort(1/(Phat), decreasing = FALSE)[nums], floods[nums], pch = 16)
# gage floods
points(sort(1/(Phat), decreasing = FALSE)[-nums], floods[-nums], pch = 16, col = 'gray')

legend('topleft', legend = c('Gage data', 'Historical floods', 'Paleoflood'), 
       col = c('gray', 'black', 'brown'), lty = c(NA,NA,NA), 
       pch = c(16,16,16), lwd = c(NA,NA,NA), bg = 'white')
dev.off()
dev.off()
######################## (2) return period plot with data + state of the art ############################
max_return_period <- 1000
q <- seq(min(data),2*paleo_flow_upper, 100)
par(mfrow=c(1,1))
cols <- brewer.pal(6,"Dark2")
pdf("Figures/Thesis_Figure_BestEstimateLP3.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
par(cex = 1.3)
par(mar=c(5,5,4,2))
plot((length(data)+1)/(1:length(data)), sort(data, decreasing = TRUE), log = "x", xlim = c(1.1, max_return_period), 
     ylim = c(min(data), paleo_flow_upper), pch = 16, cex = 0, xlab = 'Return period (years)', ylab = expression('Annual peak discharge   '*'(m'^3/s*')'),
     xaxt = 'n', yaxp = c(0, 5000, 5))
axis(1, lwd = 1.5, at=10^(seq(-1,log10(max_return_period), by = 1)), label=parse(text=paste("10^", seq(-1,log10(max_return_period), by = 1), sep="")))

# add labels
abline(v = 11, lty = 2, lwd = 2)
#text(4,2000, labels = paste('Snowmelt?'))
#text(25,2000, labels = paste('Rainfall?'))


# paleo data
arrows((mean(paleo_age)+1)/(1:length(paleo)), sort(paleo_flow_lower, decreasing = TRUE), (mean(paleo_age)+1)/(1:length(paleo)), sort(paleo_flow_upper, decreasing = TRUE), length=0.05, angle=90, code=3)
arrows((paleo_age_lower+1)/(1:length(paleo)), sort(paleo, decreasing = TRUE), (paleo_age_upper+1)/(1:length(paleo)), sort(paleo, decreasing = TRUE), length=0.05, angle=90, code=3)
points((mean(paleo_age)+1)/(1:length(paleo)), sort(paleo, decreasing = TRUE), pch = 16, cex = 1, col = 'brown')

# historical floods
nums <- c(which(floods == hist[1]),which(floods == hist[2]),which(floods == hist[3]))

# uncertianty
arrows(1/(Phat), sort(floods, decreasing = TRUE)-1.96*sort(total_uncertainty, decreasing = TRUE), 1/(Phat), sort(floods, decreasing = TRUE)+1.96*sort(total_uncertainty, decreasing = TRUE), length=0.05, angle=90, code=3)

points(sort(1/(Phat), decreasing = FALSE)[nums], floods[nums], pch = 16)
# gage floods
points(sort(1/(Phat), decreasing = FALSE)[-nums], floods[-nums], pch = 16, col = 'gray')


# LP3 fit
PLP3 <- LP3prob(q, MLE_LP3_optim$optim$bestmem)
lines(1/(1-PLP3), q, type = 'l', col = cols[2], lty = 1, lwd = 2)


legend('topleft', legend = c('Gage data', 'Historical floods', 'Paleoflood', 
                             'Log Pearson III'), 
       col = c('gray', 'black', 'brown',
               cols[2]), lty = c(NA,NA,NA,1), 
       pch = c(16,16,16,NA), lwd = c(NA,NA,NA,2), bg = 'white')
dev.off()
dev.off()
######################## (3) Best Estimate Return Period Plot for All Models ############################
max_return_period <- 1000
q <- seq(min(data),2*paleo_flow_upper, 100)
par(mfrow=c(1,1))
cols <- brewer.pal(3,"Dark2")
pdf("Figures/Thesis_Figure_BestEstimates.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
par(cex = 1.3)
par(mar=c(5,5,4,2))
plot((length(data)+1)/(1:length(data)), sort(data, decreasing = TRUE), log = "x", xlim = c(1.1, max_return_period), 
     ylim = c(min(data), paleo_flow_upper), pch = 16, cex = 0, xlab = 'Return period (years)', ylab = expression('Annual peak discharge   '*'(m'^3/s*')'),
     xaxt = 'n', yaxp = c(0, 5000, 5))
axis(1, lwd = 1.5, at=10^(seq(-1,log10(max_return_period), by = 1)), label=parse(text=paste("10^", seq(-1,log10(max_return_period), by = 1), sep="")))

# add labels
abline(v = 11, lty = 2, lwd = 2)
#text(4,800, labels = paste('Snowmelt?'))
#text(25,2000, labels = paste('Rainfall?'))

# paleo data
arrows((mean(paleo_age)+1)/(1:length(paleo)), sort(paleo_flow_lower, decreasing = TRUE), (mean(paleo_age)+1)/(1:length(paleo)), sort(paleo_flow_upper, decreasing = TRUE), length=0.05, angle=90, code=3)
arrows((paleo_age_lower+1)/(1:length(paleo)), sort(paleo, decreasing = TRUE), (paleo_age_upper+1)/(1:length(paleo)), sort(paleo, decreasing = TRUE), length=0.05, angle=90, code=3)
points((mean(paleo_age)+1)/(1:length(paleo)), sort(paleo, decreasing = TRUE), pch = 16, col = 'brown')

#gage data

# historical floods
nums <- c(which(floods == hist[1]),which(floods == hist[2]),which(floods == hist[3]))

# uncertianty
arrows(1/(Phat), sort(floods, decreasing = TRUE)-1.96*sort(total_uncertainty, decreasing = TRUE), 1/(Phat), sort(floods, decreasing = TRUE)+1.96*sort(total_uncertainty, decreasing = TRUE), length=0.05, angle=90, code=3)

points(sort(1/(Phat), decreasing = FALSE)[nums], floods[nums], pch = 16)
# gage floods
points(sort(1/(Phat), decreasing = FALSE)[-nums], floods[-nums], pch = 16, col = 'gray')

# LN2 fit
PLN2 <- plnorm(q, MLE_LN2_optim$optim$bestmem[1], MLE_LN2_optim$optim$bestmem[2])
lines(1/(1-PLN2), q, type = 'l', col = cols[1], lty = 2, lwd = 2)

# LP3 fit
PLP3 <- LP3prob(q, MLE_LP3_optim$optim$bestmem)
lines(1/(1-PLP3), q, type = 'l', col = cols[2], lty = 2, lwd = 2)

# GEV fit
PGEV <- pgev(q, MLE_GEV_optim$optim$bestmem[1], MLE_GEV_optim$optim$bestmem[3], MLE_GEV_optim$optim$bestmem[2])
lines(1/(1-PGEV), q, type = 'l', col = cols[3], lty = 2, lwd = 2)

# TCEV fit
PTCEV <- ptcev(q,MLE_TCEV_optim$optim$bestmem)
lines(1/(1-PTCEV), q, type = 'l', col = cols[1], lty = 1, lwd = 2)

# MLP3 fit
PMLP3 <- pmlp3(q, MLE_MLP3_optim$optim$bestmem)
lines(1/(1-PMLP3), q, type = 'l', col = cols[2], lty = 1, lwd = 2)

# MGEV fit
PMGEV <- pmgev(q, c(MLE_MGEV_optim$optim$bestmem[4:6], MLE_MGEV_optim$optim$bestmem[1:3],MLE_MGEV_optim$optim$bestmem[7]))
lines(1/(1-PMGEV), q, type = 'l', col = cols[3], lty = 1, lwd = 2)

legend('topleft', legend = c('Gage data', 'Historical floods', 'Paleoflood', 
                             'LN2', 'LP3', 'GEV',
                             'TCEV','Mixed LP3',"Mixed GEV"), 
       col = c('gray', 'black', 'brown',
               cols[1],cols[2],cols[3],cols[1],cols[2],cols[3]), lty = c(NA,NA,NA,2,2,2,1,1,1), 
       pch = c(16,16,16,NA,NA,NA,NA,NA,NA), lwd = c(NA,NA,NA,2,2,2,2,2,2), bg = 'white', cex = 0.9)

dev.off()
dev.off()
######################## (4) Best Estimate Return Period Plot for LP3 MGEV ############################
max_return_period <- 1000
q <- seq(min(data),2*paleo_flow_upper, 100)
par(mfrow=c(1,1))
cols <- brewer.pal(3,"Dark2")
pdf("Figures/Thesis_Figure_BestEstimates2.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
par(cex = 1.3)
par(mar=c(5,5,4,2))
plot((length(data)+1)/(1:length(data)), sort(data, decreasing = TRUE), log = "x", xlim = c(1.1, max_return_period), 
     ylim = c(min(data), paleo_flow_upper), pch = 16, cex = 0, xlab = 'Return period (years)', ylab = expression('Annual peak discharge   '*'(m'^3/s*')'),
     xaxt = 'n', yaxp = c(0, 5000, 5))
axis(1, lwd = 1.5, at=10^(seq(-1,log10(max_return_period), by = 1)), label=parse(text=paste("10^", seq(-1,log10(max_return_period), by = 1), sep="")))

# add labels
abline(v = 11, lty = 2, lwd = 2)
#text(4,2000, labels = paste('Snowmelt?'))
#text(25,2000, labels = paste('Rainfall?'))

# paleo data
arrows((mean(paleo_age)+1)/(1:length(paleo)), sort(paleo_flow_lower, decreasing = TRUE), (mean(paleo_age)+1)/(1:length(paleo)), sort(paleo_flow_upper, decreasing = TRUE), length=0.05, angle=90, code=3)
arrows((paleo_age_lower+1)/(1:length(paleo)), sort(paleo, decreasing = TRUE), (paleo_age_upper+1)/(1:length(paleo)), sort(paleo, decreasing = TRUE), length=0.05, angle=90, code=3)
points((mean(paleo_age)+1)/(1:length(paleo)), sort(paleo, decreasing = TRUE), pch = 16, col = 'brown')

# historical floods
nums <- c(which(floods == hist[1]),which(floods == hist[2]),which(floods == hist[3]))

# uncertianty
arrows(1/(Phat), sort(floods, decreasing = TRUE)-1.96*sort(total_uncertainty, decreasing = TRUE), 1/(Phat), sort(floods, decreasing = TRUE)+1.96*sort(total_uncertainty, decreasing = TRUE), length=0.05, angle=90, code=3)

points(sort(1/(Phat), decreasing = FALSE)[nums], floods[nums], pch = 16)
# gage floods
points(sort(1/(Phat), decreasing = FALSE)[-nums], floods[-nums], pch = 16, col = 'gray')

# LP3 fit
PLP3 <- LP3prob(q, MLE_LP3_optim$optim$bestmem)
lines(1/(1-PLP3), q, type = 'l', col = cols[2], lty = 1, lwd = 2)

# MGEV fit
PMGEV <- pmgev(q, c(MLE_MGEV_optim$optim$bestmem[4:6], MLE_MGEV_optim$optim$bestmem[1:3], MLE_MGEV_optim$optim$bestmem[7]))
lines(1/(1-PMGEV), q, type = 'l', col = cols[3], lty = 1, lwd = 2)

legend('topleft', legend = c('Gage data', 'Historical floods', 'Paleoflood', 
                             'LP3',"Mixed GEV"), 
       col = c('gray', 'black', 'brown',
               cols[2],cols[3]), lty = c(NA,NA,NA,1,1), 
       pch = c(16,16,16,NA,NA), lwd = c(NA,NA,NA,2,2), bg = 'white', cex = 0.9)

dev.off()
dev.off()

###################### (1) Calculating Return Period vs Reservoir Elevation ##################################
res_elev <- (c(4898.7,4919,4925,4931.27)-4898.7)*0.3048 # convert to meters

# Peak flows and reservoir elevations for each hydrograph
elev_Q_PMF_nosnow <- c(peakflows_4898.7[1,2],peakflows_4919[1,2],
                       peakflows_4925[1,2],peakflows_4931.27[1,2])
elev_Q_1921 <- c(peakflows_4898.7[2,2],peakflows_4919[2,2],
                 peakflows_4925[2,2],peakflows_4931.27[2,2])
elev_Q_TREX_hydrograph <- c(peakflows_4898.7[3,2],peakflows_4919[3,2],
                            peakflows_4925[3,2],peakflows_4931.27[3,2])

# Compute Likelihood of overtopping flow based on the Log Normal Distribution
probsLN2 <- matrix(data = NA, nrow = 3, ncol = 4)
probsLN2[1,] <- 1-plnorm(elev_Q_PMF_nosnow, meanlog = MLE_LN2_optim$optim$bestmem[1],
                         sdlog = MLE_LN2_optim$optim$bestmem[2])
probsLN2[2,] <- 1-plnorm(elev_Q_1921, meanlog = MLE_LN2_optim$optim$bestmem[1],
                         sdlog = MLE_LN2_optim$optim$bestmem[2]) 
probsLN2[3,] <- 1-plnorm(elev_Q_TREX_hydrograph, meanlog = MLE_LN2_optim$optim$bestmem[1],
                         sdlog = MLE_LN2_optim$optim$bestmem[2]) 

# Compute Likelihood of overtopping flow based on the Log Pearson Type 3 Distribution
probsLP3 <- matrix(data = NA, nrow = 3, ncol = 4)
probsLP3[1,] <- 1-LP3prob(elev_Q_PMF_nosnow, MLE_LP3_optim$optim$bestmem)
probsLP3[2,] <- 1-LP3prob(elev_Q_1921, MLE_LP3_optim$optim$bestmem)
probsLP3[3,] <- 1-LP3prob(elev_Q_TREX_hydrograph, MLE_LP3_optim$optim$bestmem)

# Compute Likelihood of overtopping flow based on the Generalized Extreme Value Distribution
probsGEV <- matrix(data = NA, nrow = 3, ncol = 4)
probsGEV[1,] <- 1-pgev(elev_Q_PMF_nosnow, xi = MLE_GEV_optim$optim$bestmem[1],
                       beta = MLE_GEV_optim$optim$bestmem[2], mu = MLE_GEV_optim$optim$bestmem[3])
probsGEV[2,] <- 1-pgev(elev_Q_1921, xi = MLE_GEV_optim$optim$bestmem[1],
                       beta = MLE_GEV_optim$optim$bestmem[2], mu = MLE_GEV_optim$optim$bestmem[3])
probsGEV[3,] <- 1-pgev(elev_Q_TREX_hydrograph, xi = MLE_GEV_optim$optim$bestmem[1],
                       beta = MLE_GEV_optim$optim$bestmem[2], mu = MLE_GEV_optim$optim$bestmem[3])

# Compute Likelihood of overtopping flow based on the Log Pearson Type 3 Distribution
probsTCEV <- matrix(data = NA, nrow = 3, ncol = 4)
probsTCEV[1,] <- 1-(exp(-MLE_TCEV_optim$optim$bestmem[1]*
                          exp(-(elev_Q_PMF_nosnow)/MLE_TCEV_optim$optim$bestmem[2])
                        -MLE_TCEV_optim$optim$bestmem[3]*
                          exp(-(elev_Q_PMF_nosnow)/MLE_TCEV_optim$optim$bestmem[4])))
probsTCEV[2,] <- 1-(exp(-MLE_TCEV_optim$optim$bestmem[1]*
                          exp(-(elev_Q_1921)/MLE_TCEV_optim$optim$bestmem[2])
                        -MLE_TCEV_optim$optim$bestmem[3]*
                          exp(-(elev_Q_1921)/MLE_TCEV_optim$optim$bestmem[4])))
probsTCEV[3,] <- 1-(exp(-MLE_TCEV_optim$optim$bestmem[1]*
                          exp(-(elev_Q_TREX_hydrograph)/MLE_TCEV_optim$optim$bestmem[2])
                        -MLE_TCEV_optim$optim$bestmem[3]*
                          exp(-(elev_Q_TREX_hydrograph)/MLE_TCEV_optim$optim$bestmem[4])))

# Compute Likelihood of overtopping flow based on the mixed Log Pearson Type 3 Distribution
probsMLP3 <- matrix(data = NA, nrow = 3, ncol = 4)
probsMLP3[1,] <- 1-pmlp3(elev_Q_PMF_nosnow, MLE_MLP3_optim$optim$bestmem)
probsMLP3[2,] <- 1-pmlp3(elev_Q_1921, MLE_MLP3_optim$optim$bestmem)
probsMLP3[3,] <- 1-pmlp3(elev_Q_TREX_hydrograph, MLE_MLP3_optim$optim$bestmem)

# Compute Likelihood of overtopping flow based on the mixed Generlized Extreme Value Distribution
probsMGEV <- matrix(data = NA, nrow = 3, ncol = 4)
probsMGEV[1,] <- 1-pmgev(elev_Q_PMF_nosnow, c(MGEVrun1$optim$bestmem[4:6], MGEVrun1$optim$bestmem[1:3], MGEVrun1$optim$bestmem[7]))
probsMGEV[2,] <- 1-pmgev(elev_Q_1921, c(MGEVrun1$optim$bestmem[4:6], MGEVrun1$optim$bestmem[1:3], MGEVrun1$optim$bestmem[7]))
probsMGEV[3,] <- 1-pmgev(elev_Q_TREX_hydrograph, c(MGEVrun1$optim$bestmem[4:6], MGEVrun1$optim$bestmem[1:3], MGEVrun1$optim$bestmem[7]))

############################### (2) Plotting all model fits ###################################
return_period_var <- 1000000000
par(mfrow=c(1,1))
cols <- brewer.pal(3,"Dark2")
cols <- brewer.pal(6,"YlOrRd")
cols2 <- brewer.pal(3,"Set1")
pdf("Figures/Thesis_Figure_Res_Return_Periodv1.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
par(cex = 1.2)
plot(probsLN2[1,], res_elev, log = "x", type = "n", xlim = c(100, return_period_var), 
     #ylim = c(min(dat, q), max(dat, q + 1.96 * sqrt(v[1:length(q)]))), 
     ylim = c(0, 11.5),
     #ylim = c(-2000, 50000),
     xaxt = 'n',
     xlab = "Return period (years)", 
     ylab = "Height above emergency spillway crest (m)")
axis(1, lwd = 1, at=10^(seq(-1,log10(return_period_var), by = 1)), label=parse(text=paste("10^", seq(-1,log10(return_period_var), by = 1), sep="")))
axis(2, lwd = 1)

# Overtopping return level data for PMF hydrograph
abline(h = max(res_elev), col = 'black', lty = 2, lwd = 2)
polygon(c(1/probsTCEV[3,],rev(1/probsLN2[2,])), c(res_elev, rev(res_elev)), 
        col = adjustcolor(cols[1], alpha.f=0.7), border = NA)
polygon(c(1/probsLN2[3,],rev(1/probsLN2[2,])), c(res_elev, rev(res_elev)), 
        col = adjustcolor(cols[2], alpha.f=0.7), border = NA)
polygon(c(1/probsLP3[3,],rev(1/probsLP3[2,])), c(res_elev, rev(res_elev)), 
        col = adjustcolor(cols[4], alpha.f=0.7), border = NA)
polygon(c(1/probsMLP3[3,],rev(1/probsMLP3[2,])), c(res_elev, rev(res_elev)), 
        col = cols[4], border = NA)
polygon(c(1/probsGEV[3,],rev(1/probsGEV[2,])), c(res_elev, rev(res_elev)), 
        col = adjustcolor(cols[5], alpha.f=0.7), border = NA)
polygon(c(1/probsMGEV[3,],rev(1/probsMGEV[2,])), c(res_elev, rev(res_elev)), 
        col = adjustcolor(cols[6], alpha.f=0.7), border = NA)

text(10*10^6,11,cex = 1,labels = "Meets regulations")
arrows(3*10^5,10.5,10*10^8,10.5,length = 0.25, code = 2, lwd = 3, col = cols2[3])
text(2*10^3,11,cex = 1,labels = "Not meet regulations")
arrows(2*10^5,10.5,6*10^1,10.5,length = 0.25, code = 2, lwd = 3, col = cols2[1])
arrows(131000,10.5,376000,10.5,length = 0.25, code = 0, lwd = 3, col = 'orange')

legend('bottomright', legend = c('Two Component Extreme Value','Log Normal', 'Log Pearson III',
                                 'Mixed Log Pearson III', 'Generalized Extreme Value', 
                                 'Mixed Generalized Extreme value',
                                 'Dam Crest Elevation', 'USBR Regulation Return Period'), 
       col = c(cols[1],cols[2],cols[4], cols[3], cols[5], cols[6], 'black', 'orange'), 
       lty = c(NA,NA,NA,NA,NA,NA,2,1), pch = c(
         15,15,15,15,15,15,NA,NA), lwd = c(
           2,2,2,2,2,2,2,2), bg = 'white', cex = 0.9)
dev.off()
dev.off()
############################### (3) Plotting only LP3 and MGEV ###################################
return_period_var <- 100000000
par(mfrow=c(1,1))
cols <- brewer.pal(3,"Dark2")
cols <- brewer.pal(6,"YlOrRd")
cols2 <- brewer.pal(3,"Set1")
pdf("Figures/Thesis_Figure_Res_Return_Periodv2.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
par(cex = 1.3)
plot(probsLP3[1,], res_elev, log = "x", type = "n", xlim = c(100, return_period_var), 
     ylim = c(0,11.5),
     xaxt = 'n',
     xlab = "Return period (years)", 
     ylab = "Height above emergency spillway (m)")
axis(1, lwd = 1.5, at=10^(seq(-1,log10(return_period_var), by = 1)), label=parse(text=paste("10^", seq(-1,log10(return_period_var), by = 1), sep="")))
axis(2, lwd = 1.5)

abline(h = max(res_elev), col = 'black', lty = 2, lwd = 2)

polygon(c(1/probsLP3[3,],rev(1/probsLP3[2,])), c(res_elev, rev(res_elev)), 
        col = adjustcolor(cols[4], alpha.f=1), border = NA)
polygon(c(1/probsMGEV[3,],rev(1/probsMGEV[2,])), c(res_elev, rev(res_elev)), 
        col = adjustcolor(cols[6], alpha.f=1), border = NA)

text(10*10^5.9,11,cex = 1,labels = "Meets regulations")
arrows(2*10^5,10.5,10*10^7.2,10.5,length = 0.25, code = 2, lwd = 3, col = cols2[3])
text(2.5*10^3,11,cex = 1,labels = "Not meet regulations")
arrows(2*10^5,10.5,6*10^1,10.5,length = 0.25, code = 2, lwd = 3, col = cols2[1])
arrows(131000,10.5,376000,10.5,length = 0.25, code = 0, lwd = 3, col = 'orange')

legend('bottomright', legend = c('LP3','Mixed GEV', 
                                 'Dam Crest',
                                 'USBR Regulations'), 
       col = c(cols[4], cols[6], 'black', 'orange'), 
       lty = c(NA,NA,2,1), pch = c(15,15,NA,NA), lwd = c(
         2,2,2,2), bg = 'white', cex = 1)
dev.off()
dev.off()

############################### (3) Plotting only LP3 and MGEV LINES ###################################
return_period_var <- 1000000000
par(mfrow=c(1,1))
cols <- brewer.pal(3,"Dark2")
cols <- brewer.pal(6,"YlOrRd")
cols2 <- brewer.pal(3,"Set1")
pdf("Figures/Thesis_Figure_Res_Return_Periodv2lines.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
par(cex = 1.3)
plot(probsLP3[1,], res_elev, log = "x", type = "n", xlim = c(100, return_period_var), 
     ylim = c(0, 11.5),
     xaxt = 'n',
     xlab = "Return period (years)", 
     ylab = "Height above emergency spillway crest (m)")
axis(1, lwd = 1.5, at=10^(seq(-1,log10(return_period_var), by = 1)), label=parse(text=paste("10^", seq(-1,log10(return_period_var), by = 1), sep="")))
axis(2, lwd = 1.5)

abline(h = max(res_elev), col = 'black', lty = 2, lwd = 2)

lines(1/probsLP3[2,], res_elev, col = adjustcolor(cols[4], alpha.f=1), type = 'l', lwd = 3)
lines(1/probsLP3[3,], res_elev, col = adjustcolor(cols[4], alpha.f=1), type = 'l', lwd = 3)
lines(1/probsLP3[1,], res_elev, col = adjustcolor(cols[4], alpha.f=1), type = 'l', lwd = 3)

lines(1/probsMGEV[1,], res_elev, col = adjustcolor(cols[6], alpha.f=1), type = 'l', lwd = 3)
lines(1/probsMGEV[2,], res_elev, col = adjustcolor(cols[6], alpha.f=1), type = 'l', lwd = 3)
lines(1/probsMGEV[3,], res_elev, col = adjustcolor(cols[6], alpha.f=1), type = 'l', lwd = 3)

text(10*10^6,11,cex = 1,labels = "Meets regulations")
arrows(3*10^5,10.5,10*10^8,10.5,length = 0.25, code = 2, lwd = 3, col = cols2[3])
text(2*10^3,11,cex = 1,labels = "Not meet regulations")
arrows(2*10^5,10.5,6*10^1,10.5,length = 0.25, code = 2, lwd = 3, col = cols2[1])
arrows(131000,10.5,376000,10.5,length = 0.25, code = 0, lwd = 3, col = 'orange')

legend('bottomright', legend = c('Log Pearson III','Mixed Generalized Extreme value', 
                                 'Dam Crest Elevation',
                                 'USBR Regulation Return Period'), 
       col = c(cols[4], cols[6], 'black', 'orange'), 
       lty = c(NA,NA,2,1), pch = c(15,15,NA,NA), lwd = c(
         2,2,2,2,2,2,5), bg = 'white', cex = 1)
dev.off()
dev.off()

################################ Hydrographs plot 2 (Overlap) ##################################
# run scaled hydrgraphs script before making this plot
TREXtimeseries_nosnow <- abs(as.numeric(unlist(d[which(abs(TREX_list_max[1:15]-peakflow)==min(abs(TREX_list_max[1:15]-peakflow)))])))
time <- 1:720
par(mfrow = c(1,1))
peakflow <- 100
cols <- brewer.pal(3,"Paired")
pdf("Figures/Thesis_Figure_Hydrographs_overlap.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
par(cex = 1.2)
plot(time/2, PMFhydrograph_nosnow[1:720,2]*(peakflow/max(PMFhydrograph_nosnow[1:720,2])), type = 'n', lty = 1,
     xlab = 'Time (hours)', ylab = 'Discharge (percent of peak)', xlim = c(50, 250))
lines(time-50, PMFhydrograph_nosnow[1:720,2]*(peakflow/max(PMFhydrograph_nosnow[1:720,2])), type = 'l', lty = 1, col = cols[2], lwd = 2)
#lines(time+100, TREXtimeseries_nosnow*(peakflow/TREX_peakflow_nosnow)+snowtimeseries, type = 'l', lty = 3)
lines(time-50, TREXtimeseries_nosnow*(peakflow/max(TREXtimeseries_nosnow)), type = 'l', lty = 1, col = cols[3], lwd = 2)
lines(time-50, floodofrecord[,2]*(peakflow/max(floodofrecord[,2])), type = 'l', lty = 1, col = cols[1], lwd = 2)

legend('topright', legend = c('scaled PMF', 'scaled TREX', "scaled 1921 flood"), 
       col = c(cols[2], cols[3], cols[1]), lty = c(1,1,1), lwd = c(2,2,2))
dev.off()
dev.off()

############################### Hydrographs plot 3 Cumulative Volume #################################
# run scaled hydrgraphs script before making this plot
# calculate cumulative volume delivered
timestep <- 60*30
vol_PMF <- cumsum(PMFhydrograph_nosnow[1:720,2]*(peakflow/max(PMFhydrograph_nosnow[1:720,2])))*timestep
vol_TREX <- cumsum(TREXtimeseries_nosnow*(peakflow/max(TREXtimeseries_nosnow)))*timestep
vol_1921 <- cumsum(floodofrecord[,2]*(peakflow/max(floodofrecord[,2])))*timestep

par(mfrow = c(1,1))
peakflow <- peakflow
cols <- brewer.pal(3,"Paired")
pdf("Figures/Thesis_Figure_Hydrographs_cumulative.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
par(cex = 1.2)
plot(time, vol_TREX/max(vol_TREX)*100, type = 'n', lty = 1,
     xlab = 'Time (hours)', ylab = 'Cumulative volume (percent)', xlim = c(50, 250))
lines(time-50, vol_PMF/max(vol_TREX)*100, type = 'l', lty = 1, col = cols[2], lwd = 2)
#lines(time+100, TREXtimeseries_nosnow*(peakflow/TREX_peakflow_nosnow)+snowtimeseries, type = 'l', lty = 3)
lines(time-50, vol_TREX/max(vol_TREX)*100, type = 'l', lty = 1, col = cols[3], lwd = 2)
lines(time-50, vol_1921/max(vol_TREX)*100, type = 'l', lty = 1, col = cols[1], lwd = 2)

legend('bottomright', legend = c('scaled PMF', 'scaled TREX', "scaled 1921 flood"), 
       col = c(cols[2], cols[3], cols[1]), lty = c(1,1,1), lwd = c(2,2,2))
dev.off()
dev.off()


################################### Convergence of MLP3 DEoptim ################################

cols <- brewer.pal(3,"Paired")
iterations <- length(MLP3run1$member$bestvalit)
par(mfrow = c(1,1))
pdf("Figures/Thesis_Figure_MLP3_convergence.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
par(cex = 1.2)
plot(1:iterations, -MLP3run1$member$bestvalit, type = 'l', col = cols[2], lwd = 2, ylim = c(-500,-650), 
     xlim = c(0, 2000), xlab = 'DEoptim iterations', ylab = 'Log likelihood')
lines(1:iterations, -MLP3run2$member$bestvalit, type = 'l', col = cols[3], lwd = 2)
legend('topright', legend = c('Mixed LP3 run 1', 'Mixed LP3 run 2'), 
       col = c(cols[2], cols[3]), lty = c(1,1), lwd = c(2,2))
dev.off()
dev.off()

################################### Convergence of MGEV DEoptim ################################

cols <- brewer.pal(5,"Paired")
par(mfrow = c(1,1))
pdf("Figures/Thesis_Figure_MGEV_convergence.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
par(cex = 1.2)
plot(1:length(MGEVrun1$member$bestvalit), -MGEVrun1$member$bestvalit, type = 'l', col = cols[1], lwd = 2, ylim = c(-500,-550), 
     xlim = c(0, 2500), xlab = 'DEoptim iterations', ylab = 'Log likelihood')
lines(1:length(MGEVrun2$member$bestvalit), -MGEVrun2$member$bestvalit, type = 'l', col = cols[2], lwd = 2)
lines(1:length(MGEVrun4$member$bestvali), -MGEVrun4$member$bestvalit, type = 'l', col = cols[3], lwd = 2)
lines(1:length(MGEVrun5$member$bestvalit), -MGEVrun5$member$bestvalit, type = 'l', col = cols[4], lwd = 2)

legend('topright', legend = c('Mixed GEV run 1', 'Mixed GEV run 2', 'Mixed GEV run 3', 'Mixed GEV run 4'), 
       col = c(cols[1], cols[2], cols[3], cols[4]), lty = c(1,1,1,1), lwd = c(2,2,2,2))
dev.off()
dev.off()


################################## Return period plot of MGEV MLEs ################################
max_return_period <- 1000
q <- seq(min(data),2*paleo_flow_upper, 100)
par(mfrow=c(1,1))
cols <- brewer.pal(4,"Dark2")
pdf("Figures/Thesis_Figure_rtMGEVruns.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
par(cex = 1.3)
plot((length(data)+1)/(1:length(data)), sort(data, decreasing = TRUE), log = "x", xlim = c(1.1, max_return_period), 
     ylim = c(min(data), paleo_flow_upper), pch = 16, cex = 0, xlab = 'Return period (years)', ylab = 'Peak flow (cms)',
     xaxt = 'n')
axis(1, lwd = 1.5, at=10^(seq(-1,log10(max_return_period), by = 1)), label=parse(text=paste("10^", seq(-1,log10(max_return_period), by = 1), sep="")))


# paleo data
arrows((mean(paleo_age+1))/(1:length(paleo)), sort(paleo_flow_lower, decreasing = TRUE), (mean(paleo_age+1)+1)/(1:length(paleo)), sort(paleo_flow_upper, decreasing = TRUE), length=0.05, angle=90, code=3)
arrows((paleo_age_lower+1)/(1:length(paleo)), sort(paleo, decreasing = TRUE), (paleo_age_upper+1)/(1:length(paleo)), sort(paleo, decreasing = TRUE), length=0.05, angle=90, code=3)
points((mean(paleo_age+1)+1)/(1:length(paleo)), sort(paleo, decreasing = TRUE), pch = 16, cex = 1, col = 'brown')

# historical floods
nums <- c(which(floods == hist[1]),which(floods == hist[2]),which(floods == hist[3]))

# uncertianty
arrows(1/(Phat), sort(floods, decreasing = TRUE)-1.96*sort(total_uncertainty, decreasing = TRUE), 1/(Phat), sort(floods, decreasing = TRUE)+1.96*sort(total_uncertainty, decreasing = TRUE), length=0.05, angle=90, code=3)

points(sort(1/(Phat), decreasing = FALSE)[nums], floods[nums], pch = 16)
# gage floods
points(sort(1/(Phat), decreasing = FALSE)[-nums], floods[-nums], pch = 16, col = 'gray')


# MGEV fit
PMGEV1 <- pmgev(q, c(MGEVrun1$optim$bestmem[4:6], MGEVrun1$optim$bestmem[1:3], MGEVrun1$optim$bestmem[7]))
lines(1/(1-PMGEV1), q, type = 'l', col = cols[1], lty = 1, lwd = 2)
PMGEV2 <- pmgev(q, c(MGEVrun2$optim$bestmem[4:6], MGEVrun2$optim$bestmem[1:3], MGEVrun2$optim$bestmem[7]))
lines(1/(1-PMGEV2), q, type = 'l', col = cols[2], lty = 1, lwd = 2)
PMGEV3 <- pmgev(q, c(MGEVrun4$optim$bestmem[4:6], MGEVrun4$optim$bestmem[1:3], MGEVrun4$optim$bestmem[7]))
lines(1/(1-PMGEV3), q, type = 'l', col = cols[3], lty = 1, lwd = 2)
PMGEV4 <- pmgev(q, c(MGEVrun5$optim$bestmem[4:6], MGEVrun5$optim$bestmem[1:3], MGEVrun5$optim$bestmem[7]))
lines(1/(1-PMGEV4), q, type = 'l', col = cols[4], lty = 1, lwd = 2)

legend('topleft', legend = c('Gage data', 'Historical floods', 'Paleoflood', 
                             'Mixed GEV fit 1','Mixed GEV fit 2','Mixed GEV fit 3','Mixed GEV fit 4'), 
       col = c('gray', 'black', 'brown',
               cols[1],cols[2],cols[3],cols[4]), lty = c(NA,NA,NA,1,1,1,1), 
       pch = c(16,16,16,NA,NA,NA,NA), lwd = c(NA,NA,NA,2,2,2,2), bg = 'white')
dev.off()
dev.off()

################ Reservoir surface elevation for 1921 flood #########################
# Plot reservoir surface elevation for each hydrograph scaled to 1921 flood
par(mfrow = c(1,1))
pdf("Figures/Thesis_Figure_reservoirresponse.pdf",width = 187*0.0393701, height = 187*(100/141)*0.0393701)
par(cex = 1.2)
plot(duration[1:length(floodofrecord_ts)]/2, floodofrecord_ts/3.2808399, type = 'l', col = cols[1], lwd = 2,
     ylim = c(4880, 4931.27)/3.2808399, ylab = 'Reservoir surface elevation (m)', xlab = 'Flood duration (hours)')
lines(duration[1:length(floodofrecord_ts)]/2, PMFscaled_ts/3.2808399, type = 'l', col = cols[2], lwd = 2)
lines(duration[1:length(floodofrecord_ts)]/2, TREXscaled_ts/3.2808399, type = 'l', col = cols[3], lwd = 2)
abline(h=4931.27/3.2808399, lty = 2, lwd = 2, col = 'red')
abline(h=4898.7/3.2808399, lty = 2, lwd = 2)

legend('topright', legend = c('1921 flood', 'scaled PMF hydrograph', 'scaled TREX hydrograph', 'Top of the flood pool', 'Reservoir crest'), 
       col = c(cols[1], cols[2], cols[3], 'black', 'red'), lty = c(1,1,1,2,2), lwd = c(2,2,2,2,2), bg = 'white')
dev.off()
dev.off()

