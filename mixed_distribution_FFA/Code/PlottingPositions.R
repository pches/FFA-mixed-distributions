# plotting positions for gage, historical, and paleo data


# method from:
# Hirsch, R. M., & Stedinger, J. R. (1987). Plotting positions for historical floods and their precision.
# Water Resources Research, 23(4), 715-727. https://doi.org/10.1029/WR023i004p00715

n = hist_record+length(data) # total record length

g = length(c(data,hist)) # total number of floods

X0 = mean(c(min(hist),max(data[which(data<min(hist))])))

k = length(which(c(data,hist)>X0)) # total number of floods above threshold X0

e = length(which(c(data)>X0)) # number of systematic record floods (gage record)

s = length(data) # length of systematic record (length gage record)



floods <- sort(c(data,hist), decreasing = FALSE)
rank <- 1:k
# flood error
total_uncertainty <- sort(c(gage_uncertainty,hist_uncertainty), decreasing = TRUE)

Phat<-vector(mode = 'numeric', length = g)

for(i in 1:k){
  Phat[i] <- (k/n)*(i/(k+1)) # for i=1....k
}

for(i in (k+1):g){
  Phat[i] <- (k/n)+((n-k)/n)*((i-k)/(s-e+1)) # for i=k+1...g
}

