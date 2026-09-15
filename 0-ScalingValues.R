
#############################################################################
#### Load in the raw data downloaded from 
#### Dryad Digital Repository at https://doi.org/10.5061/dryad.g1jwstqw7
#############################################################################

data <- read.csv("USDMDataAvg.csv")
data <- data[,-1] ## remove leading column

## apcp
min(data$apcp) ## minimum is zero
mean(is.na(data$apcp)) ## no missing

## ssrun
min(data$ssrun) ## minimum is zero
mean(is.na(data$ssrun)) ## no missing

## snod, smom, smowc, weasd not worth re-expressing

## check VPD, make transformed.
hist(data$vpd)
data$vpd.tr = log(data$vpd + 1)
hist(data$vpd.tr)

## Raising to (1/4) power seems like a reasonable addition
data$apcp.tr = (data$apcp)^(1/4)
data$ssrun.tr = (data$ssrun)^(1/4)

## Make a 28-day streamflow variable here
stream28 = data[,20]
stream28[is.na(stream28)] = data[is.na(stream28),33]
stream28[is.na(stream28)] = data[is.na(stream28),35]
stream28[is.na(stream28)] = data[is.na(stream28),37]
data = cbind(data,stream28)

rm(stream28)

n.cols <- dim(data)[2]
colnames(data)

par(mfrow=c(2,2))
for (i in 6:n.cols){
  hist(data[,i],main=colnames(data)[i])
}

means <- apply(data[,6:n.cols],2,mean, na.rm=TRUE)
sds <- apply(data[,6:n.cols],2,sd, na.rm=TRUE)
means
sds

scalingvalues = data.frame(means, sds)
rownames(scalingvalues) = colnames(data[6:n.cols])

write.csv(scalingvalues, "scalingvalues.csv")
write.csv(data, "USDMData.csv")
