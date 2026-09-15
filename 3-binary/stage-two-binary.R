
# Stage Two Binary Model

# Libraries
library(invgamma)
library(mvtnorm)

st <- Sys.time()

# load in all samples from the stage one model
# consolidate all stage one output into arrays 

load(paste('StageOneOutput/MCMCout.',1,'.Rda',sep=""))

n = 1354
T = ncol(MCMCout$Z)
M = nrow(MCMCout$Z)

Z.samp = array(NA,c(n,T,M))
beta.samp = array(NA,c(n,dim(MCMCout$beta)[2],M))
rho.Z.samp = matrix(NA,n,M)

rm('MCMCout')

for(q in 1:n){
  
  load(paste('StageOneOutput/MCMCout.',q,'.Rda',sep=""))
  
  Z.samp[q,,] = t(MCMCout$Z)
  beta.samp[q,,]=t(MCMCout$beta)
  rho.Z.samp[q,]=MCMCout$rho.Z
  rm('MCMCout')
}

# convert rho.Z to gamma samples
gamma.samp = log(rho.Z.samp/(1-rho.Z.samp))

## Set up MCMC for stage two

# set initial values
# take one draw for each site-specific parameter
Z = Z.samp[,,M]
beta = beta.samp[,,M]
gamma = gamma.samp[,M]

bp=dim(beta.samp)[2]

tausq.b = rep(1,bp)
tausq.g = 1

# Load data
# These csv files were generated previously by running 0-ScalingValues.R
data <- read.csv("USDMData.csv")
scalingvalues <- read.csv("scalingvalues.csv", row.names = 1)

bdrought <- 1*data$drought %in% c("D1","D2","D3","D4")
data <- data.frame(data, bdrought)
rm("bdrought")

# remove these 4 locations whose drought level never changes because they are over water
data = data[-which(data$grid %in% c("N78","W98","GG14","WW88")),] 

data = data[which(data$lon < -103.00),]

data = data[which(data$time>20190101),] 
training <- data[which(data$time<20210100),] 

training$timeID = as.numeric(as.factor(training$time))
gridID = training$grid[which(training$timeID==1)]

# need adjacency matrix for grid cells and number of neighbors of each
# A = adjacency matrix
A <- A.grid
A <- A[which(colnames(A.grid) %in% gridID),which(colnames(A.grid) %in% gridID)]

numnns = rowSums(A)
D = matrix(0,n,n)
diag(D)=numnns

M.iter = 100000
M.burn = 10000
M.thin = 5
M.out = (M.iter-M.burn)/M.thin

# create storage matrices
Z.out = matrix(NA,n,M.out)
beta.out = array(NA,c(n,bp,M.out))
rho.Z.out = matrix(NA,n,M.out)
tausq.b.out = matrix(NA,bp,M.out)
tausq.g.out = rep(NA,M.out)

progress_bar = txtProgressBar(min=0, max=M.iter, style = 3, char="=")

for(m in 1:M.iter){
  
  ## update tausq.b for each covariate
  for(p in 1:bp){
    tausq.b[p] = rinvgamma(1,0.01+n/2,0.01+1/2*t(beta[,p])%*%(D-A)%*%beta[,p])   
  }
  
  ## update tausq.g
  tausq.g = rinvgamma(1,0.01+n/2,0.01+1/2*t(gamma)%*%(D-A)%*%gamma)
  
  ## update the site-specific parameters for each location
  ## update full vector of site-specific parameters jointly
  for(i in 1:n){
    mi = sample(1:M,1)
    gamma.new = gamma.samp[i,mi]
    gammam = 1/numnns[i]*A[i,]%*%gamma
    betap.new = beta.samp[i,,mi]
    betam = 1/numnns[i]*A[i,]%*%beta
    ## make gamma prior ICAR
    R.new = dnorm(gamma.new,gammam,sd=sqrt(tausq.g/numnns[i]),log=TRUE)+sum(dnorm(betap.new,betam,sqrt(tausq.b/numnns[i]),log=TRUE))-dlogis(gamma.new,0,1,log=TRUE)-sum(dnorm(betap.new,0,sd=3,log=TRUE))
    R.old = dnorm(gamma[i],gammam,sd=sqrt(tausq.g/numnns[i]),log=TRUE)+sum(dnorm(beta[i,],betam,sqrt(tausq.b/numnns[i]),log=TRUE ))-dlogis(gamma[i],0,1,log=TRUE)-sum(dnorm(beta[i,],0,sd=3,log=TRUE))
    if(log(runif(1))<(R.new-R.old)){
      gamma[i] = gamma.new
      Z[i,] = Z.samp[i,,mi]
      beta[i,] = betap.new
    }
  }
  
  if(m>=M.burn & m/M.thin==floor(m/M.thin)){
    Z.out[,(m-M.burn)/M.thin] = Z[,T]
    beta.out[,,(m-M.burn)/M.thin] = beta
    rho.Z.out[,(m-M.burn)/M.thin] = exp(gamma)/(1+exp(gamma))
    tausq.b.out[,(m-M.burn)/M.thin] = tausq.b
    tausq.g.out[(m-M.burn)/M.thin] = tausq.g
  }
  
  setTxtProgressBar(progress_bar, value = m)
  
}
close(progress_bar)
time.out <- Sys.time()-st

save(beta.out,rho.Z.out,tausq.b.out,tausq.g.out,time.out,file="StageTwoOutput.Rda")



