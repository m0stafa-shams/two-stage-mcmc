
# Stage One Binary Model

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

Q = length(gridID)
Tobs = nrow(training)/Q

# Create the design matrix
vars <- c("apcp.tr","soilm","tsoil")
Xcov = cbind(training[,vars])
J = ncol(Xcov)

Xm = rep(0,J)
Xs = rep(0,J)
for(j in 1:J){
  Xm[j] = scalingvalues$means[which(rownames(scalingvalues)==vars[j])]
  Xs[j] = scalingvalues$sds[which(rownames(scalingvalues)==vars[j])]
}

Xfull = matrix(NA,Q*Tobs,J)
for(j in 1:(J)){
  Xfull[,j] = (Xcov[,j]-Xm[j])/Xs[j]
}

yfull = as.numeric(training$bdrought)
rm('Xcov','training','scalingvalues')

# run as an array job so that the MCMC for each location is performed in parallel
args = Sys.getenv('SLURM_ARRAY_TASK_ID')
q = as.numeric(args[1])

II = seq(q,(Q*Tobs),by=Q)
X = cbind(rep(1,Tobs),Xfull[II,])
Y = yfull[II]

D = 1  ### 2 total drought levels

# set initial values for MCMC
alpha = c(-Inf,0,Inf)
tau.z = 1
Z = 1*(Y - 0.5)

if(var(Z)==0) {
  beta <- rep(NA,ncol(X))
  beta[1] = Z[1]
  beta[2:ncol(X)]=0
  rho.Z=0.98
} else {
  fit = arima(Z,order=c(1,0,0),xreg=X[,2:4], method="ML") ### use AR(1) fit to get initial values
  rho.Z = fit$coef[1]
  beta = fit$coef[2:5]
}

if(is.na(sum(beta))){
  beta[1] = Z[1]
  beta[2:ncol(X)]=0
  rho.Z=.98
}

bp = length(beta)

# Set up MCMC
M.iter = 100000
M.burn = 50000
M.thin = 10

library(nimble)
library(coda)

st <- Sys.time()

mod_data=list(Y=Y, X=X)
mod_constants=list(Tobs=Tobs, bp=bp, cut=0)
mod_inits=list(beta=beta, rho.z=rho.Z, Z=Z)

model_code=nimbleCode({
  #Drought Variable
  for(t in 1:Tobs){
    Y[t] ~ dinterval(Z[t], cut)
  }
  
  #Latent Gaussian Variable
  mu[1] <- inprod(X[1,1:bp],beta[1:bp])
  Z[1] ~ dnorm(mu[1], tau = 1)
  for(t in 2:Tobs){
    mu[t] <- inprod(X[t,1:bp],beta[1:bp]) +
      rho.z*(Z[(t-1)] - (inprod(X[(t-1),1:bp], beta[1:bp])))
    Z[t] ~ dnorm(mu[t], tau = 1)
  }
  
  #prior distribution for beta, b indexes each element of beta
  for (b in 1:bp){
    beta[b] ~ dnorm(0, tau = 1/9)
  }
  
  ##prior distribution for rho
  gamma ~ dlogis(0,1)
  rho.z <- exp(gamma)/(1+exp(gamma))
  
} ## closes nimble code
)

nimble_model <- nimbleModel(model_code, mod_constants, mod_data, mod_inits)
compiled_model <- compileNimble(nimble_model,resetFunctions = TRUE)
mcmc_conf <- configureMCMC(nimble_model, monitors=c('beta','rho.z','Z'),
                           control=list(adaptive=TRUE,scale=0.1,adaptInterval=100,sliceMaxSteps=100000,maxContractions=100000,sliceWidth=1),
                           useConjugacy = TRUE)

Zvec = rep(0,Tobs)
for(j in 1:Tobs){
  Zvec[j]=paste("Z[",j, "]",sep="")
}
mcmc_conf$removeSamplers(Zvec)
mcmc_conf$addSampler(target=Zvec,type='AF_slice',control=list(adaptive=TRUE,sliceWidths=rep(.5,Tobs),sliceMaxSteps=20000,maxContractions=50000))

nimble_mcmc<-buildMCMC(mcmc_conf)
compiled_mcmc<-compileNimble(nimble_mcmc, project = nimble_model,resetFunctions = TRUE)

samples=runMCMC(compiled_mcmc,inits=mod_inits,
                nchains = 1, nburnin=M.burn,niter = M.iter,samplesAsCodaMCMC = TRUE,thin=M.thin,
                summary = FALSE, WAIC = FALSE, progressBar=TRUE)

time.out = Sys.time()-st

zl = which(colnames(samples)=="Z[1]")
zu = which(colnames(samples)==paste("Z[",Tobs,"]",sep=""))
bl = which(colnames(samples) == "beta[1]")
bu = which(colnames(samples)== paste("beta[",bp,"]",sep=""))
rl = which(colnames(samples)=="rho.z")

MCMCout <- list("Z"=samples[,zl:zu],"beta"=samples[,bl:bu],"rho.Z"=samples[,rl])

save(MCMCout,time.out, file=paste("StageOneOutput/MCMCout.",q,".Rda",sep=""))


