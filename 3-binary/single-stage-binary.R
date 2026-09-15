
# Single-Stage Binary Model

# Libraries
library(nimble)

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
rm('Xcov','scalingvalues')

st1<- Sys.time()

X = cbind(rep(1,Tobs),Xfull)
Y = matrix(yfull,Q,Tobs)

var.check=rep(0,Q)
for (i in 1:Q){
  var.check[i] <- var(Y[i,])
}
var.check

bp = dim(X)[2]

## MCMC in NIMBLE

library(nimble)
library(coda)

# set initial values
alpha = c(-Inf,0,Inf)
tau.z = rep(1,Q)
Z = 1*(Y - 0.5)

beta = matrix(NA,Q,bp) #### use lm estimates for initial values
for(q in 1:Q){
  beta[q,] = lm(Z[q,]~X[seq(q,(Q*Tobs),by=Q),]-1)$coefficients
  beta[q,which(is.na(beta[q,]))]=0
}

rho.Z = rep(.95,Q)

# need adjacency matrix
# A = adjacency matrix
A <- A.grid
A <- A[which(colnames(A.grid) %in% gridID),which(colnames(A.grid) %in% gridID)]

num<-colSums(A)

## Vectorize for NIMBLE
adj<-NULL
for(j in 1:Q){
  adj<-c(adj,which(A[j,]==1))
}
adj<-as.vector(adj)
num<-as.vector(num)
weights<-1+0*adj

## Set up MCMC
M.iter = 5000
M.burn = 0
M.thin = 5

mod_data=list(Y=Y, X=X)
mod_constants=list(Tobs=Tobs, Q=Q, bp=bp, cut=0, adj=adj, num=num, weights=weights)
mod_inits=list(beta=beta, l.rho.z=log(rho.Z/(1-rho.Z)), Z=1*(Y-.5))


model_code=nimbleCode({
  #Drought Variable
  for(q in 1:Q){
    #Latent Gaussian Variable
    mu[q,1] <- inprod(X[q,1:bp],beta[q,1:bp])
    Z[q,1] ~ dnorm(mu[q,1], tau = 1)
    Y[q,1] ~ dinterval(Z[q,1], cut)
    for(t in 2:Tobs){
      mu[q,t] <- inprod(X[Q*(t-1)+q,1:bp],beta[q,1:bp]) +
        rho.z[q]*(Z[q,(t-1)] - (inprod(X[Q*(t-2)+q,1:bp], beta[q,1:bp])))
      Z[q,t] ~ dnorm(mu[q,t], tau = 1)
      Y[q,t] ~ dinterval(Z[q,t], cut)
    }
  }
  
  ### ICAR prior distribution for each beta
  for(p in 1:bp){
    beta[1:Q,p] ~ dcar_normal(adj[], weights[], num[], tau.b[p], zero_mean=0)
    tau.b[p] ~ dgamma(.01,.01)
  }
  
  ### ICAR prior for logit(rho.z)
  l.rho.z[1:Q] ~ dcar_normal(adj[], weights[], num[], tau.lr, zero_mean=0)
  tau.lr ~ dgamma(.01,.01)
  
  for(q in 1:Q){
    rho.z[q] <- exp(l.rho.z[q])/(1+exp(l.rho.z[q]))
  }
} ## closes nimble code
)

nimble_model <- nimbleModel(model_code, mod_constants, mod_data, mod_inits)
compiled_model <- compileNimble(nimble_model,resetFunctions = TRUE)
mcmc_conf <- configureMCMC(nimble_model,monitors=c('beta','rho.z','Z'),
                           control=list(adaptive=TRUE,scale=0.1,adaptInterval=100,sliceMaxSteps=100000,maxContractions=100000,sliceWidth=1),
                           useConjugacy = TRUE)

### change Z sampler from RW to AF_slice for each location
for(q in 1:Q){
  print(q)
  Zvec = rep(0,Tobs)
  for(j in 1:Tobs){
    Zvec[j]=paste("Z[",q,", ",j, "]",sep="")
  }
  mcmc_conf$removeSamplers(Zvec)
  mcmc_conf$addSampler(target=Zvec,type='AF_slice',control=list(adaptive=TRUE,sliceWidths=rep(.5,Tobs),sliceMaxSteps=20000,maxContractions=50000))
}

nimble_mcmc<-buildMCMC(mcmc_conf)
compiled_mcmc<-compileNimble(nimble_mcmc, project = nimble_model,resetFunctions = TRUE)

st2 <- Sys.time()

samples=runMCMC(compiled_mcmc,inits=mod_inits,
                nchains = 1, nburnin=M.burn, niter = M.iter, samplesAsCodaMCMC = TRUE, thin=M.thin,
                summary = FALSE, WAIC = FALSE, progressBar=TRUE)

time.out.total <- Sys.time() - st1
time.out.run <- Sys.time() - st2

save(samples,time.out.total,time.out.run, file="SingleStageResults/MCMCoutputSingle.Rda")

if(0){
  # Run in blocks of and save
  total_iters <- 50000
  block_size <- 5000
  num_blocks <- total_iters / block_size
  
  for(i in 1:num_blocks) {
    # Run 5000 iterations, picking up where it left off
    compiled_mcmc$run(niter = block_size, thin = M.thin, reset = FALSE)
    
    # Get samples and convert to matrix
    samples <- as.matrix(compiled_mcmc$mvSamples)
    
    S <- floor(nrow(samples)/2)
    subs <- samples[S:nrow(samples),]
    
    # Save the block
    save(subs, file = paste0("SingleStageResults/MCMCoutputSingle_", i, ".RData"))
  }
}






