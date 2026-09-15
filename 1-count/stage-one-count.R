
# Stage One Count Model

code_start <- Sys.time()

sink(file = "runtime.txt", append = TRUE)
cat(paste("code started at: ", Sys.time(), "\n"))
cat("\n")
sink()

set.seed(1234)

st1 <- Sys.time()

# Load data

data <- read.csv("county_data.csv")

rm('county_data')

tab <- table(data$GEOID)

if(length(unique(tab)) != 1){
  stop("Counties do not all have the same number of time points.")
}

countyID = unique(data$GEOID)
Q = length(countyID)  # Q = 227
Tobs = nrow(data)/Q   # Tobs = 1174

# Create the design matrix

colnames(data)
means <- apply(data[, !colnames(data) %in% c("time", "GEOID", "drought_c", "seas")],2,mean, na.rm=TRUE)
sds <- apply(data[, !colnames(data) %in% c("time", "GEOID", "drought_c", "seas")],2,sd, na.rm=TRUE)

scaling = data.frame(means, sds)
rownames(scaling) = colnames(data[, !colnames(data) %in% c("time", "GEOID", "drought_c", "seas")])

vars <- c("apcp_c","soilm_c","tsoil_c")

if(!all(vars %in% colnames(data))){
  stop("At least one covariate in vars is missing from data.")
}

Xcov <- as.matrix(data[, vars])
J <- ncol(Xcov)

Xm = rep(0,J)
Xs = rep(0,J)
for(j in 1:J){
  Xm[j] = scaling$means[which(rownames(scaling)==vars[j])]
  Xs[j] = scaling$sds[which(rownames(scaling)==vars[j])]
}

Xfull = matrix(NA,Q*Tobs,J)
for(j in 1:(J)){
  Xfull[,j] = (Xcov[,j]-Xm[j])/Xs[j]
}

# Full response and population vectors
yfull <- data$cases_c
pop.full <- data$pop_c

# Pooled log case rate for the intercept prior
pooled_log_rate <- log(
  (sum(yfull) + 0.5) /
    sum(pop.full)
)

cat("Pooled log case rate =", pooled_log_rate, "\n")

rm('Xcov','scaling')

# run as an array job so that the MCMC for each county is performed in parallel

args <- Sys.getenv("SLURM_ARRAY_TASK_ID")

if(args == ""){
  stop("SLURM_ARRAY_TASK_ID is empty.")
}

q <- as.numeric(args)

if(is.na(q) || q < 1 || q > Q){
  stop("Invalid array index q.")
}

set.seed(1234 + q)

idx <- which(data$GEOID == countyID[q])

X = cbind(rep(1,Tobs),Xfull[idx,])

Y = yfull[idx]
Tobs <- length(Y)

pop <- as.numeric(pop.full[idx])

if(any(is.na(pop))) stop("Missing population values.")
if(any(pop <= 0)) stop("Population must be positive.")

if(any(is.na(X))) stop("Missing values detected in the covariate design matrix (X).")

bp = ncol(X)

# Set up MCMC

M.iter <- 100000
M.burn <- 20000
M.thin <- 20

library(nimble)
library(coda)

# Define custom ZINB distribution for NIMBLE
dZINB <- nimbleFunction(
  run = function(x = integer(0), mu = double(0), size = double(0), 
                 pi = double(0), log = integer(0, default = 0)) {
    returnType(double(0))
    
    if(mu <= 0 | size <= 0 | pi < 0 | pi > 1){
      if(log) return(-Inf)
      else return(0.0)
    }
    
    # NB prob parameter
    prob <- size / (size + mu)
    
    logNB <- dnbinom(x, prob = prob, size = size, log = TRUE)
    
    if (x == 0L) {
      # Log-Likelihood for zero: P(Y=0) = pi + (1-pi) * NB(0)
      logLik <- log(pi + (1 - pi) * exp(logNB))
    } else {
      # Log-Likelihood for non-zero: P(Y=y>0) = (1-pi) * NB(y)
      logLik <- log(1 - pi) + logNB
    }
    
    if (log) return(logLik)
    else return(exp(logLik))
  }
)

# Define the random number generator
rZINB <- nimbleFunction(
  run = function(n = integer(0), mu = double(0), 
                 size = double(0), pi = double(0)) {
    returnType(integer(0))
    
    if (n != 1) stop("rZINB only supports n = 1")
    
    isZero <- rbinom(1, size = 1, prob = pi)
    
    if (isZero == 1) {
      return(0L)
    } else {
      prob <- size / (size + mu)
      return(as.integer(rnbinom(1, size = size, prob = prob)))
    }
  }
)

# Register the distribution so it can be used in nimbleCode
registerDistributions(list(
  dZINB = list(
    BUGSdist = "dZINB(mu, size, pi)",
    discrete = TRUE,
    types = c("value = integer(0)", "mu = double(0)", 
              "size = double(0)", "pi = double(0)")
  )
))


if(any(is.na(Y))) stop("Missing response values in Y.")
if(any(Y < 0)) stop("Counts must be nonnegative.")
if(any(Y != floor(Y))) stop("Counts must be integers.")

Y <- as.integer(Y)

Ystar <- pmax(Y, 0.5)
logYstar <- log(Ystar)

beta_init <- rep(0, bp)
beta_init[1] <- pooled_log_rate

pi_init <- 0.05

# NIMBLE model

mod_data=list(Y=Y)

mod_constants=list(X=X, pop=pop, 
                   Tobs=Tobs, bp=bp,
                   logYstar=logYstar,
                   pooled_log_rate = pooled_log_rate)

mod_inits <- list(
  beta = beta_init,
  pi = pi_init,
  size = 1,
  gamma = atanh(0.3)
)

model_code <- nimbleCode({
  
  gamma ~ dlogis(location = 0, scale = 0.5)
  rho <- tanh(gamma)

  eta[1] <- log(pop[1]) + inprod(X[1, 1:bp], beta[1:bp])
  mu[1] <- exp(eta[1])
  Y[1] ~ dZINB(mu = mu[1], size = size, pi = pi)
  
  for(t in 2:Tobs){
    eta[t] <- log(pop[t]) + inprod(X[t, 1:bp], beta[1:bp])
    mu[t] <- exp(eta[t] + rho * (logYstar[t-1] - eta[t-1]))
    
    Y[t] ~ dZINB(mu = mu[t], size = size, pi = pi)
  }
  
  beta[1] ~ dnorm(
    pooled_log_rate,
    tau = 1 / 25
  )
  
  for(b in 2:bp){
    beta[b] ~ dnorm(
      0,
      tau = 1 / 4
    )
  }
  
  pi ~ dbeta(1, 1)
  size ~ dgamma(2, 2)

})

nimble_model <- nimbleModel(code=model_code, 
                            constants=mod_constants, 
                            data=mod_data, 
                            inits=mod_inits)

compiled_model <- compileNimble(nimble_model,resetFunctions = TRUE)

monitors <- c("beta", "pi", "size", "rho", "gamma")

mcmc_conf <- configureMCMC(
  nimble_model,
  monitors = monitors,
  useConjugacy = TRUE
)

# samplers 

beta_nodes <- paste0(
  "beta[",
  seq_len(bp),
  "]"
)

mcmc_conf$removeSamplers(
  c(
    beta_nodes,
    "gamma",
    "pi",
    "size"
  )
)

mcmc_conf$addSampler(
  target = "pi",
  type = "slice",
  control = list(
    sliceWidth = 0.1,
    lower = 1e-8,
    upper = 1 - 1e-8,
    sliceMaxSteps = 500,
    maxContractions = 500
  )
)

mcmc_conf$addSampler(
  target = "size",
  type = "slice",
  control = list(
    sliceWidth = 0.5,
    lower = 1e-8,
    sliceMaxSteps = 1000,
    maxContractions = 1000
  )
)

mcmc_conf$addSampler(
  target = c(
    beta_nodes,
    "gamma"
  ),
  type = "AF_slice"
)

mcmc_conf$printSamplers()

nimble_mcmc<-buildMCMC(mcmc_conf)
compiled_mcmc<-compileNimble(nimble_mcmc, project = nimble_model,
                             resetFunctions = TRUE)

st2 <- Sys.time()

samples<-runMCMC(compiled_mcmc,inits=mod_inits,
                nchains = 1, nburnin=M.burn, niter = M.iter,
                thin=M.thin, samplesAsCodaMCMC = TRUE,
                summary = FALSE, WAIC = FALSE, progressBar=TRUE)

st3 <- Sys.time()

time.out.total <- st3 - st1
time.out.run <- st3 - st2

# convert samples to matrix

samples_mat <- as.matrix(samples)

beta_cols  <- grep("^beta\\[", colnames(samples_mat))
gamma_col  <- which(colnames(samples_mat) == "gamma")
rho_col    <- which(colnames(samples_mat) == "rho")
pi_col     <- which(colnames(samples_mat) == "pi")
size_col   <- which(colnames(samples_mat) == "size")

# Diagnostics

ess_all <- coda::effectiveSize(samples)

mcmc_summary <- summary(samples)

# Save stage one output

MCMCout <- list(
  beta = samples_mat[, beta_cols, drop = FALSE],
  gamma = samples_mat[, gamma_col],
  rho = samples_mat[, rho_col],
  pi = samples_mat[, pi_col],
  size = samples_mat[, size_col],
  
  ess_all = ess_all,
  mcmc_summary = mcmc_summary,

  M.iter = M.iter,
  M.burn = M.burn,
  M.thin = M.thin,
  
  time.out.total = time.out.total,
  time.out.run = time.out.run,

  Y = Y,
  q = q,
  county = countyID[q],
  vars = vars,
  bp = bp,
  Tobs = Tobs
)

if(!dir.exists("stage1_outputs_count")){
  dir.create("stage1_outputs_count")
}

save(
  MCMCout,
  file = paste0(
    "stage1_outputs_count/stage1_county_",
    q,
    ".Rda"
  )
)

# Runtime
sink(file = "runtime.txt", append = TRUE)
cat(paste("code ended at: ", Sys.time(), "\n"))
print(difftime(Sys.time(), code_start))
cat("\n")
sink()


