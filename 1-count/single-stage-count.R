
# Single-Stage Count Model

code_start <- Sys.time()

sink(file = "runtime.txt", append = TRUE)
cat(paste("single-stage count code started at:", Sys.time(), "\n\n"))
sink()

library(nimble)
library(coda)

set.seed(1234)

st1 <- Sys.time()

################################################################################
# Load data
################################################################################

data <- read.csv("county_data.csv")

rm('county_data')

tab <- table(data$GEOID)

if(length(unique(tab)) != 1){
  stop("Counties do not all have the same number of time points.")
}

countyID <- unique(data$GEOID)

n <- length(countyID)             # number of counties
Tobs <- as.numeric(unique(tab))   # number of time points

cat("n =", n, "\n")
cat("Tobs =", Tobs, "\n")

################################################################################
# Build covariate design matrix
################################################################################

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

if(any(is.na(Xm)) || any(is.na(Xs))){
  stop("Missing covariate scaling values.")
}

if(any(Xs == 0)){
  stop("At least one covariate has zero standard deviation.")
}

Xfull <- matrix(NA, n * Tobs, J)
for(j in 1:(J)){
  Xfull[,j] = (Xcov[,j]-Xm[j])/Xs[j]
}

if(any(is.na(Xfull))){
  stop("Missing values detected in standardized covariates.")
}

rm('Xcov','scaling')

################################################################################
# Create arrays Y, pop, logYstar, and X for NIMBLE
################################################################################

yfull <- data$cases_c
pop.full <- as.numeric(data$pop_c)

if(any(is.na(yfull))) stop("Missing response values.")
if(any(yfull < 0)) stop("Counts must be nonnegative.")
if(any(yfull != floor(yfull))) stop("Counts must be integers.")

if(any(is.na(pop.full))) stop("Missing population values.")
if(any(pop.full <= 0)) stop("Population must be positive.")

############################################################
# Pooled log case rate
############################################################

pooled_log_rate <- log(
  (sum(yfull, na.rm = TRUE) + 0.5) /
    sum(pop.full, na.rm = TRUE)
)

cat("Pooled log case rate =", pooled_log_rate, "\n")

Y <- matrix(NA_integer_, n, Tobs)
pop <- matrix(NA, n, Tobs)
logYstar <- matrix(NA, n, Tobs)

# X is a 3D array: county x time x covariate including intercept
bp <- J + 1
X <- array(NA, dim = c(n, Tobs, bp))

for(i in 1:n){
  
  idx <- which(data$GEOID == countyID[i])
  
  Y[i, ] <- as.integer(yfull[idx])
  pop[i, ] <- pop.full[idx]
  
  Ystar_i <- pmax(Y[i, ], 0.5)
  logYstar[i, ] <- log(Ystar_i)
  
  X[i, , 1] <- 1
  for(j in 1:J){
    X[i, , j + 1] <- Xfull[idx, j]
  }
}

if(any(is.na(X))) stop("Missing values detected in X array.")
if(any(is.na(Y))) stop("Missing values detected in Y matrix.")
if(any(is.na(pop))) stop("Missing values detected in pop matrix.")

################################################################################
# Prepare adjacency matrix for ICAR priors
################################################################################

# A = adjacency matrix
A <- A.county

if(!all(dim(A) == c(n, n))){
  stop("A.county dimension does not match number of counties.")
}

if(any(diag(A) != 0)){
  stop("Adjacency matrix should have zero diagonal.")
}

if(!isTRUE(all.equal(A, t(A)))){
  stop("Adjacency matrix is not symmetric.")
}

numnns <- rowSums(A) # number of neighbors for each site

if(any(numnns == 0)){
  stop("At least one county has no neighbors. ICAR conditional update will fail.")
}

# Vectorize adjacency for NIMBLE dcar_normal
adj <- NULL
for(i in 1:n){
  adj <- c(adj, which(A[i, ] == 1))
}

adj <- as.integer(adj)
num <- as.integer(numnns)
weights <- rep(1, length(adj))

################################################################################
# Define custom ZINB distribution for NIMBLE
################################################################################

dZINB <- nimbleFunction(
  run = function(x = integer(0),
                 mu = double(0),
                 size = double(0),
                 pi = double(0),
                 log = integer(0, default = 0)) {
    
    returnType(double(0))
    
    if(mu <= 0 | size <= 0 | pi < 0 | pi > 1){
      if(log) return(-Inf)
      else return(0.0)
    }
    
    prob <- size / (size + mu)
    
    logNB <- dnbinom(x, prob = prob, size = size, log = TRUE)
    
    if(x == 0L){
      logLik <- log(pi + (1 - pi) * exp(logNB))
    } else {
      logLik <- log(1 - pi) + logNB
    }
    
    if(log) return(logLik)
    else return(exp(logLik))
  }
)

rZINB <- nimbleFunction(
  run = function(n = integer(0),
                 mu = double(0),
                 size = double(0),
                 pi = double(0)) {
    
    returnType(integer(0))
    
    if(n != 1) stop("rZINB only supports n = 1")
    
    isZero <- rbinom(1, size = 1, prob = pi)
    
    if(isZero == 1){
      return(0L)
    } else {
      prob <- size / (size + mu)
      return(as.integer(rnbinom(1, size = size, prob = prob)))
    }
  }
)

registerDistributions(list(
  dZINB = list(
    BUGSdist = "dZINB(mu, size, pi)",
    discrete = TRUE,
    types = c(
      "value = integer(0)",
      "mu = double(0)",
      "size = double(0)",
      "pi = double(0)"
    )
  )
))

################################################################################
# Initial values
################################################################################

beta_init <- matrix(
  0,
  nrow = n,
  ncol = bp
)

beta_init[, 1] <- pooled_log_rate

pi_init <- rep(0.05, n)

size_init <- rep(1, n)

gamma_init <- rep(atanh(0.3), n)

# Precision parameters for ICAR.
# tau.b = 1 / tau_beta^2, tau.g = 1 / tau_gamma^2
tau.b_init <- rep(1, bp)
tau.g_init <- 1

mod_inits <- list(
  beta = beta_init,
  gamma = gamma_init,
  pi = pi_init,
  size = size_init,
  tau.b = tau.b_init,
  tau.g = tau.g_init
)

################################################################################
# NIMBLE model
################################################################################

L <- length(adj)

mod_data <- list(
  Y = Y
)

mod_constants <- list(
  n = n,
  Tobs = Tobs,
  bp = bp,
  X = X,
  pop = pop,
  logYstar = logYstar,
  adj = adj,
  weights = weights,
  num = num
)

model_code <- nimbleCode({
  
  ########################################################################
  # ZINB likelihood with observation-driven AR(1) structure
  ########################################################################
  
  for(i in 1:n){
    
    rho[i] <- tanh(gamma[i])
    
    eta[i, 1] <- log(pop[i, 1]) + inprod(X[i, 1, 1:bp], beta[i, 1:bp])
    mu[i, 1] <- exp(eta[i, 1])
    
    Y[i, 1] ~ dZINB(mu = mu[i, 1], size = size[i], pi = pi[i])
    
    for(t in 2:Tobs){
      
      eta[i, t] <- log(pop[i, t]) + inprod(X[i, t, 1:bp], beta[i, 1:bp])
      
      mu[i, t] <- exp(
        eta[i, t] +
          rho[i] * (logYstar[i, t - 1] - eta[i, t - 1])
      )
      
      Y[i, t] ~ dZINB(mu = mu[i, t], size = size[i], pi = pi[i])
    }
  }
  
  ########################################################################
  # ICAR priors for regression coefficients beta[, p]
  ########################################################################
  
  for(p in 1:bp){
    
    beta[1:n, p] ~ dcar_normal(
      adj[],
      weights[],
      num[],
      tau.b[p],
      zero_mean = 0
    )
    
    tau.b[p] ~ dgamma(0.5, 0.5)
    tausq.b[p] <- 1 / tau.b[p]
  }
  
  ########################################################################
  # ICAR prior for gamma = atanh(rho)
  ########################################################################
  
  gamma[1:n] ~ dcar_normal(
    adj[],
    weights[],
    num[],
    tau.g,
    zero_mean = 0
  )
  
  tau.g ~ dgamma(0.5, 0.5)
  tausq.g <- 1 / tau.g
  
  ########################################################################
  # County-specific ZINB priors
  ########################################################################
  
  for(i in 1:n){
    pi[i] ~ dbeta(1, 1)
    size[i] ~ dgamma(2, 2)
  }
})

################################################################################
# Build and compile model
################################################################################

nimble_model <- nimbleModel(
  code = model_code,
  constants = mod_constants,
  data = mod_data,
  inits = mod_inits
)

compiled_model <- compileNimble(nimble_model, resetFunctions = TRUE)

################################################################################
# Configure MCMC
################################################################################

M.iter <- 100000
M.burn <- 20000
M.thin <- 20

monitors <- c(
  "beta",
  "gamma",
  "rho",
  "pi",
  "size",
  "tau.b",
  "tau.g",
  "tausq.b",
  "tausq.g"
)

mcmc_conf <- configureMCMC(
  nimble_model,
  monitors = monitors,
  useConjugacy = TRUE
)

################################################################################
# samplers
################################################################################

pi_nodes <- paste0("pi[", seq_len(n), "]")
size_nodes <- paste0("size[", seq_len(n), "]")

mcmc_conf$removeSamplers(
  c(
    "beta",
    "gamma",
    pi_nodes,
    size_nodes
  )
)

for(i in seq_len(n)){
  
  mcmc_conf$addSampler(
    target = pi_nodes[i],
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
    target = size_nodes[i],
    type = "slice",
    control = list(
      sliceWidth = 0.5,
      lower = 1e-8,
      sliceMaxSteps = 1000,
      maxContractions = 1000
    )
  )
  
  county_block <- c(
    paste0(
      "beta[",
      i,
      ",",
      seq_len(bp),
      "]"
    ),
    paste0(
      "gamma[",
      i,
      "]"
    )
  )
  
  mcmc_conf$addSampler(
    target = county_block,
    type = "AF_slice"
  )
}

mcmc_conf$printSamplers()

################################################################################
# Build, compile, and run MCMC
################################################################################

nimble_mcmc <- buildMCMC(mcmc_conf)

compiled_mcmc <- compileNimble(
  nimble_mcmc,
  project = nimble_model,
  resetFunctions = TRUE
)

st2 <- Sys.time()

samples <- runMCMC(
  compiled_mcmc,
  inits = mod_inits,
  nchains = 1,
  nburnin = M.burn,
  niter = M.iter,
  thin = M.thin,
  samplesAsCodaMCMC = TRUE,
  summary = FALSE,
  WAIC = FALSE,
  progressBar = TRUE
)

st3 <- Sys.time()

time.out.total <- st3 - st1
time.out.run <- st3 - st2

################################################################################
# Diagnostics and organize single-stage output
################################################################################

samples_mat <- as.matrix(samples)

beta_cols    <- grep("^beta\\[", colnames(samples_mat))
gamma_cols   <- grep("^gamma\\[", colnames(samples_mat))
rho_cols     <- grep("^rho\\[", colnames(samples_mat))
pi_cols      <- grep("^pi\\[", colnames(samples_mat))
size_cols    <- grep("^size\\[", colnames(samples_mat))

tau_b_cols   <- grep("^tau\\.b\\[", colnames(samples_mat))
tau_g_col    <- which(colnames(samples_mat) == "tau.g")

tausq_b_cols <- grep("^tausq\\.b\\[", colnames(samples_mat))
tausq_g_col  <- which(colnames(samples_mat) == "tausq.g")

# Diagnostics

ess_all <- coda::effectiveSize(samples)

mcmc_summary <- summary(samples)


############################################################
# Save single-stage output
############################################################

SingleStageOut <- list(
  
  ##########################################################
  # Parameter samples
  ##########################################################
  
  beta = samples_mat[, beta_cols, drop = FALSE],
  gamma = samples_mat[, gamma_cols, drop = FALSE],
  rho = samples_mat[, rho_cols, drop = FALSE],
  pi = samples_mat[, pi_cols, drop = FALSE],
  size = samples_mat[, size_cols, drop = FALSE],
  
  tau.b = samples_mat[, tau_b_cols, drop = FALSE],
  tau.g = samples_mat[, tau_g_col, drop = FALSE],
  tausq.b = samples_mat[, tausq_b_cols, drop = FALSE],
  tausq.g = samples_mat[, tausq_g_col, drop = FALSE],
  
  ##########################################################
  # Diagnostics
  ##########################################################
  
  ess_all = ess_all,
  mcmc_summary = mcmc_summary,
  
  ##########################################################
  # Runtime and MCMC settings
  ##########################################################
  
  M.iter = M.iter,
  M.burn = M.burn,
  M.thin = M.thin,
  
  time.out.total = time.out.total,
  time.out.run = time.out.run,
  
  ##########################################################
  # Data/model metadata
  ##########################################################
  
  vars = vars,
  n = n,
  Tobs = Tobs,
  bp = bp
)

if(!dir.exists("single_stage_output")){
  dir.create("single_stage_output")
}

save(
  SingleStageOut,
  file = "single_stage_output/single_stage_count.Rda"
)

################################################################################
# Runtime
################################################################################

sink(file = "runtime.txt", append = TRUE)
cat(paste("single-stage count code ended at:", Sys.time(), "\n"))
print(difftime(Sys.time(), code_start))
cat("\n")
sink()


