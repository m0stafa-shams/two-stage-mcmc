
# Stage One Continuous Model

code_start <- Sys.time()

sink(file = "runtime.txt", append = TRUE)
cat(paste("code started at: ", Sys.time(), "\n"))
cat("\n")
sink()

# Libraries
library(nimble)
library(MASS)
library(coda)
library(data.table)

set.seed(123)

st1 <- Sys.time()

############################################################
# Load data
############################################################

# This csv file was generated previously by running 0-ScalingValues.R
data <- fread("USDMData.csv")

data <- as.data.frame(data)
data <- data[, -1]

means <- apply(data[, 6:ncol(data)], 2, mean, na.rm = TRUE)
sds   <- apply(data[, 6:ncol(data)], 2, sd, na.rm = TRUE)

scalingvalues <- data.frame(means, sds)
rownames(scalingvalues) <- colnames(data[, 6:ncol(data)])

# Remove grids with no change in drought status
data <- data[-which(data$grid %in% c("N78", "W98", "GG14", "WW88")), ]

# Keep selected region and time period
training <- data[
  data$lon < -112 &
    data$lat > 31 &
    data$lat < 38 &
    data$time > 20190101 &
    data$time < 20210100,
]

training$timeID <- as.numeric(as.factor(training$time))

gridID <- training$grid[training$timeID == 1]

Q <- length(gridID)

coords <- training[training$timeID == 1, c("lon", "lat")]

Tobs <- nrow(training) / Q

cat("Q =", Q, "\n")
cat("Tobs =", Tobs, "\n")

############################################################
#  Standardize covariates
############################################################

# Response:
# soilm = weekly average soil moisture content
#
# Covariates:
# apcp  = weekly total precipitation
# pevap = weekly average potential evaporation
# tsoil = weekly average soil temperature

vars <- c("apcp", "pevap", "tsoil")

Xcov <- cbind(training[, vars])

J <- ncol(Xcov)

Xm <- rep(0, J)
Xs <- rep(0, J)

for(j in 1:J){
  Xm[j] <- scalingvalues$means[which(rownames(scalingvalues) == vars[j])]
  Xs[j] <- scalingvalues$sds[which(rownames(scalingvalues) == vars[j])]
}

Xfull <- matrix(NA, Q * Tobs, J)

for(j in 1:J){
  Xfull[, j] <- (Xcov[, j] - Xm[j]) / Xs[j]
}

# Response
yfull <- training[, "soilm"]

cat("\nSummary of soilm:\n")
print(summary(yfull))

if(any(yfull <= 0, na.rm = TRUE)){
  stop("Gamma response soilm must be strictly positive.")
}

############################################################
# run as an array job so that the MCMC for each location is performed in parallel
############################################################

args <- Sys.getenv("SLURM_ARRAY_TASK_ID")

q <- as.numeric(args[1])

if(is.na(q) || q < 1 || q > Q){
  stop("Invalid q. Check SLURM_ARRAY_TASK_ID.")
}

cat("Running Stage 1 for location q =", q, "\n")
cat("Grid ID =", gridID[q], "\n")

# Since training is ordered by timeID and grid, location q appears
# at rows q, q + Q, q + 2Q, ..., q + (Tobs - 1)Q.
II <- seq(q, Q * Tobs, by = Q)

X <- cbind(1, Xfull[II, ])

Y <- yfull[II]

logY <- log(Y)

bp <- ncol(X)

############################################################
#  Initial values
############################################################

beta_init <- rep(0, bp)

gamma_rho_init <- 0

rho_init <- 0

alpha_init <- 2

############################################################
#  MCMC settings
############################################################

M.iter <- 200000
M.burn <- 50000
M.thin <- 50

############################################################
#  NIMBLE data, constants, and initial values
############################################################

mod_data <- list(
  Y = Y
)

mod_constants <- list(
  Tobs = Tobs,
  bp = bp,
  X = X,
  logY = logY
)

mod_inits <- list(
  beta = beta_init,
  gamma_rho = gamma_rho_init,
  alpha = alpha_init
)

rm(data, Xcov, training, scalingvalues)

############################################################
#  NIMBLE model
############################################################

model_code <- nimbleCode({
  
  ##########################################################
  # Priors
  ##########################################################
  
  for(p in 1:bp){
    beta[p] ~ dnorm(0, sd = 10)
  }
  
  # Transformed AR parameter
  gamma_rho ~ dnorm(0, sd = 10)
  
  rho <- 2 / (1 + exp(-2 * gamma_rho)) - 1
  
  # Gamma shape parameter
  alpha ~ dgamma(0.01, 0.01)
  
  ##########################################################
  # Data model
  ##########################################################
  
  # Time t = 1
  xb[1] <- inprod(X[1, 1:bp], beta[1:bp])
  
  eta_mean[1] <- xb[1]
  
  mu[1] <- exp(eta_mean[1])
  
  rateY[1] <- alpha / mu[1]
  
  Y[1] ~ dgamma(shape = alpha, rate = rateY[1])
  
  # Time t = 2,...,Tobs
  for(t in 2:Tobs){
    
    xb[t] <- inprod(X[t, 1:bp], beta[1:bp])
    
    resid_prev[t] <- logY[t - 1] - xb[t - 1]
    
    eta_mean[t] <- xb[t] + rho * resid_prev[t]
    
    mu[t] <- exp(eta_mean[t])
    
    rateY[t] <- alpha / mu[t]
    
    Y[t] ~ dgamma(shape = alpha, rate = rateY[t])
  }
})

############################################################
#  Build and compile model
############################################################

nimble_model <- nimbleModel(
  code = model_code,
  constants = mod_constants,
  data = mod_data,
  inits = mod_inits
)

compiled_model <- compileNimble(
  nimble_model,
  resetFunctions = TRUE
)

############################################################
#  Configure MCMC
############################################################

mcmc_conf <- configureMCMC(
  nimble_model,
  monitors = c(
    "beta",
    "gamma_rho",
    "rho",
    "alpha"
  ),
  useConjugacy = TRUE
)

############################################################
#  print samplers
############################################################

cat("\nSamplers for Stage 1 location", q, ":\n")
mcmc_conf$printSamplers()

############################################################
#  Build and compile MCMC
############################################################

nimble_mcmc <- buildMCMC(mcmc_conf)

compiled_mcmc <- compileNimble(
  nimble_mcmc,
  project = nimble_model,
  resetFunctions = TRUE
)

st2 <- Sys.time()

############################################################
# Run MCMC
############################################################

mcmc_out <- runMCMC(
  compiled_mcmc,
  inits = mod_inits,
  nchains = 1,
  nburnin = M.burn,
  niter = M.iter,
  thin = M.thin,
  samplesAsCodaMCMC = TRUE,
  summary = TRUE,
  WAIC = FALSE,
  progressBar = TRUE
)

samples <- mcmc_out$samples
mcmc_summary <- mcmc_out$summary

############################################################
#  Extract posterior samples
############################################################

samples_mat <- as.matrix(samples)

beta_cols <- grep("^beta\\[", colnames(samples_mat))

gamma_rho_col <- which(colnames(samples_mat) == "gamma_rho")

rho_col <- which(colnames(samples_mat) == "rho")

alpha_col <- which(colnames(samples_mat) == "alpha")

ess_all <- effectiveSize(samples)

accept_rates_all <- apply(
  samples_mat,
  2,
  function(x) mean(diff(x) != 0)
)

st3 <- Sys.time()

############################################################
# Save stage one outputs
############################################################

MCMCout <- list(
  beta = samples_mat[, beta_cols, drop = FALSE],
  gamma_rho = samples_mat[, gamma_rho_col],
  rho = samples_mat[, rho_col],
  alpha = samples_mat[, alpha_col],
  ess_all = ess_all,
  accept_rates = accept_rates_all,
  mcmc_summary = mcmc_summary,
  time.out = st3 - st1,
  time.out.run = st3 - st2,
  q = q,
  grid = gridID[q],
  coords = coords[q, , drop = FALSE],
  vars = vars,
  bp = bp,
  Tobs = Tobs
)

if(!dir.exists("stage1_outputs")){
  dir.create("stage1_outputs")
}

save(
  MCMCout,
  file = paste0(
    "stage1_outputs/stage1_continuous_location_",
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




