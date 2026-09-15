
# Single-Stage Continuous Model

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

gridID <- training$grid[which(training$timeID == 1)]

# Number of locations
Q <- length(gridID)

# Coordinates for the Q locations
coords <- training[which(training$timeID == 1), c("lon", "lat")]

# Number of time points
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

# Design matrix with intercept
X <- cbind(rep(1, Q * Tobs), Xfull)

# Number of regression coefficients including intercept
bp <- ncol(X)

############################################################
# Convert design matrix and response into arrays
############################################################

# Xarr[i, t, p] = pth covariate for location i at time t
Xarr <- array(NA, dim = c(Q, Tobs, bp))

for(t in 1:Tobs){
  rows_t <- ((t - 1) * Q + 1):(t * Q)
  Xarr[, t, ] <- X[rows_t, ]
}

# Response matrix:
# Y[i, t] = response at location i and time t
Y <- matrix(yfull, nrow = Q, ncol = Tobs)

# logY is used in the centered AR(1) term
logY <- log(Y)

rm(data, Xcov, training, scalingvalues)

############################################################
# Distance matrix and phi bounds
############################################################

Dmat <- as.matrix(dist(coords))

nonzeroD <- Dmat[Dmat > 0]
minD <- min(nonzeroD)
maxD <- max(Dmat)

# phi in exponential correlation matrix.
phi_lower <- 0.01 / maxD
phi_upper <- 10 / minD

############################################################
# NIMBLE model
############################################################

nimble_code <- nimbleCode({
  
  ##########################################################
  #  Spatial hyperpriors
  ##########################################################
  
  # For betas
  for(p in 1:bp){
    mu_beta[p] ~ dnorm(0, sd = 10)
    sigma2_beta[p] ~ dinvgamma(2, 1)
    phi_beta[p] ~ dunif(phi_lower, phi_upper)
  }
  
  # For transformed rho 
  mu_gamma_rho ~ dnorm(0, sd = 10)
  sigma2_gamma_rho ~ dinvgamma(2, 1)
  phi_gamma_rho ~ dunif(phi_lower, phi_upper)
  
  # Gamma shape parameters
  for(i in 1:Q){
    alpha[i] ~ dgamma(0.01, 0.01)
  }
  
  ##########################################################
  # Means and Covariances
  ##########################################################
  
  for(i in 1:Q){
    
    one[i] <- 1
    
    mean_gamma_rho[i] <- mu_gamma_rho * one[i]
    
    for(p in 1:bp){
      mean_beta[i, p] <- mu_beta[p] * one[i]
    }
    
    for(j in 1:Q){
      
      Sigma_gamma_rho[i, j] <-
        sigma2_gamma_rho * exp(-phi_gamma_rho * Dmat[i, j])
      
      for(p in 1:bp){
        Sigma_beta[i, j, p] <-
          sigma2_beta[p] * exp(-phi_beta[p] * Dmat[i, j])
      }
    }
  }
  
  ##########################################################
  # Spatial MVN Priors
  ##########################################################
  
  for(p in 1:bp){
    beta[1:Q, p] ~ dmnorm(
      mean = mean_beta[1:Q, p],
      cov = Sigma_beta[1:Q, 1:Q, p]
    )
  }
  
  gamma_rho[1:Q] ~ dmnorm(
    mean = mean_gamma_rho[1:Q],
    cov = Sigma_gamma_rho[1:Q, 1:Q]
  )
  
  ##########################################################
  # Data Model
  ##########################################################
  
  for(i in 1:Q){
    
    rho[i] <- 2 / (1 + exp(-2 * gamma_rho[i])) - 1
    
    ########################################################
    # Time t = 1
    ########################################################
    
    xb[i, 1] <- inprod(Xarr[i, 1, 1:bp], beta[i, 1:bp])
    
    eta_mean[i, 1] <- xb[i, 1]
    
    mu[i, 1] <- exp(eta_mean[i, 1])
    
    rateY[i, 1] <- alpha[i] / mu[i, 1]
    
    Y[i, 1] ~ dgamma(shape = alpha[i], rate = rateY[i, 1])
    
    ########################################################
    # Time t = 2,...,Tobs
    ########################################################
    
    for(t in 2:Tobs){
      
      xb[i, t] <- inprod(Xarr[i, t, 1:bp], beta[i, 1:bp])
      
      resid_prev[i, t] <- logY[i, t - 1] - xb[i, t - 1]
      
      eta_mean[i, t] <- xb[i, t] + rho[i] * resid_prev[i, t]
      
      mu[i, t] <- exp(eta_mean[i, t])
      
      rateY[i, t] <- alpha[i] / mu[i, t]
      
      Y[i, t] ~ dgamma(shape = alpha[i], rate = rateY[i, t])
    }
  }
})

############################################################
# Constants and data
############################################################

constants <- list(
  Q = Q,
  Tobs = Tobs,
  bp = bp,
  Xarr = Xarr,
  logY = logY,
  Dmat = Dmat,
  phi_lower = phi_lower,
  phi_upper = phi_upper
)

data_list <- list(
  Y = Y
)

############################################################
# Initial values
############################################################

phi_init <- 3 / maxD

if(phi_init <= phi_lower || phi_init >= phi_upper){
  phi_init <- 0.25 * phi_upper
}

# Initial beta
beta_init <- matrix(
  0,
  nrow = Q,
  ncol = bp
)

# Initial transformed rho
gamma_rho_init <- rep(0, Q)

inits <- list(
  
  # beta
  beta = beta_init,
  
  # transformed rho
  gamma_rho = gamma_rho_init,
  
  # Gamma shape parameters
  alpha = rep(2, Q),
  
  # beta hyperparameters
  mu_beta = rep(0, bp),
  sigma2_beta = rep(1, bp),
  phi_beta = rep(phi_init, bp),
  
  # gamma_rho hyperparameters
  mu_gamma_rho = 0,
  sigma2_gamma_rho = 1,
  phi_gamma_rho = phi_init
)

############################################################
# Build and compile model
############################################################

Rmodel <- nimbleModel(
  code = nimble_code,
  constants = constants,
  data = data_list,
  inits = inits
)

Cmodel <- compileNimble(Rmodel)

############################################################
# Configure MCMC
############################################################

monitors <- c(
  "beta",
  "rho",
  "gamma_rho",
  "alpha",
  "mu_beta",
  "sigma2_beta",
  "phi_beta",
  "mu_gamma_rho",
  "sigma2_gamma_rho",
  "phi_gamma_rho"
)

conf <- configureMCMC(
  Rmodel,
  monitors = monitors,
  useConjugacy = TRUE
)

# Remove the default RW_block samplers for the spatial fields
conf$removeSamplers(c("beta", "gamma_rho"))

# Add AF_slice samplers
for(p in 1:bp) {
  conf$addSampler(target = paste0("beta[1:", Q, ", ", p, "]"), type = "AF_slice")
}
conf$addSampler(target = paste0("gamma_rho[1:", Q, "]"), type = "AF_slice")

cat("\nAll samplers:\n")
conf$printSamplers()

Rmcmc <- buildMCMC(conf)

Cmcmc <- compileNimble(
  Rmcmc,
  project = Rmodel
)

st2 <- Sys.time()

############################################################
# Run MCMC
############################################################

M.iter <- 10000
M.burn <- 4000
M.thin <- 2

mcmc_out <- runMCMC(
  Cmcmc,
  niter = M.iter,
  nburnin = M.burn,
  thin = M.thin,
  nchains = 1,
  samplesAsCodaMCMC = TRUE,
  summary = TRUE,
  setSeed = TRUE
)

samples <- mcmc_out$samples
mcmc_summary <- mcmc_out$summary

samples_mat <- as.matrix(samples)

############################################################
# Extract posterior columns
############################################################

beta_cols <- grep("^beta\\[", colnames(samples_mat))

rho_cols <- grep("^rho\\[", colnames(samples_mat))

gamma_rho_cols <- grep("^gamma_rho\\[", colnames(samples_mat))

alpha_cols <- grep("^alpha\\[", colnames(samples_mat))

mu_beta_cols <- grep("^mu_beta\\[", colnames(samples_mat))

sigma2_beta_cols <- grep("^sigma2_beta\\[", colnames(samples_mat))

phi_beta_cols <- grep("^phi_beta\\[", colnames(samples_mat))

mu_gamma_rho_col <- which(colnames(samples_mat) == "mu_gamma_rho")

sigma2_gamma_rho_col <- which(colnames(samples_mat) == "sigma2_gamma_rho")

phi_gamma_rho_col <- which(colnames(samples_mat) == "phi_gamma_rho")

############################################################
# Diagnostics
############################################################

ess_all <- effectiveSize(samples)

global_names <- c(
  colnames(samples_mat)[mu_beta_cols],
  colnames(samples_mat)[sigma2_beta_cols],
  colnames(samples_mat)[phi_beta_cols],
  "mu_gamma_rho",
  "sigma2_gamma_rho",
  "phi_gamma_rho"
)

ess_global <- ess_all[global_names]

accept_rates_all <- apply(
  samples_mat,
  2,
  function(x) mean(diff(x) != 0)
)

st3 <- Sys.time()

############################################################
# Save output
############################################################

SingleOut <- list(
  
  # Location-specific posterior samples
  beta = samples_mat[, beta_cols, drop = FALSE],
  rho = samples_mat[, rho_cols, drop = FALSE],
  gamma_rho = samples_mat[, gamma_rho_cols, drop = FALSE],
  alpha = samples_mat[, alpha_cols, drop = FALSE],
  
  # Spatial hyperparameter posterior samples for beta
  mu_beta = samples_mat[, mu_beta_cols, drop = FALSE],
  sigma2_beta = samples_mat[, sigma2_beta_cols, drop = FALSE],
  phi_beta = samples_mat[, phi_beta_cols, drop = FALSE],
  
  # Spatial hyperparameter posterior samples for gamma_rho
  mu_gamma_rho = samples_mat[, mu_gamma_rho_col],
  sigma2_gamma_rho = samples_mat[, sigma2_gamma_rho_col],
  phi_gamma_rho = samples_mat[, phi_gamma_rho_col],
  
  # Diagnostics
  ess_all = ess_all,
  ess_global = ess_global,
  accept_rates = accept_rates_all,
  mcmc_summary = mcmc_summary,
  
  # metadata
  vars = vars,
  bp = bp,
  Q = Q,
  Tobs = Tobs,
  phi_lower = phi_lower,
  phi_upper = phi_upper,
  
  # Timing
  time.out = st3 - st1,
  time.out.run = st3 - st2
)

if(!dir.exists("single_output")){
  dir.create("single_output")
}

save(
  SingleOut,
  file = "single_output/single_continuous_output.Rda"
)

# Runtime
sink(file = "runtime.txt", append = TRUE)
cat(paste("code ended at: ", Sys.time(), "\n"))
print(difftime(Sys.time(), code_start))
cat("\n")
sink()




