
# Stage Two Continuous Model

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
library(mvtnorm)

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

yfull <- training[, "soilm"]

if(any(yfull <= 0, na.rm = TRUE)){
  stop("Gamma response soilm must be strictly positive.")
}

rm(data, Xcov, training, scalingvalues)

############################################################
#  Load stage one posterior samples
############################################################

stage1_files <- paste0(
  "stage1_outputs/stage1_continuous_location_",
  1:Q,
  ".Rda"
)

load_stage1_MCMCout <- function(file_name){
  
  tmp_env <- new.env()
  
  load(file_name, envir = tmp_env)
  
  if(!exists("MCMCout", envir = tmp_env)){
    stop(paste("File does not contain object named MCMCout:", file_name))
  }
  
  return(tmp_env$MCMCout)
}

first_out <- load_stage1_MCMCout(stage1_files[1])

M <- nrow(first_out$beta)

bp <- ncol(first_out$beta)

beta.samp <- array(NA, dim = c(Q, bp, M))

gamma_rho.samp <- matrix(NA, nrow = Q, ncol = M)

rho.samp <- matrix(NA, nrow = Q, ncol = M)

alpha.samp <- matrix(NA, nrow = Q, ncol = M)

rm(first_out)

for(i in 1:Q){
  
  MCMCout_i <- load_stage1_MCMCout(stage1_files[i])
  
  beta.samp[i, , ] <- t(MCMCout_i$beta)
  
  gamma_rho.samp[i, ] <- MCMCout_i$gamma_rho
  
  rho.samp[i, ] <- MCMCout_i$rho
  
  alpha.samp[i, ] <- MCMCout_i$alpha
  
  rm(MCMCout_i)
}

############################################################
#  Distance matrix and phi bounds
############################################################

Dmat <- as.matrix(dist(coords))

nonzeroD <- Dmat[Dmat > 0]

minD <- min(nonzeroD)

maxD <- max(Dmat)

phi_lower <- 0.01 / maxD

phi_upper <- 10 / minD

############################################################
#  Hyperprior settings
############################################################

# Spatial variance priors:
# sigma2_beta[p] ~ IG(a_sigma2_beta, b_sigma2_beta)
# sigma2_gamma_rho ~ IG(a_sigma2_gamma_rho, b_sigma2_gamma_rho)

a_sigma2_beta <- 2
b_sigma2_beta <- 1

a_sigma2_gamma_rho <- 2
b_sigma2_gamma_rho <- 1

# Global mean priors:
# mu_beta[p] ~ N(0, sd_mu_beta^2)
# mu_gamma_rho ~ N(0, sd_mu_gamma_rho^2)

sd_mu_beta <- 10
sd_mu_gamma_rho <- 10

############################################################
#  MCMC settings
############################################################

M.iter <- 200000
M.burn <- 50000
M.thin <- 50

M.out <- (M.iter - M.burn) / M.thin

############################################################
#  Initial values
############################################################

# Initialize location-specific parameters from last stage one draw
beta <- beta.samp[, , M]

gamma_rho <- gamma_rho.samp[, M]

rho <- rho.samp[, M]

alpha <- alpha.samp[, M]

# Initialize spatial hyperparameters
mu_beta <- apply(beta, 2, mean)

sigma2_beta <- apply(beta, 2, var)

sigma2_beta[!is.finite(sigma2_beta) | sigma2_beta <= 0] <- 1

phi_init <- 3 / maxD

if(phi_init <= phi_lower || phi_init >= phi_upper){
  phi_init <- 0.25 * phi_upper
}

phi_beta <- rep(phi_init, bp)

mu_gamma_rho <- mean(gamma_rho)

sigma2_gamma_rho <- var(gamma_rho)

if(!is.finite(sigma2_gamma_rho) || sigma2_gamma_rho <= 0){
  sigma2_gamma_rho <- 1
}

phi_gamma_rho <- phi_init

############################################################
# Storage
############################################################

beta.out <- array(NA, dim = c(Q, bp, M.out))

gamma_rho.out <- matrix(NA, nrow = Q, ncol = M.out)

rho.out <- matrix(NA, nrow = Q, ncol = M.out)

alpha.out <- matrix(NA, nrow = Q, ncol = M.out)

mu_beta.out <- matrix(NA, nrow = bp, ncol = M.out)

sigma2_beta.out <- matrix(NA, nrow = bp, ncol = M.out)

phi_beta.out <- matrix(NA, nrow = bp, ncol = M.out)

mu_gamma_rho.out <- rep(NA, M.out)

sigma2_gamma_rho.out <- rep(NA, M.out)

phi_gamma_rho.out <- rep(NA, M.out)

accept.site <- rep(0, Q)

accept.phi.beta <- rep(0, bp)

accept.phi.gamma_rho <- 0

accept.phi.beta.post <- rep(0, bp)

accept.phi.gamma_rho.post <- 0

n.phi.post <- 0


############################################################
# functions for M-H update for location specific parameters
############################################################

make_spatial_K_ld <- function(phi, Dmat, phi_lower, phi_upper){
  
  if(phi <= phi_lower || phi >= phi_upper){
    stop("phi is outside the allowed interval.")
  }
  
  R <- exp(-phi * Dmat)
  
  cholR <- chol(R)
  
  K <- chol2inv(cholR)
  
  ld <- 2 * sum(log(diag(cholR)))
  
  list(
    K = K,
    ld = ld
  )
}

log_mvn_spatial_fast <- function(x, mu_scalar, sigma2, K, ld){
  
  if(sigma2 <= 0){
    return(-Inf)
  }
  
  Q <- length(x)
  
  z <- x - mu_scalar
  
  Kz <- as.numeric(K %*% z)
  
  qf <- sum(z * Kz)
  
  out <-
    -0.5 * Q * log(2 * pi) -
    0.5 * Q * log(sigma2) -
    0.5 * ld -
    0.5 * qf / sigma2
  
  return(out)
}

log_beta_spatial_all_fast <- function(beta, mu_beta, sigma2_beta,
                                      K_beta_list, ld_beta_vec){
  
  bp <- ncol(beta)
  
  out <- 0
  
  for(p in 1:bp){
    
    out <- out +
      log_mvn_spatial_fast(
        x = beta[, p],
        mu_scalar = mu_beta[p],
        sigma2 = sigma2_beta[p],
        K = K_beta_list[[p]],
        ld = ld_beta_vec[p]
      )
  }
  
  return(out)
}

log_gamma_rho_spatial_fast <- function(gamma_rho, mu_gamma_rho,
                                       sigma2_gamma_rho,
                                       K_gamma, ld_gamma){
  
  log_mvn_spatial_fast(
    x = gamma_rho,
    mu_scalar = mu_gamma_rho,
    sigma2 = sigma2_gamma_rho,
    K = K_gamma,
    ld = ld_gamma
  )
}


# Stage one independent priors used in the proposal correction.
# Stage one priors:
#   beta[p] ~ N(0, 10^2)
#   gamma_rho ~ N(0, 10^2)
#
# alpha has the same prior in stage one and stage two, so it cancels.

log_stage1_param_prior <- function(beta_i, gamma_rho_i){
  
  sum(dnorm(beta_i, mean = 0, sd = 10, log = TRUE)) +
    dnorm(gamma_rho_i, mean = 0, sd = 10, log = TRUE)
}

############################################################
# Local spatial-prior update functions
############################################################

make_spatial_state <- function(x, mu_scalar, sigma2, K){
  
  z <- x - mu_scalar
  
  Kz <- as.numeric(K %*% z)
  
  list(
    z = z,
    Kz = Kz,
    Kdiag = diag(K),
    sigma2 = sigma2
  )
}

local_logprior_change <- function(delta, i, state){
  
  delta_qf <-
    2 * delta * state$Kz[i] +
    delta^2 * state$Kdiag[i]
  
  -0.5 * delta_qf / state$sigma2
}

update_spatial_state_one <- function(state, K, delta, i){
  
  if(delta == 0){
    return(state)
  }
  
  state$z[i] <- state$z[i] + delta
  
  state$Kz <- state$Kz + delta * K[, i]
  
  state
}


############################################################
# functions for phi transformation and M-H update for phi
############################################################

phi_to_z <- function(phi, phi_lower, phi_upper){
  log((phi - phi_lower) / (phi_upper - phi))
}

### plogis(z) = exp(z) / (1 + exp(z))
z_to_phi <- function(z, phi_lower, phi_upper){
  phi_lower + (phi_upper - phi_lower) * plogis(z)
}

log_phi_jacobian <- function(phi, phi_lower, phi_upper){
  
  if(phi <= phi_lower || phi >= phi_upper){
    return(-Inf)
  }
  
  log(phi - phi_lower) +
    log(phi_upper - phi) -
    log(phi_upper - phi_lower)
}


############################################################
# M-H update for phi
############################################################

update_phi_MH_z_fast <- function(phi_current, K_current, ld_current,
                                 x, mu_scalar, sigma2, Dmat,
                                 phi_lower, phi_upper,
                                 proposal_sd = 0.5){
  
  z_old <- phi_to_z(
    phi = phi_current,
    phi_lower = phi_lower,
    phi_upper = phi_upper
  )
  
  z_new <- rnorm(
    n = 1,
    mean = z_old,
    sd = proposal_sd
  )
  
  phi_new <- z_to_phi(
    z = z_new,
    phi_lower = phi_lower,
    phi_upper = phi_upper
  )
  
  ##########################################################
  # Old target
  ##########################################################
  
  log_old <-
    log_mvn_spatial_fast(
      x = x,
      mu_scalar = mu_scalar,
      sigma2 = sigma2,
      K = K_current,
      ld = ld_current
    ) +
    log_phi_jacobian(
      phi = phi_current,
      phi_lower = phi_lower,
      phi_upper = phi_upper
    )
  
  ##########################################################
  # New target: phi changed, so we must compute new K and ld
  ##########################################################
  
  new_cache <- make_spatial_K_ld(
    phi = phi_new,
    Dmat = Dmat,
    phi_lower = phi_lower,
    phi_upper = phi_upper
  )
  
  log_new <-
    log_mvn_spatial_fast(
      x = x,
      mu_scalar = mu_scalar,
      sigma2 = sigma2,
      K = new_cache$K,
      ld = new_cache$ld
    ) +
    log_phi_jacobian(
      phi = phi_new,
      phi_lower = phi_lower,
      phi_upper = phi_upper
    )
  
  log_accept <- log_new - log_old
  
  if(log(runif(1)) < log_accept){
    
    return(
      list(
        phi = phi_new,
        K = new_cache$K,
        ld = new_cache$ld,
        accepted = 1
      )
    )
    
  } else {
    
    return(
      list(
        phi = phi_current,
        K = K_current,
        ld = ld_current,
        accepted = 0
      )
    )
  }
}


############################################################
# Gibbs updates for mu and sigma2
############################################################

update_mu_spatial_fast <- function(x, sigma2, K, sd_mu = 10){
  
  Q <- length(x)
  
  one_vec <- rep(1, Q)
  
  K_one <- as.numeric(K %*% one_vec)
  
  precision_mu <-
    1 / sd_mu^2 +
    sum(one_vec * K_one) / sigma2
  
  variance_mu <- 1 / precision_mu
  
  mean_mu <-
    variance_mu *
    sum(x * K_one) / sigma2
  
  rnorm(1, mean = mean_mu, sd = sqrt(variance_mu))
}

update_sigma2_spatial_fast <- function(x, mu_scalar, K,
                                       a_sigma2, b_sigma2){
  
  Q <- length(x)
  
  z <- x - mu_scalar
  
  Kz <- as.numeric(K %*% z)
  
  qf <- sum(z * Kz)
  
  shape <- a_sigma2 + Q / 2
  
  rate <- b_sigma2 + 0.5 * qf
  
  1 / rgamma(1, shape = shape, rate = rate)
}


############################################################
# Adaptive tuning settings for phi
############################################################

phi_beta_proposal_sd <- rep(0.5, bp)

phi_gamma_rho_proposal_sd <- 0.5

block_size <- 50

accept_phi_beta_block <- rep(0, bp)

accept_phi_gamma_rho_block <- 0

############################################################
# Initialize spatial inverse matrices and log determinants
############################################################

K_beta_list <- vector("list", bp)

ld_beta_vec <- numeric(bp)

for(p in 1:bp){
  
  tmp_cache <- make_spatial_K_ld(
    phi = phi_beta[p],
    Dmat = Dmat,
    phi_lower = phi_lower,
    phi_upper = phi_upper
  )
  
  K_beta_list[[p]] <- tmp_cache$K
  
  ld_beta_vec[p] <- tmp_cache$ld
}

tmp_gamma_cache <- make_spatial_K_ld(
  phi = phi_gamma_rho,
  Dmat = Dmat,
  phi_lower = phi_lower,
  phi_upper = phi_upper
)

K_gamma <- tmp_gamma_cache$K

ld_gamma <- tmp_gamma_cache$ld


############################################################
# Stage Two MCMC
############################################################

st2 <- Sys.time()

progress_bar <- txtProgressBar(
  min = 0,
  max = M.iter,
  style = 3,
  char = "="
)

for(m in 1:M.iter){
  
  ##########################################################
  # Update beta spatial hyperparameters
  ##########################################################
  
  for(p in 1:bp){
    
    ########################################################
    # Update mu_beta[p] using current K
    ########################################################
    
    mu_beta[p] <- update_mu_spatial_fast(
      x = beta[, p],
      sigma2 = sigma2_beta[p],
      K = K_beta_list[[p]],
      sd_mu = sd_mu_beta
    )
    
    ########################################################
    # Update sigma2_beta[p] using current K
    ########################################################
    
    sigma2_beta[p] <- update_sigma2_spatial_fast(
      x = beta[, p],
      mu_scalar = mu_beta[p],
      K = K_beta_list[[p]],
      a_sigma2 = a_sigma2_beta,
      b_sigma2 = b_sigma2_beta
    )
    
    ########################################################
    # Update phi_beta[p]; this also updates K and ld
    ########################################################
    
    phi_update <- update_phi_MH_z_fast(
      phi_current = phi_beta[p],
      K_current = K_beta_list[[p]],
      ld_current = ld_beta_vec[p],
      x = beta[, p],
      mu_scalar = mu_beta[p],
      sigma2 = sigma2_beta[p],
      Dmat = Dmat,
      phi_lower = phi_lower,
      phi_upper = phi_upper,
      proposal_sd = phi_beta_proposal_sd[p]
    )
    
    phi_beta[p] <- phi_update$phi
    
    K_beta_list[[p]] <- phi_update$K
    
    ld_beta_vec[p] <- phi_update$ld
    
    accept.phi.beta[p] <- accept.phi.beta[p] + phi_update$accepted
    
    accept_phi_beta_block[p] <-
      accept_phi_beta_block[p] + phi_update$accepted
    
    if(m > M.burn){
      accept.phi.beta.post[p] <-
        accept.phi.beta.post[p] + phi_update$accepted
    }
  }
  
  
  ##########################################################
  # Update gamma_rho spatial hyperparameters
  ##########################################################
  
  mu_gamma_rho <- update_mu_spatial_fast(
    x = gamma_rho,
    sigma2 = sigma2_gamma_rho,
    K = K_gamma,
    sd_mu = sd_mu_gamma_rho
  )
  
  sigma2_gamma_rho <- update_sigma2_spatial_fast(
    x = gamma_rho,
    mu_scalar = mu_gamma_rho,
    K = K_gamma,
    a_sigma2 = a_sigma2_gamma_rho,
    b_sigma2 = b_sigma2_gamma_rho
  )
  
  phi_update_gamma <- update_phi_MH_z_fast(
    phi_current = phi_gamma_rho,
    K_current = K_gamma,
    ld_current = ld_gamma,
    x = gamma_rho,
    mu_scalar = mu_gamma_rho,
    sigma2 = sigma2_gamma_rho,
    Dmat = Dmat,
    phi_lower = phi_lower,
    phi_upper = phi_upper,
    proposal_sd = phi_gamma_rho_proposal_sd
  )
  
  phi_gamma_rho <- phi_update_gamma$phi
  
  K_gamma <- phi_update_gamma$K
  
  ld_gamma <- phi_update_gamma$ld
  
  accept.phi.gamma_rho <-
    accept.phi.gamma_rho + phi_update_gamma$accepted
  
  accept_phi_gamma_rho_block <-
    accept_phi_gamma_rho_block + phi_update_gamma$accepted
  
  if(m > M.burn){
    
    accept.phi.gamma_rho.post <-
      accept.phi.gamma_rho.post + phi_update_gamma$accepted
    
    n.phi.post <- n.phi.post + 1
  }
  
  ##########################################################
  # Adaptive tuning for phi proposals during burn-in
  ##########################################################
  
  if(m <= M.burn && m %% block_size == 0){
    
    for(p in 1:bp){
      
      acc_rate_phi_beta_block <- accept_phi_beta_block[p] / block_size
      
      if(acc_rate_phi_beta_block > 0.50){
        phi_beta_proposal_sd[p] <- phi_beta_proposal_sd[p] * 1.5
      }
      
      if(acc_rate_phi_beta_block < 0.25){
        phi_beta_proposal_sd[p] <- phi_beta_proposal_sd[p] * 0.75
      }
      
      accept_phi_beta_block[p] <- 0
    }
    
    acc_rate_phi_gamma_block <-
      accept_phi_gamma_rho_block / block_size
    
    if(acc_rate_phi_gamma_block > 0.50){
      phi_gamma_rho_proposal_sd <- phi_gamma_rho_proposal_sd * 1.5
    }
    
    if(acc_rate_phi_gamma_block < 0.25){
      phi_gamma_rho_proposal_sd <- phi_gamma_rho_proposal_sd * 0.75
    }
    
    accept_phi_gamma_rho_block <- 0
  }
  
  
  ##########################################################
  # Build local spatial states before location updates
  # These are used to compute local spatial-prior changes.
  ##########################################################
  
  beta_state_list <- vector("list", bp)
  
  for(p in 1:bp){
    
    beta_state_list[[p]] <- make_spatial_state(
      x = beta[, p],
      mu_scalar = mu_beta[p],
      sigma2 = sigma2_beta[p],
      K = K_beta_list[[p]]
    )
  }
  
  gamma_state <- make_spatial_state(
    x = gamma_rho,
    mu_scalar = mu_gamma_rho,
    sigma2 = sigma2_gamma_rho,
    K = K_gamma
  )
  
  
  ##########################################################
  # Location-specific M-H updates using local spatial-prior
  # changes instead of full spatial log-prior
  ##########################################################
  
  for(i in 1:Q){
    
    mi <- sample.int(M, 1)
    
    beta.new.i <- beta.samp[i, , mi]
    gamma_rho.new.i <- gamma_rho.samp[i, mi]
    rho.new.i <- rho.samp[i, mi]
    alpha.new.i <- alpha.samp[i, mi]
    
    ########################################################
    # Compute local spatial-prior change
    ########################################################
    
    log_spatial_change <- 0
    
    for(p in 1:bp){
      
      delta_beta <- beta.new.i[p] - beta[i, p]
      
      log_spatial_change <- log_spatial_change +
        local_logprior_change(
          delta = delta_beta,
          i = i,
          state = beta_state_list[[p]]
        )
    }
    
    delta_gamma <- gamma_rho.new.i - gamma_rho[i]
    
    log_spatial_change <- log_spatial_change +
      local_logprior_change(
        delta = delta_gamma,
        i = i,
        state = gamma_state
      )
    
    ########################################################
    # Stage one prior
    ########################################################
    
    log_stage1_prior_old <- log_stage1_param_prior(
      beta_i = beta[i, ],
      gamma_rho_i = gamma_rho[i]
    )
    
    log_stage1_prior_new <- log_stage1_param_prior(
      beta_i = beta.new.i,
      gamma_rho_i = gamma_rho.new.i
    )
    
    ########################################################
    # M-H log acceptance ratio
    ########################################################
    
    log_accept <-
      log_spatial_change +
      log_stage1_prior_old -
      log_stage1_prior_new
    
    ########################################################
    # Accept/reject
    ########################################################
    
    if(log(runif(1)) < log_accept){
      
      ######################################################
      # Update beta states and beta values
      ######################################################
      
      for(p in 1:bp){
        
        delta_beta <- beta.new.i[p] - beta[i, p]
        
        beta_state_list[[p]] <- update_spatial_state_one(
          state = beta_state_list[[p]],
          K = K_beta_list[[p]],
          delta = delta_beta,
          i = i
        )
        
        beta[i, p] <- beta.new.i[p]
      }
      
      ######################################################
      # Update gamma state and gamma/rho/alpha values
      ######################################################
      
      gamma_state <- update_spatial_state_one(
        state = gamma_state,
        K = K_gamma,
        delta = delta_gamma,
        i = i
      )
      
      gamma_rho[i] <- gamma_rho.new.i
      rho[i] <- rho.new.i
      alpha[i] <- alpha.new.i
      
      accept.site[i] <- accept.site[i] + 1
    }
  }
  

  ##########################################################
  # Save posterior draws
  ##########################################################
  
  if(m > M.burn && m %% M.thin == 0){
    
    save_id <- (m - M.burn) / M.thin
    
    beta.out[, , save_id] <- beta
    
    gamma_rho.out[, save_id] <- gamma_rho
    
    rho.out[, save_id] <- rho
    
    alpha.out[, save_id] <- alpha
    
    mu_beta.out[, save_id] <- mu_beta
    
    sigma2_beta.out[, save_id] <- sigma2_beta
    
    phi_beta.out[, save_id] <- phi_beta
    
    mu_gamma_rho.out[save_id] <- mu_gamma_rho
    
    sigma2_gamma_rho.out[save_id] <- sigma2_gamma_rho
    
    phi_gamma_rho.out[save_id] <- phi_gamma_rho
  }
  
  if(m %% 100 == 0){
    setTxtProgressBar(progress_bar, value = m)
  }
}

setTxtProgressBar(progress_bar, value = M.iter)

close(progress_bar)

############################################################
# Acceptance rates
############################################################

accept.site.rate <- accept.site / M.iter

accept.phi.beta.rate <- accept.phi.beta / M.iter

accept.phi.gamma_rho.rate <- accept.phi.gamma_rho / M.iter

accept.phi.beta.post.rate <- accept.phi.beta.post / n.phi.post

accept.phi.gamma_rho.post.rate <-
  accept.phi.gamma_rho.post / n.phi.post


############################################################
# Final stage two output
############################################################

st3 <- Sys.time()

Stage2out <- list(
  
  beta = beta.out,
  gamma_rho = gamma_rho.out,
  rho = rho.out,
  alpha = alpha.out,
  
  mu_beta = mu_beta.out,
  sigma2_beta = sigma2_beta.out,
  phi_beta = phi_beta.out,
  
  mu_gamma_rho = mu_gamma_rho.out,
  sigma2_gamma_rho = sigma2_gamma_rho.out,
  phi_gamma_rho = phi_gamma_rho.out,
  
  phi_lower = phi_lower,
  phi_upper = phi_upper,
  
  phi_beta_proposal_sd_final = phi_beta_proposal_sd,
  phi_gamma_rho_proposal_sd_final = phi_gamma_rho_proposal_sd,
  
  accept.site.rate = accept.site.rate,
  accept.phi.beta.rate = accept.phi.beta.rate,
  accept.phi.gamma_rho.rate = accept.phi.gamma_rho.rate,
  accept.phi.beta.post.rate = accept.phi.beta.post.rate,
  accept.phi.gamma_rho.post.rate = accept.phi.gamma_rho.post.rate,
  
  vars = vars,
  bp = bp,
  Q = Q,
  Tobs = Tobs,
  
  a_sigma2_beta = a_sigma2_beta,
  b_sigma2_beta = b_sigma2_beta,
  a_sigma2_gamma_rho = a_sigma2_gamma_rho,
  b_sigma2_gamma_rho = b_sigma2_gamma_rho,
  sd_mu_beta = sd_mu_beta,
  sd_mu_gamma_rho = sd_mu_gamma_rho,
  
  time.out = st3 - st1,
  time.out.run = st3 - st2
)

if(!dir.exists("stage2_output")){
  dir.create("stage2_output")
}

save(
  Stage2out,
  file = "stage2_output/stage2_continuous_output.Rda"
)

# Runtime
sink(file = "runtime.txt", append = TRUE)
cat(paste("code ended at: ", Sys.time(), "\n"))
print(difftime(Sys.time(), code_start))
cat("\n")
sink()




