
# Stage Two Count Model 

code_start <- Sys.time()

sink(file = "runtime.txt", append = TRUE)
cat(paste("code started at: ", Sys.time(), "\n"))
cat("\n")
sink()

library(invgamma)

set.seed(1234)

st1 <- Sys.time()

# Load data

data <- read.csv("county_data.csv")

rm('county_data')

countyID = unique(data$GEOID)

n = length(countyID) # n = 227
T = nrow(data)/n.    # T = 1174

tab <- table(data$GEOID)

if(length(unique(tab)) != 1){
  stop("Counties do not all have the same number of time points.")
}

yfull <- data$cases_c
pop.full <- data$pop_c

if(anyNA(yfull)){
  stop("Missing values detected in the full response.")
}

if(anyNA(pop.full) || any(pop.full <= 0)){
  stop("Population values must be positive and nonmissing.")
}

pooled_log_rate <- log(
  (sum(yfull) + 0.5) /
    sum(pop.full)
)

# Load the stage one samples

load(file.path(
  "stage1_outputs_count",
  paste0("stage1_county_", 1, ".Rda")
))

M <- nrow(MCMCout$beta)
bp <- ncol(MCMCout$beta)

# Storage for stage one posterior samples
pi.samp    <- matrix(NA, n, M)
size.samp  <- matrix(NA, n, M)
gamma.samp <- matrix(NA, n, M)
beta.samp  <- array(NA, c(n, bp, M))

rm('MCMCout')

# Consolidate stage one outputs

for(q in seq_len(n)){
  
  load(file.path(
    "stage1_outputs_count",
    paste0("stage1_county_", q, ".Rda")
  ))
  
  if(nrow(MCMCout$beta) != M){
    stop("Stage one sample size differs for county ", q)
  }
  
  if(ncol(MCMCout$beta) != bp){
    stop("Number of beta parameters differs for county ", q)
  }
  
  if(length(MCMCout$pi) != M ||
     length(MCMCout$size) != M ||
     length(MCMCout$gamma) != M){
    stop("Stage one sample lengths are inconsistent for county ", q)
  }
  
  pi.samp[q, ]       <- MCMCout$pi
  size.samp[q, ]     <- MCMCout$size
  gamma.samp[q, ]    <- MCMCout$gamma
  beta.samp[q, , ]   <- t(MCMCout$beta)
  
  rm(MCMCout)
}

# Adjacency matrix

# A = adjacency matrix
A <- A.county

if(!all(dim(A) == c(n, n))){
  stop("A.county dimension does not match number of counties.")
}

if(any(diag(A) != 0)){
  stop("Adjacency matrix should have zero diagonal.")
}

numnns <- rowSums(A) # number of neighbors for each site

if(any(numnns == 0)){
  stop("At least one county has no neighbors. ICAR conditional update will fail.")
}

D <- diag(numnns)
Qmat <- D - A
rank_Q <- qr(Qmat)$rank

##########################################################################
# Set up MCMC for stage two
##########################################################################

# initial values
# take one draw for each county-specific parameter

pi    <- pi.samp[, M]
size  <- size.samp[, M]
beta  <- beta.samp[, , M]
gamma <- gamma.samp[, M]

tausq.b <- rep(1, bp)
tausq.g <- 1

# MCMC settings

M.iter <- 200000
M.burn <- 50000
M.thin <- 50
M.out  <- (M.iter - M.burn) / M.thin

# Storage

pi.out      <- matrix(NA, n, M.out)
size.out    <- matrix(NA, n, M.out)
beta.out    <- array(NA, c(n, bp, M.out))
rho.out     <- matrix(NA, n, M.out)
gamma.out <- matrix(NA, n, M.out)

tausq.b.out <- matrix(NA, bp, M.out)
tausq.g.out <- rep(NA, M.out)

accept.site <- rep(0, n)
accept.site.post <- rep(0, n)

############################################################
# Stage two MCMC
############################################################

st2 <- Sys.time()

progress_bar <- txtProgressBar(min = 0, max = M.iter, style = 3, char = "=")

for(m in 1:M.iter){

  ##########################################################
  # Update ICAR variance (tausq.b) for beta
  ##########################################################
  
  for(p in 1:bp){
    
    quad.b <- as.numeric(t(beta[, p]) %*% Qmat %*% beta[, p])
    
    tausq.b[p] <- invgamma::rinvgamma(
      n = 1,
      shape = 0.5 + rank_Q / 2,
      rate = 0.5 + 0.5 * quad.b
    )
  }
  
  ##########################################################
  # Update ICAR variance (tausq.g) for gamma
  ##########################################################
  
  quad.g <- as.numeric(t(gamma) %*% Qmat %*% gamma)
  
  tausq.g <- invgamma::rinvgamma(
    n = 1,
    shape = 0.5 + rank_Q / 2,
    rate = 0.5 + 0.5 * quad.g
  )
  
  ##########################################################
  # Update county-specific parameters using stage one draws
  ##########################################################
  
  for(i in 1:n){
    
    mi <- sample.int(M, 1)
    
    beta.new  <- beta.samp[i, , mi]
    gamma.new <- gamma.samp[i, mi]
    pi.new    <- pi.samp[i, mi]
    size.new  <- size.samp[i, mi]
    
    gamma.mean <- as.numeric(A[i, ] %*% gamma) / numnns[i]
    beta.mean  <- as.numeric(A[i, ] %*% beta) / numnns[i]
  
    log.full.new <-
      dnorm(
        gamma.new,
        mean = gamma.mean,
        sd = sqrt(tausq.g / numnns[i]),
        log = TRUE
      ) +
      sum(dnorm(
        beta.new,
        mean = beta.mean,
        sd = sqrt(tausq.b / numnns[i]),
        log = TRUE
      ))
    
    log.full.old <-
      dnorm(
        gamma[i],
        mean = gamma.mean,
        sd = sqrt(tausq.g / numnns[i]),
        log = TRUE
      ) +
      sum(dnorm(
        beta[i, ],
        mean = beta.mean,
        sd = sqrt(tausq.b / numnns[i]),
        log = TRUE
      ))
    
    log.stage1.new <-
      dlogis(
        gamma.new,
        location = 0,
        scale = 0.5,
        log = TRUE
      ) +
      dnorm(
        beta.new[1],
        mean = pooled_log_rate,
        sd = 5,
        log = TRUE
      ) +
      sum(
        dnorm(
          beta.new[2:bp],
          mean = 0,
          sd = 2,
          log = TRUE
        )
      )
    
    log.stage1.old <-
      dlogis(
        gamma[i],
        location = 0,
        scale = 0.5,
        log = TRUE
      ) +
      dnorm(
        beta[i, 1],
        mean = pooled_log_rate,
        sd = 5,
        log = TRUE
      ) +
      sum(
        dnorm(
          beta[i, 2:bp],
          mean = 0,
          sd = 2,
          log = TRUE
        )
      )
    
    logR <- (log.full.new - log.stage1.new) -
      (log.full.old - log.stage1.old)
    
    if(is.na(logR)){
      stop(
        "NA logR at iteration ",
        m,
        ", county ",
        i
      )
    }
    
    if(log(runif(1)) < min(0, logR)){
      
      beta[i, ] <- beta.new
      gamma[i]  <- gamma.new
      pi[i]     <- pi.new
      size[i]   <- size.new
      
      accept.site[i] <- accept.site[i] + 1
      
      if(m > M.burn){
        accept.site.post[i] <- accept.site.post[i] + 1
      }
    }
  }
    
  ##########################################################
  # Save posterior draws
  ##########################################################
  
  if(m > M.burn && m/M.thin==floor(m/M.thin)){
    
    s <- (m - M.burn) / M.thin
    
    pi.out[, s]      <- pi
    size.out[, s]    <- size
    beta.out[, , s]  <- beta
    gamma.out[, s] <- gamma
    rho.out[, s]   <- tanh(gamma)
    
    tausq.b.out[, s] <- tausq.b
    tausq.g.out[s]   <- tausq.g
  }
  
  setTxtProgressBar(progress_bar, value = m)
}

close(progress_bar)
  
accept.site.rate <- accept.site / M.iter
accept.site.post.rate <- accept.site.post / (M.iter - M.burn)

cat("\nSite acceptance rates, all iterations:\n")
print(summary(accept.site.rate))

cat("\nSite acceptance rates, post burn-in:\n")
print(summary(accept.site.post.rate))

dim(pi.out)
dim(size.out)

dim(beta.out)
dim(tausq.b.out)

dim(rho.out)
dim(tausq.g.out)
length(tausq.g.out)

st3 <- Sys.time()

time.out.total <- st3 - st1
time.out.run <- st3 - st2

Stage2out <- list(
  pi = pi.out,
  size = size.out,
  beta = beta.out,
  rho = rho.out,
  gamma = gamma.out,
  tausq.b = tausq.b.out,
  tausq.g = tausq.g.out,
  
  accept.site = accept.site,
  accept.site.post = accept.site.post,
  accept.site.rate = accept.site.rate,
  accept.site.post.rate = accept.site.post.rate,
  accept.overall.rate = mean(accept.site.rate),
  accept.overall.post.rate = mean(accept.site.post.rate),
  
  time.out.total = time.out.total,
  time.out.run = time.out.run,
  
  M.iter = M.iter,
  M.burn = M.burn,
  M.thin = M.thin,
  n = n,
  bp = bp
)

if(!dir.exists("stage2_output")){
  dir.create("stage2_output")
}

save(
  Stage2out,
  file = "stage2_output/stage2_output_count.Rda"
)

#  Runtime
sink(file = "runtime.txt", append = TRUE)
cat(paste("code ended at: ", Sys.time(), "\n"))
print(difftime(Sys.time(), code_start))
cat("\n")
sink()



