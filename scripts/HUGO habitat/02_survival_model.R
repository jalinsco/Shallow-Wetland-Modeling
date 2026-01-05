
#############################################X
#------------ Survival Analysis -------------X
#############################################X
# ------- fits CJS survival models  ---------X
#############################################X


# Load packages ----------------------

# required R packages
required_pkgs <- c("readr", 
                   "jagsUI",
                   "loo", 
                   "tibble",
                   "ggplot2",
                   "dplyr")

# load 
missing <- required_pkgs[!vapply(required_pkgs, 
                                 requireNamespace, 
                                 logical(1), 
                                 quietly = TRUE)]

# check if missing
if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "))
}





library(tidyverse)
library(data.table)
library(jagsUI)
library(loo)
library(parallel)



# Load data ----------------------

# load capture/re-sighting matrix
return_df <- read_csv("data/HUGO habitat/survival_df.csv")
return_df %>% as_tibble() %>% print(n = 500)

# separate column describing sex
sex <- return_df %>% mutate(sex = ifelse(Sex == "M", 1, 2)) %>% pull(sex)
return_df <- return_df %>% dplyr::select(-Sex)

# fill missing years
return_df <- return_df %>% add_column(`2018` = 0, .before = '2019')
return_df <- return_df %>% add_column(`2013` = 0, .before = '2014')
return_df

# matrix for 2009-2024
# 2009 = first year of captures
# 2024 = last year of re-sightings
Y <- as.matrix(return_df[,2:17]) 
head(Y)

# first encounter index
first <- apply(Y, 1, function(x) {
  w <- which(x == 1)
  if (length(w)) w[1] else NA_integer_
})
first

# effort vector (reflects whether resightings occurred)
miss <- rep(0, length = ncol(Y))
miss[5] <- 1
miss[10] <- 1

# count of individuals and years 
n_ind <- nrow(Y)
n_years <- ncol(Y)

# load drought indices (SPEI over 180 days preceding May 1)
climate_df <- read_csv("data/HUGO habitat/survival_climate_vals.csv")
drought = as.vector(climate_df$drought_idx) 

# check
length(drought)
n_years
stopifnot(length(drought) == n_years - 1)

# set up
jags_data <- list(
  Y = Y,
  drought = drought,
  n_years = n_years,
  n_ind = n_ind,
  first = first,
  miss = miss
)

# Prepare for model fitting ----------------------

# inits for z added to avoid impossible dead/alive states
z_inits <- matrix(NA_real_, nrow(Y), ncol(Y))
for (i in seq_len(nrow(Y))) {
  f <- first[i]
  if (is.na(f) || f >= ncol(Y)) next
  last1 <- if (any(Y[i,] == 1)) max(which(Y[i,] == 1)) else f
  if (last1 > f) z_inits[i, (f+1):last1] <- 1 # for alive btw first encounter and last 1 
  if (last1 < ncol(Y)) z_inits[i, (last1+1):ncol(Y)] <- 0 # for dead after last 1 
  
}
z_inits

# model runs
n_chains = 4
n_iter = 20000
n_burn = 4000
n_adapt = 2000

# number of cores?
parallel::detectCores()


# Fit Models ----------------------------------------

# M1. SPEI survival + time-varying detection
# M2. SPEI survival + constant detection
# M3. Time-varying survival & detection
# M4. Time-varying survival & constant detection
# M5. Constant survival + time-varying detection
# M6. Constant survival + constant detection
# M7. SPEI w/ sex-specific intercepts + time-varying detection
# M8. SPEI w/ sex-specific slopes & intercepts + time-varying detection

## M1: SPEI + time-varying detection ----------

# parameters
params <- c("alpha_phi", "beta_phi", "phi", "p", "p_mean", "log_lik")

# initial values
inits <- function() {
  list(
    alpha_phi = 2,
    beta_phi = 0,
    alpha_p = rep(0, n_years),
    z = z_inits
  )
}

# fit model
fit_SPEI_Vdet <- jags(data = jags_data,
            inits = inits,
            parameters.to.save = params,
            model.file = "scripts/models/mod_cjs_SPEI_Vdet.txt",
            n.chains = n_chains,
            n.adapt = n_adapt,
            n.iter = n_iter,
            n.burnin = n_burn,
            n.thin = 1,
            parallel = TRUE)

# review
print(fit_SPEI_Vdet)

# clear env
#rm(fit_SPEI_Vdet)
#gc()

## M2: SPEI + constant detection ------------------

# parameters
params <- c("alpha_phi", "beta_phi", "phi", "alpha_p0", "log_lik")


# initial values
inits <- function() {
  list(
    alpha_phi = 2,
    beta_phi = 0,
    alpha_p0 = 0,
    z = z_inits
  )
}

# fit
fit_SPEI_Cdet <- jags(data = jags_data,
            inits = inits,
            parameters.to.save = params,
            model.file = "scripts/models/mod_cjs_SPEI_Cdet.txt",
            n.chains = n_chains,
            n.adapt = n_adapt,
            n.iter = n_iter,
            n.burnin = n_burn,
            n.thin = 1,
            parallel = TRUE)

# review
print(fit_SPEI_Cdet)

# clear env
#rm(fit_SPEI_Cdet)
#gc()


## M3: Time-varying survival & detection --------------

# parameters
params <- c("alpha_phi", "phi", "alpha_p", "p", "p_mean", "log_lik")

# inital values
inits <- function() {
  list(
    alpha_phi = rnorm(n_years - 1, 0, 0.1),
    alpha_p   = rnorm(n_years, 0, 0.1),
    z = z_inits
  )
}

# fit
fit_time_Vdet <- jags(data = jags_data,
                 inits = inits,
                 parameters.to.save = params,
                 model.file = "scripts/models/mod_cjs_Vsurv_Vdet.txt",
                 n.chains = n_chains,
                 n.adapt = n_adapt,
                 n.iter = n_iter,
                 n.burnin = n_burn,
                 n.thin = 1,
                 parallel = TRUE)

# review
print(fit_time_Vdet)

# clear env
#rm(fit_time_Vdet)
#gc()

## M4: Time-varying survival + constant detection -------------

# parameters
params <- c("alpha_phi0",  "phi", "alpha_p", "p", "p_mean", "log_lik")


# inital values
inits <- function() {
  list(
    alpha_phi = rnorm(n_years - 1, 0, 0.1),
    alpha_p   = rnorm(1, 0, 0.1),
    z = z_inits
  )
}

# fit
fit_time_Cdet <- jags(data = jags_data,
                      inits = inits,
                      parameters.to.save = params,
                      model.file = "scripts/models/mod_cjs_Vsurv_Cdet.txt",
                      n.chains = n_chains,
                      n.adapt = n_adapt,
                      n.iter = n_iter,
                      n.burnin = n_burn,
                      n.thin = 1,
                      parallel = TRUE)

# review
print(fit_time_Cdet)
summary(fit_time_Cdet)


# clear env
#rm(fit_time_Cdet)
#gc()

## M5: Constant survival + time-varying detection -----------------

# parameters
params <- c("alpha_phi", "phi", "alpha_p", "p", "p_mean", "log_lik")


# inital values
inits <- function() {
  list(
    alpha_phi = rnorm(1, 0, 0.1),
    alpha_p   = rnorm(n_years, 0, 0.1),
    z = z_inits
  )
}

# fit
fit_Csurv_Vdet <- jags(data = jags_data,
                      inits = inits,
                      parameters.to.save = params,
                      model.file = "scripts/models/mod_cjs_Csurv_Vdet.txt",
                      n.chains = n_chains,
                      n.adapt = n_adapt,
                      n.iter = n_iter,
                      n.burnin = n_burn,
                      n.thin = 1,
                      parallel = TRUE)

# review
print(fit_Csurv_Vdet)
summary(fit_Csurv_Vdet)

# clear env
#rm(fit_Csurv_Vdet)
#gc()

## M6: Constant survival + detection -----------------------

# parameters
params <- c("phi", "p_avail", "p_mean", "log_lik")


# inital values
inits <- function() {
  list(
    alpha_phi = rnorm(1, 0, 0.1),
    alpha_p   = rnorm(1, 0, 0.1),
    z = z_inits
  )
}

# fit
fit_Csurv_Cdet <- jags(data = jags_data,
                      inits = inits,
                      parameters.to.save = params,
                      model.file = "scripts/models/mod_cjs_Csurv_Cdet.txt",
                      n.chains = n_chains,
                      n.adapt = n_adapt,
                      n.iter = n_iter,
                      n.burnin = n_burn,
                      n.thin = 1,
                      parallel = TRUE)

# review
print(fit_Csurv_Cdet)
summary(fit_Csurv_Cdet)

# clear env
#rm(fit_Csurv_Cdet)
#gc()

## M7: SPEI w/ sex-specific intercepts + time-varying detection ----------

# parameters
params <- c("alpha_phi", "beta_phi", "phi", "p", "p_mean", "log_lik")

# update jags_data
jags_data <- list(
  Y       = Y,
  drought = drought,
  n_years = n_years,
  n_ind   = n_ind,
  first   = first,
  miss    = miss,
  sex     = sex,
  n_sex   = 2
)


# initial values
inits <- function() {
  list(
    alpha_phi = c(2, 2), # starting values for the two sexes
    beta_phi  = 0,
    alpha_p   = rep(0, n_years),
    z         = z_inits
  )
}

# fit model
fit_SPEI_SEX_Vdet <- jags(data = jags_data,
                      inits = inits,
                      parameters.to.save = params,
                      model.file = "scripts/models/mod_cjs_SPEI_SEX_Vdet.txt",
                      n.chains = n_chains,
                      n.adapt = n_adapt,
                      n.iter = n_iter,
                      n.burnin = n_burn,
                      n.thin = 1,
                      parallel = TRUE)

# review
print(fit_SPEI_SEX_Vdet)

# clear env
#rm(fit_SPEI_SEX_Vdet)
#gc()


## M8: SPEI w/ sex-specific intercepts/slopes + time-varying detection ----------

# parameters
params <- c("alpha_phi", "beta_phi", "phi", "p", "p_mean", "log_lik")

# update jags_data
jags_data <- list(
  Y       = Y,
  drought = drought,
  n_years = n_years,
  n_ind   = n_ind,
  first   = first,
  miss    = miss,
  sex     = sex,
  n_sex   = 2
)


# initial values
inits <- function() {
  list(
    alpha_phi = c(2, 2), # starting values for the two sexes
    beta_phi  = c(0, 0), # now one per sex
    alpha_p   = rep(0, n_years),
    z         = z_inits
  )
}

# fit model
fit_SPEI_SEX_int_Vdet <- jags(data = jags_data,
                          inits = inits,
                          parameters.to.save = params,
                          model.file = "scripts/models/mod_cjs_SPEI_SEXint_Vdet.txt",
                          n.chains = n_chains,
                          n.adapt = n_adapt,
                          n.iter = n_iter,
                          n.burnin = n_burn,
                          n.thin = 1,
                          parallel = TRUE)

# review
print(fit_SPEI_SEX_int_Vdet)

# clear env
#rm(fit_SPEI_SEX_int_Vdet)
#gc()


# Compare ----------------------------------------

# function for WAIC loo comparison
waic_loo_cjs <- function(fit) {
  # extract per-time log_lik[i,t] columns for each chain
  ll_by_chain <- lapply(fit$samples, function(ch) {
    mat <- as.matrix(ch)
    keep <- grep("^log_lik\\[", colnames(mat))
    if (!length(keep)) return(NULL)
    mat[, keep, drop = FALSE]
  })
  ll_by_chain <- Filter(Negate(is.null), ll_by_chain)
  if (!length(ll_by_chain)) stop("No log_lik columns found in any chain.")
  
  # collapse time within individual
  get_i <- function(nm) as.integer(sub("^log_lik\\[([0-9]+),.*$", "\\1", nm))
  ll_ind_by_chain <- lapply(ll_by_chain, function(mat) {
    i_idx <- get_i(colnames(mat))
    ids <- sort(unique(i_idx[is.finite(i_idx) & i_idx > 0]))
    sapply(ids, function(i) rowSums(mat[, i_idx == i, drop = FALSE]))
  })
  
  # stack chains -> draws x individuals
  log_lik_ind <- do.call(rbind, ll_ind_by_chain)
  chain_id <- rep(seq_along(ll_ind_by_chain), vapply(ll_ind_by_chain, nrow, integer(1)))
  r_eff <- relative_eff(exp(log_lik_ind), chain_id = chain_id)
  
  # both criteria
  list(
    waic = waic(log_lik_ind),
    loo  = loo(log_lik_ind, r_eff = r_eff, moment_match = TRUE)
  )
}


# WAIC & PSIS-LOO
t1 <- waic_loo_cjs(fit_SPEI_Vdet) # **
t2 <- waic_loo_cjs(fit_SPEI_Cdet) 
t3 <- waic_loo_cjs(fit_time_Vdet) # **
t4 <- waic_loo_cjs(fit_time_Cdet) 
t5 <- waic_loo_cjs(fit_Csurv_Vdet) # **
t6 <- waic_loo_cjs(fit_Csurv_Cdet) 
t7 <- waic_loo_cjs(fit_SPEI_SEX_Vdet)
t8 <- waic_loo_cjs(fit_SPEI_SEX_int_Vdet)

# clear env
rm(fit_time_Vdet, fit_SPEI_Cdet, fit_time_Vdet, 
   fit_time_Cdet, fit_Csurv_Vdet, fit_Csurv_Cdet,
   fit_SPEI_SEX_Vdet, fit_SPEI_SEX_int_Vdet)
gc()

# make ELPD table
mods <- list(M1 = t1, 
             M2 = t2, 
             M3 = t3, 
             M4 = t4, 
             M5 = t5, 
             M6 = t6,
             M7 = t7,
             M8 = t8)
loos <- lapply(mods, `[[`, "loo")
cmp <- loo::loo_compare(loos) # rows are sorted best>>>worst
cmp %>% as_tibble() %>% mutate(perf = 2*se_diff)
ord <- rownames(cmp)
elpd <- sapply(loos, function(x) x$estimates["elpd_loo","Estimate"])
se   <- sapply(loos, function(x) x$estimates["elpd_loo","SE"])
p_loo <- sapply(loos, function(z) z$estimates["p_loo","Estimate"])
k_gt_07 <- sapply(loos, function(z) mean(z$diagnostics$pareto_k > 0.7, na.rm=TRUE))
k_gt_1  <- sapply(loos, function(z) mean(z$diagnostics$pareto_k > 1.0, na.rm=TRUE))

tab <- data.frame(
  model      = ord,
  elpd_loo   = unname(elpd[ord]),
  se_elpd    = unname(se[ord]),
  # ΔELPD from best (positive = worse than best)
  delta_elpd = unname(-cmp[,"elpd_diff"]),
  se_delta   = unname(cmp[,"se_diff"]),
  p_loo        = unname(p_loo[ord]),
  pct_k_gt_0.7 = round(100 * k_gt_07[ord], 1),
  pct_k_gt_1   = round(100 * k_gt_1[ord], 1),
  rank       = seq_along(ord),
  row.names  = NULL
)
tab




