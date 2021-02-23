library(tidyverse)
library(boot)
library(coda)
library(nimble)
library(glue)

# Define spatial scale ----------------------------------------------------
spatscale <- 750 # meters

# Create dedicated output directory ---------------------------------------
spatdir <- glue("data output/modelout_{spatscale}m")

if(dir.exists(spatdir)){
  print("Directory exists")
  } else {
    dir.create(spatdir)
    print("Directory created")
  }

# Import scale-specific workspace -----------------------------------------
load(glue("data output/model_data_{spatscale}m.RData"))

# Define model ------------------------------------------------------------
mammal_model <- nimbleCode({

  ## Ecological model for latent occurrence z (process model) ##
  for(k in 1:nspec){            # Loop over species
    for (i in 1:nsite) {         # Loop over sites

      z[i,k] ~ dbern(psi[i,k])
      logit(psi[i,k]) <- lpsi[k] + beta1[k]*produ[i] + beta2[k]*pheno[i] + beta3[k]*struc[i]

    }
  }

  ## Observation model for observed data Y ##

  for(k in 1:nspec){            # Loop over species
   for (i in 1:nsite) {         # Loop over sites

      logit(p[i,k]) <- lp[k] + alpha[1]*equals(season[i],2) + alpha[2]*equals(veg[i],2) # Fixed effect of rain season and open/under story

      mup[i,k] <- z[i,k] * p[i,k]
      ysum[i,k] ~ dbin(mup[i,k], J[i])

    }
  }

  # Species richness
  for (i in 1:nsite) {            # Loop over sites
    numspp[i] <- sum(z[i,])       # Add up number of occurring species at each site
  }

  ## Priors ##
  for(k in 1:nspec){                    # Loop over species
    lpsi[k] ~ dnorm(mu.lpsi, tau.lpsi)
    lp[k] ~ dnorm(mu.lp, tau.lp)
    beta1[k] ~ dnorm(mu.beta1, tau.beta1)    # produ
    beta2[k] ~ dnorm(mu.beta2, tau.beta2)    # pheno
    beta3[k] ~ dnorm(mu.beta3, tau.beta3)    # struc
  }

  for(k in 1:2){                # 2 terms in detection model
    alpha[k] ~ dnorm(0, 0.1)
  }

  ## Random intercept priors ##
  mu.psi ~ dunif(0,1)
  mu.lpsi <- logit(mu.psi)
  sd.lpsi ~ dunif(0,5)
  tau.lpsi <- pow(sd.lpsi, -2)

  mu.p ~ dunif(0,1)
  mu.lp <- logit(mu.p)
  sd.lp ~ dunif(0,5)
  tau.lp <- pow(sd.lp, -2)

  mu.beta1 ~ dnorm(0, 0.1)
  x.beta1 ~ dt(0,1,1)
  sd.beta1 <- abs(x.beta1)
  tau.beta1 <- pow(sd.beta1, -2)

  mu.beta2 ~ dnorm(0, 0.1)
  x.beta2 ~ dt(0,1,1)
  sd.beta2 <- abs(x.beta2)
  tau.beta2 <- pow(sd.beta2, -2)

  mu.beta3 ~ dnorm(0, 0.1)
  x.beta3 ~ dt(0,1,1)
  sd.beta3 <- abs(x.beta3)
  tau.beta3 <- pow(sd.beta3, -2)

})


# Monitored parameters ----------------------------------------------------
params <-  c("mu.psi","mu.p","alpha",
             "mu.beta1","mu.beta2","mu.beta3",
             "beta1", "beta2", "beta3",
             "lpsi","lp","numspp","z")

# Assign constants --------------------------------------------------------
occ_mod_consts <- list(nspec = n_spp,
                       nsite = n_sites)

# Compile data ------------------------------------------------------------
occ_mod_data <- list(ysum = Y,
                     J = J,
                     veg = under_vec,
                     season = season_vec,
                     produ = produ_vec,
                     pheno = pheno_vec,
                     struc = struc_vec)

# Initial values ----------------------------------------------------------
zst <- array(1, dim = c(n_sites,n_spp)) # Observed occurrence as inits for z
zst[is.na(zst)] <- 1

occ_mod_inits <- list(z = zst,
                      mu.psi = runif(1),
                      mu.p = runif(1),
                      alpha = rnorm(1),
                      mu.beta1 = rnorm(1),
                      mu.beta2 = rnorm(1),
                      mu.beta3 = rnorm(1),
                      beta1 = rep(rnorm(1), occ_mod_consts$nspec),
                      beta2 = rep(rnorm(1), occ_mod_consts$nspec),
                      beta3 = rep(rnorm(1), occ_mod_consts$nspec)
                      )

# Build the NIMBLE model --------------------------------------------------
s0 <- Sys.time()
occ_model <- nimbleModel(code = mammal_model, name = "Mammal MSOM",
                         constants = occ_mod_consts,
                         data = occ_mod_data, inits = occ_mod_inits,
                         calculate = TRUE) ## set to FALSE for large models
e0 <- Sys.time()
e0-s0

## Build an MCMC object for this
s1 <- Sys.time()
occ_model_MCMC <- buildMCMC(occ_model, monitors = params, enableWAIC = TRUE)
e1 <- Sys.time()
e1-s1

## Compile
s2 <- Sys.time()
comp_occ_mod <- compileNimble(occ_model)
e2 <- Sys.time()
e2-s2

s3 <- Sys.time()
compile_occ_MCMC <- compileNimble(occ_model_MCMC, project = comp_occ_mod)
e3 <- Sys.time()
e3-s3


# Run MCMC ----------------------------------------------------------------

# Number of samples returned will be floor((niter-nburnin)/thin)
s4 <- Sys.time()

# Final run
# burn <- 15000
# iter <- 60000

# Test
burn <- 1000
iter <- 5000

runMCMC_samples <- runMCMC(compile_occ_MCMC,
                           nburnin = burn,
                           niter = iter,
                           nchains = 3,
                           thin = 5,
                           summary = TRUE,
                           samplesAsCodaMCMC = TRUE,
                           WAIC = TRUE)
e4 <- Sys.time()
e4-s4

e4-s4 + e3-s3 + e2-s2 + e1-s1 + e0-s0

# Check output ------------------------------------------------------------
head(runMCMC_samples$summary$all.chains)

## WAIC check
runMCMC_samples$WAIC # 4045.371

# Write summary to file ---------------------------------------------------
runMCMC_samples$summary$all.chains %>%
  as_tibble %>%
  mutate(param = rownames(runMCMC_samples$summary$all.chains)) %>%
  dplyr::select(param, everything()) %>%
  write_csv(glue("{spatdir}/mammal_mcmc_out_{spatscale}m.csv"))

# Write posteriors to file ------------------------------------------------
saveRDS(runMCMC_samples, glue("{spatdir}/mammal_mcmc_samples_{spatscale}m.rds"))

# MCMC plotting -----------------------------------------------------------
MCMCvis::MCMCtrace(runMCMC_samples$samples,
          pdf = TRUE,
          filename = glue("{spatdir}/nimble_traceplots_{spatscale}m.pdf"),
          ind = TRUE,
          n.eff = TRUE,
          Rhat = TRUE)
