library(tidyverse)
library(boot)
library(coda)
library(nimble)

# Load functions ----------------------------------------------------------

# Import workspace --------------------------------------------------------
load("data output/model_data.RData")

# Prepare data ------------------------------------------------------------

# Define model ------------------------------------------------------------
mammal_model <- nimbleCode({

  ## Ecological model for latent occurrence z (process model) ##
  for(k in 1:nspec){            # Loop over species
    for (i in 1:nsite) {         # Loop over sites

      z[i,k] ~ dbern(psi[i,k])
      logit(psi[i,k]) <- lpsi[k] + beta[1]*equals(veg[i],2) # Fixed effect of open under-story

    }
  }

  ## Observation model for observed data Y ##

  for(k in 1:nspec){            # Loop over species
   for (i in 1:nsite) {         # Loop over sites

      logit(p[i,k]) <- lp[k] + alpha[1]*equals(season[i],2) # Fixed effect of rainy season

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

  }

  for(k in 1:1){                # 1 term in detection model
    alpha[k] ~ dnorm(0, 0.1)
  }

  for(k in 1:1){               # 1 term in occupancy model
    beta[k] ~ dnorm(0, 0.1)
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

})


# Monitored parameters ----------------------------------------------------
params <-  c("mu.psi","mu.p","alpha","beta","lpsi","lp","numspp")

# Assign constants --------------------------------------------------------
occ_mod_consts <- list(nspec = n_spp,
                       nsite = n_sites)

# Compile data ------------------------------------------------------------
occ_mod_data <- list(ysum = Y,
                     J = J,
                     veg = under_vec,
                     season = season_vec)

# Initial values ----------------------------------------------------------
zst <- array(1, dim = c(n_sites,n_spp)) # Observed occurrence as inits for z
zst[is.na(zst)] <- 1

occ_mod_inits <- list(z = zst,
                      mu.psi = runif(1),
                      mu.p = runif(1),
                      alpha = rnorm(1),
                      beta = rnorm(1)
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
runMCMC_samples <- runMCMC(compile_occ_MCMC,
                           nburnin = 5000,
                           niter = 20000,
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

# Write summary to file ---------------------------------------------------
runMCMC_samples$summary$all.chains %>%
  as_tibble %>%
  mutate(param = rownames(runMCMC_samples$summary$all.chains)) %>%
  select(param, everything()) %>%
  write_csv("data output/mammal_mcmc_out.csv")

# Write posteriors to file ------------------------------------------------
saveRDS(runMCMC_samples, "data output/mammal_mcmc_samples.rds")

# MCMC plotting -----------------------------------------------------------
MCMCvis::MCMCtrace(runMCMC_samples$samples,
          pdf = TRUE,
          filename = "data output/nimble_traceplots.pdf",
          ind = TRUE,
          n.eff = TRUE,
          Rhat = TRUE)