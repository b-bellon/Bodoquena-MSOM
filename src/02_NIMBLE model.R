library(tidyverse)
library(boot)
library(coda)
library(nimble)
library(glue)

# Define spatial scale ----------------------------------------------------
spatscale <- 250 # meters

# Create dedicated output directory ---------------------------------------
# spatdir <- glue("data output/modelout_{spatscale}m")

# Updated analysis with new covariates
spatdir <- glue("data output/modelout_{spatscale}m_v2")

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
      logit(psi[i,k]) <- lpsi[k] + beta1[k]*produ[i] + beta2[k]*pheno[i] +
        beta3[k]*struc[i] + beta4[k]*slope[i] + beta5[k]*elev[i] +
        beta6[k]*water[i] + beta7[k]*sbnp[i] + beta8[k]*infra[i]

      # Calculate residuals
      yres[i,k]<-psi[i,k]-z[i,k]
    }
  }

  ## Observation model for observed data Y ##

  for(k in 1:nspec){            # Loop over species
    for (i in 1:nsite) {         # Loop over sites

      logit(p[i,k]) <- lp[k] + alpha[1]*equals(season[i],2) +
        alpha[2]*equals(veg[i],2) # Fixed effect of rain season and open/under story

      mup[i,k] <- z[i,k] * p[i,k]
      ysum[i,k] ~ dbin(mup[i,k], J[i])

      dev[i,k]<-(ysum[i,k]*log(mup[i,k]+0.00001)+(1-ysum[i,k])*log(1-(mup[i,k]-0.0001)))*-2

      # Predict new observations and compute deviance
      ysum.new[i,k] ~ dbin(mup[i,k], J[i])
      dev.sim[i,k]<- (ysum.new[i,k]*log(mup[i,k]+0.00001)+(1-ysum.new[i,k])*log(1-(mup[i,k]-0.0001)))*-2

      # Save this node to create plots
      p_val[i,k] <- step(dev.sim[i,k] - dev[i,k])

    }
  }

  # Estimate Bayesian p-value
  sum.dev<-sum(dev[1:nsite,1:nspec])
  sum.dev.sim<-sum(dev.sim[1:nsite,1:nspec])
  test<-step(sum.dev.sim - sum.dev)

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
    beta4[k] ~ dnorm(mu.beta4, tau.beta4)    # slope
    beta5[k] ~ dnorm(mu.beta5, tau.beta5)    # elev
    beta6[k] ~ dnorm(mu.beta6, tau.beta6)    # water
    beta7[k] ~ dnorm(mu.beta7, tau.beta7)    # sbnp
    beta8[k] ~ dnorm(mu.beta8, tau.beta8)    # infra
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

  mu.beta4 ~ dnorm(0, 0.1)
  x.beta4 ~ dt(0,1,1)
  sd.beta4 <- abs(x.beta4)
  tau.beta4 <- pow(sd.beta4, -2)

  mu.beta5 ~ dnorm(0, 0.1)
  x.beta5 ~ dt(0,1,1)
  sd.beta5 <- abs(x.beta5)
  tau.beta5 <- pow(sd.beta5, -2)

  mu.beta6 ~ dnorm(0, 0.1)
  x.beta6 ~ dt(0,1,1)
  sd.beta6 <- abs(x.beta6)
  tau.beta6 <- pow(sd.beta6, -2)

  mu.beta7 ~ dnorm(0, 0.1)
  x.beta7 ~ dt(0,1,1)
  sd.beta7 <- abs(x.beta7)
  tau.beta7 <- pow(sd.beta7, -2)

  mu.beta8 ~ dnorm(0, 0.1)
  x.beta8 ~ dt(0,1,1)
  sd.beta8 <- abs(x.beta8)
  tau.beta8 <- pow(sd.beta8, -2)

})

# Monitored parameters ----------------------------------------------------
params <-  c("mu.psi","mu.p","alpha",
             "mu.beta1","mu.beta2","mu.beta3",
             "mu.beta4","mu.beta5","mu.beta6",
             "mu.beta7","mu.beta8",
             "beta1", "beta2", "beta3",
             "beta4", "beta5", "beta6",
             "beta7", "beta8",
             "lpsi","lp","numspp","z",
             "sum.dev", "sum.dev.sim",
             "p_val",
             "ysum.new","ysum",
             "yres",
             "test"
)

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
                     struc = struc_vec,
                     slope = slope_vec,
                     elev = elev_vec,
                     water = water_vec,
                     sbnp = sbnp_vec,
                     infra = infra_vec)

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
                      mu.beta4 = rnorm(1),
                      mu.beta5 = rnorm(1),
                      mu.beta6 = rnorm(1),
                      mu.beta7 = rnorm(1),
                      mu.beta8 = rnorm(1),
                      beta1 = rep(rnorm(1), occ_mod_consts$nspec),
                      beta2 = rep(rnorm(1), occ_mod_consts$nspec),
                      beta3 = rep(rnorm(1), occ_mod_consts$nspec),
                      beta4 = rep(rnorm(1), occ_mod_consts$nspec),
                      beta5 = rep(rnorm(1), occ_mod_consts$nspec),
                      beta6 = rep(rnorm(1), occ_mod_consts$nspec),
                      beta7 = rep(rnorm(1), occ_mod_consts$nspec),
                      beta8 = rep(rnorm(1), occ_mod_consts$nspec)

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
# nthin <- 5

# Test
burn <- 1000
iter <- 5000
nthin <- 5

runMCMC_samples <- runMCMC(compile_occ_MCMC,
                           nburnin = burn,
                           niter = iter,
                           nchains = 3,
                           thin = nthin,
                           summary = TRUE,
                           samplesAsCodaMCMC = TRUE,
                           WAIC = TRUE)
e4 <- Sys.time()
e4-s4

e4-s4 + e3-s3 + e2-s2 + e1-s1 + e0-s0

# Check output ------------------------------------------------------------
head(runMCMC_samples$summary$all.chains)

## WAIC check
runMCMC_samples$WAIC # 250m = 4056.404 || 750m = 4049.205 || 1250m = 4054.412

# Write summary to file ---------------------------------------------------
runMCMC_samples$summary$all.chains %>%
  as_tibble %>%
  mutate(param = rownames(runMCMC_samples$summary$all.chains)) %>%
  dplyr::select(param, everything()) %>%
  write_csv(glue("{spatdir}/mammal_mcmc_out_{spatscale}m_BPV.csv"))

# Write posteriors to file ------------------------------------------------
saveRDS(runMCMC_samples, glue("{spatdir}/mammal_mcmc_samples_{spatscale}m_BPV.rds"))

# MCMC plotting -----------------------------------------------------------
MCMCvis::MCMCtrace(runMCMC_samples$samples,
                   pdf = TRUE,
                   filename = glue("{spatdir}/nimble_traceplots_{spatscale}m_BPV.pdf"),
                   ind = TRUE,
                   n.eff = TRUE,
                   Rhat = TRUE)
