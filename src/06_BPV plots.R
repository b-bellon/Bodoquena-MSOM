## Libraries
library(tidyverse)
library(boot)
library(coda)
library(nimble)
library(glue)

## ________________________________________________________________________

# Define spatial scale ----------------------------------------------------
spatscale <- 250 # meters

# Read in data ------------------------------------------------------------
# spatdir <- glue("data output/modelout_{spatscale}m")

# Updated analysis with new covariates
spatdir <- glue("data output/modelout_{spatscale}m_v2")

load(glue("data output/model_data_{spatscale}m.RData"))

runMCMC_samples <- readRDS(glue("{spatdir}/mammal_mcmc_samples_{spatscale}m_BPV.rds"))
mod <- as.matrix(runMCMC_samples$samples)
dim(mod)

att <- attributes(mod)$dimnames[2][[1]][]
att

# Extract BPV and plots ---------------------------------------------------
bpv_posterior <- mod[,"test"]
round(mean(bpv_posterior),2)
table(bpv_posterior)

par(mfrow = c(2,2))

dev <- mod[,"sum.dev"]
hist(dev)

dev.sim <- mod[,"sum.dev.sim"]
hist(dev.sim)

plot(dev,dev.sim, main = glue("Deviance at {spatscale}m scale"), xlim = c(7000,10000),  ylim = c(7000,10000))
abline(0,1, col = "red", lwd = 2)

