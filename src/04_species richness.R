library(boot)
library(coda)
library(tidyverse)

# Read in basic JAGS model workspace --------------------------------------
load("data output/model_data.RData")
runMCMC_samples <- readRDS("data output/mammal_mcmc_samples.rds")
mod <- as.matrix(runMCMC_samples$samples)
dim(mod)

# Read in all pentad covs, scale to create predictor ----------------------
covs_all_sc <- readxl::read_xlsx("data input/Id_pix.xlsx") %>%
  mutate_at(vars(Produ:Struc),scale) %>%
  mutate_at(vars(Produ:Struc),as.numeric)
covs_all_sc

# Create array ------------------------------------------------------------
(n_samp <- nrow(mod))                      # 18000 is too many!
select_samp <- sort(sample(1:n_samp, 100)) # Chose random sample of 100 posterior draws (can increase this)
(n_samp <- length(select_samp))
(n_pix <- nrow(covs_all_sc))

## Create posterior predictive distribution for Z for all pentads
str(zPen <- array(NA, dim = c(n_pix, n_spp, n_samp)))

# Acccess node names and extract posterior samples ------------------------
att <- attributes(mod)$dimnames[2][[1]][]
coeffs <- c("lpsi","beta1","beta2","beta3")
coeff_out <- list()

for(i in 1:length(coeffs)){
  ref <- which(grepl(coeffs[i],att) & !grepl(c("sd"),att) & !grepl(c("mu"),att)) # Exclude sd and mu
  coeff_out[[i]] <- mod[,ref]
}
names(coeff_out) <- c("LPSI","BETA1","BETA2","BETA3")

# Generate occupancy predictions for 100 samples --------------------------
start <- Sys.time()

for(i in 1:n_pix){

  print.vec <- seq(0,n_pix, by = 10000)
  if(i %in% print.vec){cat(paste("\nPixel ", i))}

  for(u in 1:length(select_samp)){
    psi <- plogis(coeff_out[["LPSI"]][select_samp[u],] +
                    coeff_out[["BETA1"]][select_samp[u],] * covs_all_sc$Produ[i] +
                    coeff_out[["BETA2"]][select_samp[u],] * covs_all_sc$Pheno[i] +
                    coeff_out[["BETA3"]][select_samp[u],] * covs_all_sc$Struc[i])
    zPen[i,,u] <- rbinom(n_spp, 1, psi)

  }
}

end <- Sys.time()

## Runtime
end - start

str(zPen)

# Compute posterior distribution of SR - collapse z array -----------------
SR <- apply(zPen, c(1,3), sum, na.rm= TRUE)   # posterior distribution
pmSR <- apply(SR, 1, mean, na.rm = TRUE)      # posterior mean
lowSR <- apply(SR, 1,function (x) quantile(x, probs = 0.025)) # posterior lower bound
uppSR <- apply(SR, 1,function (x) quantile(x, probs = 0.975)) # posterior upper bound
sdSR <- apply(SR, 1, sd, na.rm = TRUE)        # posterior standard deviation

pmSR[1:10];lowSR[1:10];uppSR[1:10]

tibble(ID_pix = covs_all_sc$ID_pix,
       spp_rich_mean = pmSR,
       spp_rich_low = lowSR,
       spp_rich_upp = uppSR,
       spp_rich_sd = sdSR) %>%
  write_csv("data output/spp_richness.csv")

# Save workspace ----------------------------------------------------------
save(list = c(ls()[!ls() %in% c("SR","covs_all_sc")]),
     file = "data output/Species richness maps.RData")


