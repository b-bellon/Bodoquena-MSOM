library(boot)
library(coda)
library(tidyverse)

# Read in basic JAGS model workspace --------------------------------------
load("data output/model_data.RData")
runMCMC_samples <- readRDS("data output/mammal_mcmc_samples.rds")
mod <- as.matrix(runMCMC_samples$samples)
dim(mod)

# Read in all pentad covs, scale to create predictor ----------------------
covs_all_sc <- LM_df %>%
  mutate_at(vars(Produ:Struc),scale) %>%
  mutate_at(vars(Produ:Struc),as.numeric)
covs_all_sc

# Create array ------------------------------------------------------------
(n_samp <- nrow(mod))                      # 18000 is too many!
select_samp <- sort(sample(1:n_samp, 100)) # Chose random sample of 100 posterior draws (can increase this)
(n_samp <- length(select_samp))
(n_pix <- nrow(covs_all_sc))

# Create posterior predictive distribution for psi for each species -------
str(psi_pen <- array(NA, dim = c(n_pix, n_samp,n_spp)))

# Acccess node names and extract posterior samples ------------------------
att <- attributes(mod)$dimnames[2][[1]][]
coeffs <- c("lpsi","beta1","beta2","beta3")
coeff_out <- list()

for(i in 1:length(coeffs)){
  ref <- which(grepl(coeffs[i],att) &
                 !grepl(c("sd"),att) & !grepl(c("mu"),att)) # Exclude sd and mu
  coeff_out[[i]] <- mod[,ref]
}
names(coeff_out) <- c("LPSI","BETA1","BETA2","BETA3")

## Check
coeff_out$LPSI[1:5,1:5]
coeff_out$BETA3[1:5,1:5]

# Occupancy prediction ----------------------------------------------------

for (s in 1:n_spp){
  print.vec <- seq(0,n_spp, by = 1)
  if(s %in% print.vec){cat(paste("\nSpecies ", s))}
  for(i in 1:n_pix){
    for(u in 1:length(select_samp)){
      psi <- plogis(coeff_out[["LPSI"]][select_samp[u],s] +
                      coeff_out[["BETA1"]][select_samp[u],s] * covs_all_sc$Produ[i] +
                      coeff_out[["BETA2"]][select_samp[u],s] * covs_all_sc$Pheno[i] +
                      coeff_out[["BETA3"]][select_samp[u],s] * covs_all_sc$Struc[i])
      psi_pen[i,u,s] <- psi
    }
  }
}

# Extract 4 measures of psi -----------------------------------------------

## psi:: mean and sd
extract_psi_meansd <- function(psifun){
  df <- apply(psi_pen, c(1,3), psifun, na.rm = TRUE)  %>%
    as_tibble  %>%
    magrittr::set_colnames(spp) %>%
    mutate(ID_pix = covs_all_sc$ID_pix) %>%
    dplyr::select(ID_pix, everything())
}
psifuncs <- list(list(mean),list(sd))
psi_spdf <- invoke_map(.f = extract_psi_meansd,
                       .x = psifuncs)
names(psi_spdf) <- c("mean","sd")

## psi:: 95% CIs
extract_psi_quart <- function(lev){
  df <- apply(psi_pen, c(1,3), function (x) quantile(x, probs = lev)) %>%
    as_tibble  %>%
    magrittr::set_colnames(spp) %>%
    mutate(ID_pix = covs_all_sc$ID_pix) %>%
    dplyr::select(ID_pix, everything())
}

lev <- c(0.025,0.975)
psi_spdf_quart <- map(.x = lev,
                      .f = extract_psi_quart)
names(psi_spdf_quart) <- c("low","upp")

## Join lists
psi_all <- c(psi_spdf,psi_spdf_quart)
names(psi_all)

# Write psi measures to file ----------------------------------------------
write_csv(psi_all[[1]],"data output/psi_mean.csv")
write_csv(psi_all[[2]],"data output/psi_sd.csv")
write_csv(psi_all[[3]],"data output/psi_low.csv")
write_csv(psi_all[[4]],"data output/psi_upp.csv")
