library(boot)
library(coda)
library(tidyverse)
library(raster)
library(grDevices)
library(glue)

# Define spatial scale ----------------------------------------------------
spatscale <- 250 # meters

# Read in basic JAGS model workspace --------------------------------------
# spatdir <- glue("data output/modelout_{spatscale}m")

# Updated analysis with new covariates
spatdir <- glue("data output/modelout_{spatscale}m_v2")

load(glue("data output/model_data_{spatscale}m.RData"))
runMCMC_samples <- readRDS(glue("{spatdir}/mammal_mcmc_samples_{spatscale}m_BPV.rds"))
mod <- as.matrix(runMCMC_samples$samples)
dim(mod)

# Read in all pentad covs, scale to create predictor ----------------------
LMs_df <- read.csv(glue("data input/Landscape_metrics_values_{spatscale}m.csv"))

covs_all_sc <- LMs_df %>%
  mutate_at(vars(Produ:Infra),scale) %>%
  mutate_at(vars(Produ:Infra),as.numeric)
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
coeffs <- c("lpsi","beta1","beta2","beta3","beta4","beta5","beta6","beta7","beta8")
coeff_out <- list()

for(i in 1:length(coeffs)){
  ref <- which(grepl(coeffs[i],att) &
                 !grepl(c("sd"),att) & !grepl(c("mu"),att)) # Exclude sd and mu
  coeff_out[[i]] <- mod[,ref]
}
names(coeff_out) <- c("LPSI","BETA1","BETA2","BETA3","BETA4","BETA5","BETA6","BETA7","BETA8")

## Check
coeff_out$LPSI[1:5,1:5]
coeff_out$BETA8[1:5,1:5]

# Occupancy prediction ----------------------------------------------------
for (s in 1:n_spp){
  print.vec <- seq(0,n_spp, by = 1)
  if(s %in% print.vec){cat(paste("\nSpecies ", s))}
  for(i in 1:n_pix){
    for(u in 1:length(select_samp)){
      psi <- plogis(
        coeff_out[["LPSI"]][select_samp[u],s] +
          coeff_out[["BETA1"]][select_samp[u],s] * covs_all_sc$Produ[i] +
          coeff_out[["BETA2"]][select_samp[u],s] * covs_all_sc$Pheno[i] +
          coeff_out[["BETA3"]][select_samp[u],s] * covs_all_sc$Struc[i] +
          coeff_out[["BETA4"]][select_samp[u],s] * covs_all_sc$Slope[i] +
          coeff_out[["BETA5"]][select_samp[u],s] * covs_all_sc$Elev[i] +
          coeff_out[["BETA6"]][select_samp[u],s] * covs_all_sc$Water[i]+
          coeff_out[["BETA7"]][select_samp[u],s] * covs_all_sc$SBNP[i] +
          coeff_out[["BETA8"]][select_samp[u],s] * covs_all_sc$Infra[i])
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
write_csv(psi_all[[1]],glue("{spatdir}/psi_mean_{spatscale}m.csv"))
write_csv(psi_all[[2]],glue("{spatdir}/psi_sd_{spatscale}m.csv"))
write_csv(psi_all[[3]],glue("{spatdir}/psi_low_{spatscale}m.csv"))
write_csv(psi_all[[4]],glue("{spatdir}/psi_upp_{spatscale}m.csv"))

# Generate rasters with results -------------------------------------------
psi_mean <- read_csv(glue("{spatdir}/psi_mean_{spatscale}m.csv"))
raster_template <- raster(glue("data input/NDVI_produ_{spatscale}.tif"))
dir.create(glue("{spatdir}/spp_psi_mean"))

breaks_seq <- seq(0, 1, by = 0.1) # Color gradient with 10 values between 0.0 and 1.0
col_occupancy <- colorRampPalette(c("yellow", "magenta", "blue", "black"))(length(breaks_seq))

for(k in 2:length(psi_mean)){

 psi_mean_mx <- as.matrix(psi_mean[,k])
 psi_mean_raster <- setValues(raster_template, psi_mean_mx)
 name_spp <- colnames(psi_mean_mx)

 writeRaster(psi_mean_raster,
             filename=glue("{spatdir}/spp_psi_mean/{name_spp}_psi_mean_{spatscale}m.tif"),
             format="GTiff", datatype="FLT4S",
             overwrite=TRUE)

 jpeg(glue("{spatdir}/spp_psi_mean/{name_spp}_psi_mean_{spatscale}m.jpg"),
      width = 1000, height = 1000, res = 200)
  plot(psi_mean_raster, breaks=breaks_seq, col=col_occupancy,
       main=paste0(name_spp, " mean occupancy probability"), cex.main = 0.8)
 dev.off()

}

saveRDS(psi_all, glue("{spatdir}/psi_all_data_{spatscale}m.RDS"))
