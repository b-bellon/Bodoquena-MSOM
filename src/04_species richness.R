library(boot)
library(coda)
library(tidyverse)
library(raster)
library(rgdal)
library(grDevices)

# Read in basic JAGS model workspace --------------------------------------
load("data output/model_data.RData")
runMCMC_samples <- readRDS("data output/mammal_mcmc_samples.rds")
mod <- as.matrix(runMCMC_samples$samples)
dim(mod)


# Read LM rasters ---------------------------------------------------------
NDVI_produ <- raster("data input/NDVI_produ_250.tif")
NDVI_pheno <- raster("data input/NDVI_pheno_250.tif")
NDVI_struc <- raster("data input/NDVI_struc_250.tif")


# Extract LM values + pixels IDs ------------------------------------------
NDVI_produ_vec <- as.vector(NDVI_produ)
NDVI_pheno_vec <- as.vector(NDVI_pheno)
NDVI_struc_vec <- as.vector(NDVI_struc)

LMs_df <- as.data.frame(cbind(NDVI_produ_vec, NDVI_pheno_vec, NDVI_struc_vec))
LMs_df$ID_pix <- seq.int(nrow(LMs_df))

colnames(LMs_df) <- c("Produ","Pheno","Struc", "ID_pix")
LMs_df <- LMs_df[, c(4, 1, 2, 3)] #reorder columns

LMs_df[,2]<- as.integer(LMs_df[,2])
LMs_df[,3]<- as.integer(LMs_df[,3])
LMs_df[,4]<- as.integer(LMs_df[,4])

write.csv(LMs_df, "data input/Landscape_metrics_values.csv", row.names = FALSE)

# Read in all pentad covs, scale to create predictor ----------------------
covs_all_sc <- LMs_df %>%
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


# Generate rasters with results -------------------------------------------
pmSR_raster <- setValues(NDVI_produ, pmSR)
sdSR_raster <- setValues(NDVI_produ, sdSR)
lowSR_raster <- setValues(NDVI_produ, lowSR)
uppSR_raster <- setValues(NDVI_produ, uppSR)

# display.brewer.all()
col_richness <- colorRampPalette(c("red", "yellow", "darkgreen", "black"))
col_richness <-col_richness(10)
col_sd <- colorRampPalette(c("white", "papayawhip", "firebrick2"))
col_sd <- col_sd(10)

jpeg("data output/spp_richness_plots.jpg", width = 1000, height = 1000, res = 200)
par(mfrow=c(2,3), mai=c(0.5, 0.3, 0.3, 0.3))
plot(pmSR_raster, col=col_richness, main='Mean spp. richness')
plot(lowSR_raster, col=col_richness, main='Lower limit')
plot(uppSR_raster, col=col_richness,main='Upper limit')
plot(sdSR_raster, col=col_sd, main='Standard Deviation')
dev.off()

writeRaster(pmSR_raster, filename="data output/Mean_spp_richness.tif", format="GTiff", datatype="FLT4S", overwrite=TRUE)
writeRaster(sdSR_raster, filename="data output/SD_spp_richness.tif", format="GTiff", datatype="FLT4S", overwrite=TRUE)
writeRaster(lowSR_raster, filename="data output/Lower_limit_spp_richness.tif", format="GTiff", datatype="FLT4S", overwrite=TRUE)
writeRaster(uppSR_raster, filename="data output/Upper_limit_spp_richness.tif", format="GTiff", datatype="FLT4S", overwrite=TRUE)


# Save workspace ----------------------------------------------------------
save(list = c(ls()[!ls() %in% c("SR","covs_all_sc")]),
     file = "data output/Species richness maps.RData")


