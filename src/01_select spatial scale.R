## Libraries
library(tidyverse)
library(glue)

# Define spatial scale ----------------------------------------------------
spatscale <- 250 # meters

# Load data ---------------------------------------------------------------
load("data output/model_data.RData")

# Process landscape variables ---------------------------------------------
prodvar <- glue("produ_{spatscale}")
phenovar <- glue("pheno_{spatscale}")
strucvar <- glue("struc_{spatscale}")

produ_vec <- as.numeric(scale(station_data %>% pull({{prodvar}}))[,1])
pheno_vec <- as.numeric(scale(station_data %>% pull({{phenovar}}))[,1])
struc_vec <- as.numeric(scale(station_data %>% pull({{strucvar}}))[,1])

slope_vec <- as.numeric(scale(station_data %>% pull(slope))[,1])
elev_vec <- as.numeric(scale(station_data %>% pull(elevation))[,1])
water_vec <- as.numeric(scale(station_data %>% pull(dist_water))[,1])
sbnp_vec <- as.numeric(scale(station_data %>% pull(dist_sbnp))[,1])
infra_vec <- as.numeric(scale(station_data %>% pull(dist_infra))[,1])

# Write workspace ---------------------------------------------------------
save.image(glue("data output/model_data_{spatscale}m.RData"))

