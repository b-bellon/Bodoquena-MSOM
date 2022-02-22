library(tidyverse)
library(nimble)
library(glue)
library(bayesboot)
library(bayestestR)

# Define spatial scale ----------------------------------------------------
spatscale <- 250 # meters

# Load data ---------------------------------------------------------------
# spatdir <- glue("data output/modelout_{spatscale}m")

# Updated analysis with new covariates
spatdir <- glue("data output/modelout_{spatscale}m_v2")

runMCMC_samples <- readRDS(glue("{spatdir}/mammal_mcmc_samples_{spatscale}m_BPV.rds"))
mod <- as.matrix(runMCMC_samples$samples)
load(glue("data output/model_data_{spatscale}m.RData"))

# Common name search ------------------------------------------------------

# spp_common <- taxize::sci2comm(spp, simplify = TRUE)
# spp_list <- simplify(spp_common) %>%
#   enframe()
#
# spp_list <- spp %>%
#   enframe() %>%
#   left_join(spp_list, by = c("value" = "name")) %>%
#   rename(spp_sci = value, spp_com = value.y)
#
# spp_list %>%
#   print(n = 29)

# Create plotting data frames ---------------------------------------------
plot_data <- runMCMC_samples$summary$all.chains %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(param = row.names(runMCMC_samples$summary$all.chains)) %>%
  dplyr::select(param, mean, x95_percent_ci_low, x95_percent_ci_upp) %>%
  rename(lower = x95_percent_ci_low, upper = x95_percent_ci_upp)

mup <- plot_data %>%
  filter(str_detect(param,"^mu.p$"))

mupsi <- plot_data %>%
  filter(str_detect(param,"mu.psi"))

mubeta <- plot_data %>%
  filter(str_detect(param, "mu.beta")) %>%
  mutate(param = as.factor(c("produ", "pheno", "struc",
                             "slope", "elev", "water", "sbnp", "infra"))) %>%
  mutate(param = fct_relevel(param, "produ", "pheno", "struc",
                             "slope", "elev", "water", "sbnp", "infra"))

psi_spp <- plot_data %>%
  filter(str_detect(param,"lpsi")) %>%
  mutate_at(vars(mean:upper), boot::inv.logit) %>%
  mutate(spp = spp)

p_spp <- plot_data %>%
  filter(str_detect(param,"lp\\[")) %>%
  mutate_at(vars(mean:upper), boot::inv.logit) %>%
  mutate(spp = spp)

alphas <- plot_data %>%
  filter(str_detect(param,"alpha")) %>%
  mutate(param = c("Rain season","Open understory"))

lpsi <- plot_data %>%
  filter(str_detect(param, "^lpsi")) %>%
  pull(mean)

beta_extract <- function(x){
  plot_data %>%
    filter(str_detect(param, x)) %>%
    pull(mean)
}

beta_list <- map(c("^beta1","^beta2","^beta3",
                   "^beta4","^beta5","^beta6",
                   "^beta7","^beta8"),
                 beta_extract)

# Detection coefficients --------------------------------------------------
ggplot(alphas, aes(x = param, y = mean, ymin = lower, ymax = upper)) +
  geom_pointrange() +
  scale_y_continuous(limits = c(-0.5,0.5))+
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  ylab("Estimate (+- 95% CI)") + xlab("")+
  theme_classic()+
  labs(title = "Detection")+
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"))
ggsave(glue("{spatdir}/alphas_{spatscale}m.jpg"))

# Occupancy coefficients --------------------------------------------------
ggplot(mubeta, aes(x = fct_rev(param), y = mean, ymin = lower, ymax = upper)) +
  geom_pointrange() +
  scale_y_continuous(limits = c(-1,1))+
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  ylab("Estimate (+- 95% CI)") + xlab("")+
  theme_classic()+
  labs(title = glue("Community occupancy at {spatscale}m scale"))+
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"))
ggsave(glue("{spatdir}/mubeta_{spatscale}m.jpg"))


# HDI beta & alpha plots --------------------------------------------------
produ <- mod[,"mu.beta1"]
pheno <- mod[,"mu.beta2"]
struc <- mod[,"mu.beta3"]
slope <- mod[,"mu.beta4"]
elev <- mod[,"mu.beta5"]
water <- mod[,"mu.beta6"]
sbnp <- mod[,"mu.beta7"]
infra <- mod[,"mu.beta8"]

jpeg(glue("{spatdir}/HDI_mubeta_{spatscale}m.jpg"),width = 2000, height = 2200, res = 200)
par(mfrow = c(4,2))
bayesboot::plotPost(produ, credMass = 0.89, xlim = c(-0.5,0.8))
bayesboot::plotPost(pheno, credMass = 0.89, xlim = c(-0.5,0.8))
bayesboot::plotPost(struc, credMass = 0.89, xlim = c(-0.5,0.8))
bayesboot::plotPost(elev, credMass = 0.89, xlim = c(-0.5,0.8))
bayesboot::plotPost(slope, credMass = 0.89, xlim = c(-0.5,0.8))
bayesboot::plotPost(water, credMass = 0.89, xlim = c(-0.5,0.8))
bayesboot::plotPost(sbnp, credMass = 0.89, xlim = c(-0.5,0.8))
bayesboot::plotPost(infra, credMass = 0.89, xlim = c(-0.5,0.8))
dev.off()

rain_season <- mod[,"alpha[1]"]
open_understory <- mod[,"alpha[2]"]

jpeg(glue("{spatdir}/HDI_alpha_{spatscale}m.jpg"),width = 1000, height = 1100, res = 200)
par(mfrow = c(2,1))
bayesboot::plotPost(rain_season, credMass = 0.89, xlim = c(-0.3,0.3))
bayesboot::plotPost(open_understory, credMass = 0.89, xlim = c(-0.3,0.3))
dev.off()

# Species occupancy probability -------------------------------------------
psi_spp %>%
  ggplot(aes(x = fct_reorder(spp,mean), y = mean, ymin = lower, ymax = upper)) +
  geom_pointrange(size = 0.6) + theme_bw() + scale_y_continuous(limits = c(0,1))+
  coord_flip() + geom_hline(yintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = mupsi$mean, color = "red", size = 1) +
  geom_hline(yintercept = mupsi$lower, color = alpha("blue",0.3), size = 1) +
  geom_hline(yintercept = mupsi$upper, color = alpha("blue",0.3), size = 1) +
  ylab("Occupancy probability") + xlab("")+
  theme(axis.text.y=element_text(size=13),
        axis.text.x=element_text(size=13),
        axis.title.x = element_text(size=13))
ggsave(glue("{spatdir}/occprob_{spatscale}m.jpg"))

# Species detection probability -------------------------------------------
p_spp %>%
  ggplot(aes(x = fct_reorder(spp,mean), y = mean, ymin = lower, ymax = upper)) +
  geom_pointrange(size = 0.6) + theme_bw() + scale_y_continuous(limits = c(0,1))+
  coord_flip() + geom_hline(yintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = mup$mean, color = "red", size = 1) +
  geom_hline(yintercept = mup$lower, color = alpha("blue",0.3), size = 1) +
  geom_hline(yintercept = mup$upper, color = alpha("blue",0.3), size = 1) +
  ylab("Detection probability") + xlab("")+
  theme(axis.text.y=element_text(size=13),
        axis.text.x=element_text(size=13),
        axis.title.x = element_text(size=13))
ggsave(glue("{spatdir}/detecprob_{spatscale}m.jpg"))

# Species richness --------------------------------------------------------
plot_data %>%
  filter(str_detect(param,"numspp")) %>%
  ggplot(aes(x = mean))+
  geom_histogram(binwidth = 0.5,
                 fill = "grey",
                 color = "black")+
  xlab("Species richness")

ggsave(glue("{spatdir}/specrich_{spatscale}m.jpg"))

# Random effect curve plots -----------------------------------------------

## Used to back-transform predictions
prodvar <- glue("produ_{spatscale}")
phenovar <- glue("pheno_{spatscale}")
strucvar <- glue("struc_{spatscale}")

beta1_sc <- scale(station_data %>% pull({{prodvar}}))
beta2_sc <- scale(station_data %>% pull({{phenovar}}))
beta3_sc <- scale(station_data %>% pull({{strucvar}}))
beta4_sc <- scale(station_data %>% pull(slope))
beta5_sc <- scale(station_data %>% pull(elevation))
beta6_sc <- scale(station_data %>% pull(dist_water))
beta7_sc <- scale(station_data %>% pull(dist_sbnp))
beta8_sc <- scale(station_data %>% pull(dist_infra))

## Predict over these range of values
x <- seq(-2,2,0.01)

## Function to extract occupancy probability over x range
pred_data_extract <- function(sppvec){


  x1 <- x * attr(beta1_sc, 'scaled:scale') + attr(beta1_sc, 'scaled:center')
  pred1 <- plogis(lpsi[sppvec] + beta_list[[1]][sppvec] * x)

  x2 <- x * attr(beta2_sc, 'scaled:scale') + attr(beta2_sc, 'scaled:center')
  pred2 <- plogis(lpsi[sppvec] + beta_list[[2]][sppvec] * x)

  x3 <- x * attr(beta3_sc, 'scaled:scale') + attr(beta3_sc, 'scaled:center')
  pred3 <- plogis(lpsi[sppvec] + beta_list[[3]][sppvec] * x)

  x4 <- x * attr(beta4_sc, 'scaled:scale') + attr(beta4_sc, 'scaled:center')
  pred4 <- plogis(lpsi[sppvec] + beta_list[[4]][sppvec] * x)

  x5 <- x * attr(beta5_sc, 'scaled:scale') + attr(beta5_sc, 'scaled:center')
  pred5 <- plogis(lpsi[sppvec] + beta_list[[5]][sppvec] * x)

  x6 <- x * attr(beta6_sc, 'scaled:scale') + attr(beta6_sc, 'scaled:center')
  pred6 <- plogis(lpsi[sppvec] + beta_list[[6]][sppvec] * x)

  x7 <- x * attr(beta7_sc, 'scaled:scale') + attr(beta7_sc, 'scaled:center')
  pred7 <- plogis(lpsi[sppvec] + beta_list[[7]][sppvec] * x)

  x8 <- x * attr(beta8_sc, 'scaled:scale') + attr(beta8_sc, 'scaled:center')
  pred8 <- plogis(lpsi[sppvec] + beta_list[[8]][sppvec] * x)

  data <- tibble(mean = c(pred1,pred2,pred3,pred4,pred5,pred6,pred7,pred8),
                 betas = c(
                   rep("a) Productivity", length(x)),
                   rep("b) Phenology", length(x)),
                   rep("c) Structure", length(x)),
                   rep("d) Slope", length(x)),
                   rep("e) Elevation", length(x)),
                   rep("f) Dist_water", length(x)),
                   rep("g) Dist_SBNP", length(x)),
                   rep("h) Dist_infra", length(x))
                   ),
                 spp = as.character(sppvec),
                 xdat = c(x1,x2,x3,x4,x5,x6,x7,x8)
                 )
  return(data)
}

pred_data <- map_df(1:n_spp,
                 pred_data_extract)

## Check
pred_data %>%
  count(betas, spp)

pred_data %>%
  ggplot(aes(x = xdat, y = mean, color = spp))+
  geom_line() +
  scale_color_manual(values=c(rep("black",n_spp))) +
  theme(legend.position = "none") +
  xlab("")+
  ylab("Occupancy probability")+
  facet_wrap(~ betas, scales = "free_x")

ggsave(glue("{spatdir}/Occ prob curves_{spatscale}m.jpg"))


## Extracting community level predictions to add to plots above
## Not complete, still a work in progress...
# mu.lpsi <- mupsi %>%
#   pull(mean) %>%
#   logit(.)
#
# mubetas <- mubeta %>%
#   pull(mean)
#
# mu.lpsi_upp <- mupsi %>%
#   pull(upper) %>%
#   logit(.)
#
# mubetas_upp <- mubeta %>%
#   pull(upper)
#
# mubetas_lower <- mubeta %>%
#   pull(lower)
#
# hist(pheno_vec)
# hist(struc_vec)
# hist(produ_vec)
#
# x <- seq(-2,2,0.01)
# pred1 <- plogis(mu.lpsi[1] + mubetas[1] * x)
# pred2 <- plogis(mu.lpsi[1] + mubetas[2] * x)
# pred3 <- plogis(mu.lpsi[1] + mubetas[3] * x)


# Write workspace ---------------------------------------------------------
save.image(glue("{spatdir}/plotting_data_{spatscale}m.RData"))
