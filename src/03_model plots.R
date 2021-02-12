library(tidyverse)
library(nimble)

# Load data ---------------------------------------------------------------
runMCMC_samples <- readRDS("data output/mammal_mcmc_samples.rds")
load("data output/model_data.RData")

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
  select(param, mean, x95_percent_ci_low, x95_percent_ci_upp) %>%
  rename(lower = x95_percent_ci_low, upper = x95_percent_ci_upp)

mup <- plot_data %>%
  filter(str_detect(param,"^mu.p$"))

mupsi <- plot_data %>%
  filter(str_detect(param,"mu.psi"))

mubeta <- plot_data %>%
  filter(str_detect(param, "mu.beta")) %>%
  mutate(param = c("produ", "pheno", "struc"))

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

beta_list <- map(c("^beta1","^beta2","^beta3"),
                 beta_extract)

# Detection coefficients --------------------------------------------------
ggplot(alphas, aes(x = fct_reorder(param,mean), y = mean, ymin = lower, ymax = upper)) +
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
ggsave("data output/alphas.jpg")


# Occupancy coefficients --------------------------------------------------
ggplot(mubeta, aes(x = fct_reorder(param,mean), y = mean, ymin = lower, ymax = upper)) +
  geom_pointrange() +
  scale_y_continuous(limits = c(-1,1))+
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  ylab("Estimate (+- 95% CI)") + xlab("")+
  theme_classic()+
  labs(title = "Occupancy")+
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"))
ggsave("data output/mubeta.jpg")


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

ggsave("data output/occprob.jpg")


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

ggsave('data output/detecprob.jpg')


# Species richness --------------------------------------------------------
plot_data %>%
  filter(str_detect(param,"numspp")) %>%
  ggplot(aes(x = mean))+
  geom_histogram(binwidth = 0.5,
                 fill = "grey",
                 color = "black")+
  xlab("Species richness")

ggsave("data output/specrich.jpg")

# Random effect curve plots -----------------------------------------------

## Used to back-transform predictions
beta1_sc <- scale(station_data$produ_250) # beta 1
beta2_sc <- scale(station_data$pheno_250) # beta 2
beta3_sc <- scale(station_data$struc_250) # beta 3

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

  data <- tibble(mean = c(pred1,pred2,pred3),
                 betas = c(rep("a) Productivity", length(x)),
                           rep("b) Phenology", length(x)),
                           rep("c) Structure", length(x))
                           ),
                 spp = as.character(sppvec),
                 xdat = c(x1,x2,x3)
                 )
}

pred_data <- map_df(1:n_spp,
                 pred_data_extract)

## Check
pred_data %>%
  group_by(betas, spp) %>%
  tally()

pred_data %>%
  ggplot(aes(x = xdat, y = mean, color = spp))+
  geom_line() +
  scale_color_manual(values=c(rep("black",n_spp))) +
  theme(legend.position = "none") +
  xlab("")+
  ylab("Occupancy probability")+
  facet_wrap(~ betas, scales = "free_x")

ggsave('data output/Occ prob curves.jpg')

## Extracting community level predictions to add to plots above
## Not complete, still a work in progress...
mu.lpsi <- mupsi %>%
  pull(mean) %>%
  logit(.)

mubetas <- mubeta %>%
  pull(mean)

mu.lpsi_upp <- mupsi %>%
  pull(upper) %>%
  logit(.)

mubetas_upp <- mubeta %>%
  pull(upper)

mubetas_lower <- mubeta %>%
  pull(lower)

hist(pheno_vec)
hist(struc_vec)
hist(produ_vec)

x <- seq(-2,2,0.01)
pred1 <- plogis(mu.lpsi[1] + mubetas[1] * x)
pred2 <- plogis(mu.lpsi[1] + mubetas[2] * x)
pred3 <- plogis(mu.lpsi[1] + mubetas[3] * x)


# Write workspace ---------------------------------------------------------
save.image("data output/plotting_data.RData")

