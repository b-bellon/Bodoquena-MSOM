## ________________________________________________________________________

## Title:
## Purpose:
## Author:
## Date:

## Libraries
library(tidyverse)
library(nimble)
## ________________________________________________________________________


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

psi_spp <- plot_data %>%
  filter(str_detect(param,"lpsi")) %>%
  mutate_at(vars(mean:upper), boot::inv.logit) %>%
  mutate(spp = spp)

p_spp <- plot_data %>%
  filter(str_detect(param,"lp\\[")) %>%
  mutate_at(vars(mean:upper), boot::inv.logit) %>%
  mutate(spp = spp)


# Posterior plots ---------------------------------------------------------

## Detection coefficients
alphas <- plot_data %>%
  filter(str_detect(param,"alpha")) %>%
  mutate(param = c("Rain season"))

alphas

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

## Occupancy coefficients
betas <- plot_data %>%
  filter(str_detect(param,"beta")) %>%
  mutate(param = c("Open understory"))

ggplot(betas, aes(x = fct_reorder(param,mean), y = mean, ymin = lower, ymax = upper)) +
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
ggsave("data output/betas.jpg")

## Species occupancy probability
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

## Species detection probability
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

## Species richness
plot_data %>%
  filter(str_detect(param,"numspp")) %>%
  ggplot(aes(x = mean))+
  geom_histogram(binwidth = 0.5,
                 fill = "grey",
                 color = "black")+
  xlab("Species richness")

ggsave("data output/specrich.jpg")

# Write workspace ---------------------------------------------------------
save.image("data output/plotting_data.RData")

