## Libraries
library(tidyverse)
## ________________________________________________________________________

## renv::init() - run this once project code is ready for first commit
## renv::snapshot() - run when new packages are installed (to update renv.lock file)

dnd <- readxl::read_xlsx("data input/Station_MammalSpp_Count_1d_pooled.xlsx") %>%
  janitor::clean_names()

dnd

spp <- dnd %>%
  distinct(species) %>%
  pull

Y <- dnd %>%
  pivot_wider(names_from = c(species),
              values_from = count) %>%
  select(-id_station) %>%
  as.matrix()

Y
dim(Y)
dimnames(Y) <- list(NULL, spp)
dimnames(Y)
Y[1:5,1:10]
n_sites <- dim(Y)[1] # number of sites
n_spp <- dim(Y)[2] # number of species

n_spp
n_sites

freq_data <- dnd %>%
  group_by(species) %>%
  summarise(n_detections = sum(count)) %>%
  print(n = 40)

freq_data %>%
  ggplot(aes(x = species, y = n_detections))+
  geom_bar(stat = "identity")+
  coord_flip()

ggsave("data output/detec_freq.jpg")

station_data <- readxl::read_xlsx("data input/Station_data.xlsx") %>%
  janitor::clean_names()

station_data
summary(station_data)

J <-  station_data$days
J;length(J)

table(station_data$season)
table(station_data$understory)

season_vec <- as.numeric(factor(station_data$season)) # Dry = 1, Rain = 2
under_vec <- as.numeric(factor(station_data$understory)) # Dense = 1, Open = 2

produ_vec <- as.numeric(scale(station_data$produ_250)[,1])
pheno_vec <- as.numeric(scale(station_data$pheno_250)[,1])
struc_vec <- as.numeric(scale(station_data$struc_250)[,1])

save.image("data output/model_data.RData")
