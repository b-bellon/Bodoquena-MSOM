# Bodoquena-MSOM
# Bodoquena-MSOM
#R code to run the multi-species occupancy model (29 mammal spp.)
#The 1-day pooled observation data for mammals are stored in "Station_MammalSpp_Count_1d_pooled.xlsx"
#The camera trap station data is stored in "Station_data.xlsx"

#Change this to the directory where your files are
setwd("C:/Users/")

#Load the libraries
#install.packages("R2jags")
#install.packages("reshape")
#install.packages("readxl")
library(R2jags)
library(reshape)
library(readxl)

#define the data files
occfile<-"Station_MammalSpp_Count_1d_pooled.xlsx"
stationfile<-"Station_data.xlsx"

#Load the observation data
data <- read_excel(occfile, sheet = "Sheet1")
head(data)

#List the unique species
uspecies = as.character(unique(data$Species))

#Number of observed species
n=length(uspecies)

#List the  camera stations
ustations = as.character(unique(data$ID_Station))

#Number of camera stations
J=length(ustations)

#Convert the data to a matrix format (Y) with the first dimension the stations (y-axis) and the second dimension the species (x-axis)
melt.tmp=melt(data,id.var=c("Species", "ID_Station"), measure.var="Count")
Y=cast(melt.tmp, ID_Station ~ Species,sum)
Y<-Y[,-1]
Y
dim(Y)

#Load the station covariates
station.cov <- read_excel(stationfile, sheet = "Sheet1")
head(station.cov)

#Number of days each camera was operating
K=station.cov$Days #if 6-day pooled data: K=round(station.cov$Days/6,0)

#Season covariate (Rainy/Dry)
season <- as.factor(as.vector(station.cov$Season))
season.factor<- as.numeric(season) #1 = Dry; 2 = Rainy
season.levels<-length(levels(season))

#Understory covariate (Dense/Open)
understory <- as.factor(as.vector(station.cov$Understory))
understory.factor<- as.numeric(understory) #1 = Dense; 2 = Open
understory.levels<-length(levels(understory))

#LULC covariate (Mapbiomas v5 class)
LULC <- as.factor(as.numeric(station.cov$MapBiomas))
LULC.factor<- as.numeric(LULC) #1 = Forest; 2 = Savanna; 3 = Grassland; 4 = Pasture; 5 = Soybean
LULC.levels<-length(levels(LULC))

