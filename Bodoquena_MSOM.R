#R code to run the Multi-Species Occupancy Model
#The 1-day pooled observation data for mammals are stored in "Station_MammalSpp_Count_1d_pooled.xlsx"
#The camera trap station data are stored in "Station_data.xlsx"

rm(list=ls())

#Set working directory
setwd("C:/Users")

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

#Convert the data to a matrix format (Y) with the first dimension the stations (y-axis)
#and the second dimension the species (x-axis)
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

# #Survey covariate
# survey<- as.factor(as.vector(station.cov$SurveyName))
# survey.factor<- as.numeric(survey) 
# survey.levels<-length(levels(as.factor(survey)))

# #Habitat covariate
# habitat<- as.factor(as.vector(station.cov$Habitat))
# habitat.factor<- as.numeric(habitat)
# habitat.levels<- length(levels(as.factor(habitat)))
# 
# #Trail covariate (trail/off-trail_)
# trail<-station.cov$Trail
# trail<- as.factor(as.character(trail))
# trail.factor<- as.numeric(trail)
# trail.levels<- length(levels(as.factor(trail)))

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


# # Multi-session Royle-Nichols multi-species occupancy model with Cauchy random effect --------
# #This model is used to estimate abundance
#
# #Specify the model file, the file needs to be in the working directory
# modelFile='MSRoyleNicholsMSOccJAGS.txt'
# 
# #Specify the data
# occ.data = list(n=n, J=J, k=K, y=as.matrix(Y), survey.levels=survey.levels,survey.factor=survey.factor,habitat.levels=habitat.levels,habitat.factor=habitat.factor,trail.factor=trail.factor,trail.levels=trail.levels)
# 
# #Specify the parameters to be monitored
# occ.params = c('u', 'v','xi.v','xi.u','eta.survey.v','eta.survey.u','w','omega','z','v.trail','u.habitat','p.fit','p.fitnew')
# 
# #Specify the initial values
#     occ.inits = function() {
#     psi.meanGuess = runif(1, .25,1)
#     list(omega=runif(n,0.3,.9),w=matrix(1,nrow=n,ncol=survey.levels),
#                u=rnorm(n),
#                eta.survey.u=matrix(rnorm((n)*survey.levels),nrow=n,ncol=survey.levels),
#                v=rnorm(n),
#                eta.survey.v=matrix(rnorm((n)*survey.levels),nrow=n,ncol=survey.levels),
#                v.trail=matrix(c(rep(NA,n),rnorm((n)*(trail.levels-1))),nrow=n,ncol=trail.levels),
#                u.habitat=matrix(c(rep(NA,n),rnorm((n)*(habitat.levels-1))),nrow=n,ncol=habitat.levels),
#                mu.a = matrix(rbinom((n)*J, size=1, prob=1),
# 		           nrow=J, ncol=(n)),
# 		           sigma.u=runif(1,0,5),
# 		           sigma.v=runif(1,0,5)
#                )
#            }
# 
# 
# occ.inits()
# #run the model in JAGS using the R2jags package
# fit <- jags(occ.data, occ.inits, occ.params, modelFile,
#         n.chains=3, n.iter=50000, n.burnin=30000, n.thin=20)
# 
# #use the BUGS output format
# fit<-fit$BUGSoutput
# 
# #--- end of  RN model, jump to data analyses


# Multi-session multi-species occupancy model with Cauchy random effect from Tobler et al. 2015 --------

#Specify the model file, the file needs to be in the working directory
modelFile='MSMSOccJAGS.txt'

#Specify the data
occ.data = list(n=n, J=J, k=K, y=as.matrix(Y))
#occ.data = list(n=n, J=J, k=K, y=as.matrix(Y), season.levels=season.levels,season.factor=season.factor,understory.levels=understory.levels,understory.factor=understory.factor,LULC.factor=LULC.factor,LULC.levels=LULC.levels)

#Specify the parameters to be monitored
occ.params = c('u', 'v','xi.v','xi.u','eta.survey.v','eta.survey.u','w','omega','z','v.trail','u.habitat','p.fit','p.fitnew')

#Specify the initial values
    occ.inits = function() {
    psi.meanGuess = runif(1, .25,1)
    list(omega=runif(n,0.3,.9),w=matrix(1,nrow=n,ncol=survey.levels),
               u=rnorm(n),
               eta.survey.u=matrix(rnorm((n)*survey.levels),nrow=n,ncol=survey.levels),
               v=rnorm(n),
               eta.survey.v=matrix(rnorm((n)*survey.levels),nrow=n,ncol=survey.levels),
               v.trail=matrix(c(rep(NA,n),rnorm((n)*(trail.levels-1))),nrow=n,ncol=trail.levels),
               u.habitat=matrix(c(rep(NA,n),rnorm((n)*(habitat.levels-1))),nrow=n,ncol=habitat.levels),
               z = matrix(rbinom((n)*J, size=1, prob=1),
		           nrow=J, ncol=(n)),
		           sigma.u=runif(1,0,5),
		           sigma.v=runif(1,0,5)
               )
           }

occ.inits()

#run the model in JAGS using the R2jags package
fit <- jags(occ.data, occ.inits, occ.params, modelFile,
        n.chains=3, n.iter=500, n.burnin=100, n.thin=20)

#use the BUGS output format
fit<-fit$BUGSoutput


# Some results ------------------------------------------------------------

#Goodness of fit
dev.new(width=4.5, height=4.5)
op <- par(pty = "s",mai=c(0.9,0.8,0.2,0.1),cex.axis=1.1,cex.lab=1.2)
lim=c(0.9*min(c(fit$sims.list$p.fit, fit$sims.list$p.fitnew)),1.1*max(c(fit$sims.list$p.fit, fit$sims.list$p.fitnew)))
plot(fit$sims.list$p.fit, fit$sims.list$p.fitnew,main="",
     las=1, xlab=expression(paste(Chi^2,"for actual data set",sep="")), ylab=expression(Chi^2 ~" for ideal data sets"),xlim=lim,ylim=lim)
lines(c(0,max(lim)),c(0,max(lim)))
par(op)
mean(fit$sims.list$p.fit > fit$sims.list$p.fitnew) # Bayesian P-value
mean((fit$sims.list$p.fit) / (fit$sims.list$p.fitnew))  #lack-of-fit

#Probability of presence
w<-fit$sims.list$w
w.na<-w
w.na[w.na==0]<-NA
p.pres<-apply(w,c(2,3),mean)
rownames(p.pres)=sort(uspecies)
colnames(p.pres)<- levels(as.factor(survey))
round(p.pres,2)


#Estimated number of species for each survey
n.sp<-apply(w,c(1,3),sum)
colnames(n.sp)<-levels(as.factor(survey))
apply(n.sp,2,mean)


#Occupancy terra firme mean by survey
u<-fit$sims.list$u
eta.survey.u<-fit$sims.list$eta.survey.u
xi.u<-fit$sims.list$xi.u
u.habitat<-fit$sims.list$u.habitat
otf<-array(0,dim(eta.survey.u))
for(i in 1:dim(eta.survey.u)[1])
{
  for(j in 1:dim(eta.survey.u)[2])
  {
    otf[i,j,]<-u[i,j]+eta.survey.u[i,j,]*xi.u[i,j]+u.habitat[i,j,2]
  }
}
#occ.tf<-plogis(otf)*w   #simple model only
occ.tf<-exp(otf)*w   #RN model only
occ.tf<-1-exp(-occ.tf)   #RN model only
occ.tf<-apply(occ.tf,c(2,3),mean,na.rm=T)
rownames(occ.tf)=sort(uspecies)
colnames(occ.tf)<- levels(as.factor(survey))
round(occ.tf,3)


#Occupancy floodplain mean by survey
u<-fit$sims.list$u
eta.survey.u<-fit$sims.list$eta.survey.u
xi.u<-fit$sims.list$xi.u
ofp<-array(0,dim(eta.survey.u))
for(i in 1:dim(eta.survey.u)[1])
{
  for(j in 1:dim(eta.survey.u)[2])
  {
    ofp[i,j,]<-u[i,j]+eta.survey.u[i,j,]*xi.u[i,j]
  }
}
#occ.fp<-plogis(ofp)*w   #simple model only
occ.fp<-exp(ofp)*w #RN mmodel only
occ.fp<-1-exp(-occ.fp) #RN mmodel only
occ.fp<-apply(occ.fp,c(2,3),mean,na.rm=T)
rownames(occ.fp)=sort(uspecies)
colnames(occ.fp)<- levels(as.factor(survey))
round(occ.fp[,4:7],3)


#detection probabilities off-trail mean
v<-fit$sims.list$v
eta.survey.v<-fit$sims.list$eta.survey.v
xi.v<-fit$sims.list$xi.v
p<-array(0,dim(eta.survey.v))
for(i in 1:dim(eta.survey.v)[1])
{
  for(j in 1:dim(eta.survey.v)[2]) #species
  {
    p[i,j,]<-v[i,j]+eta.survey.v[i,j,]*xi.v[i,j]
  }
}
p.survey<-plogis(p)*w.na
p.survey<-apply(p.survey,c(2,3),mean,na.rm=T)
rownames(p.survey)=sort(uspecies)
colnames(p.survey)<- levels(as.factor(survey))
round(p.survey,3)


#detection probabilities trail mean
v<-fit$sims.list$v
eta.survey.v<-fit$sims.list$eta.survey.v
xi.v<-fit$sims.list$xi.v
v.trail<-fit$sims.list$v.trail
p<-array(0,dim(eta.survey.v))
for(i in 1:dim(eta.survey.v)[1])
{
  for(j in 1:dim(eta.survey.v)[2]) #species
  {
    p[i,j,]<-v[i,j]+eta.survey.v[i,j,]*xi.v[i,j] + v.trail[i,j,2]
  }
}
p.survey<-plogis(p)*w.na
p.survey<-apply(p.survey,c(2,3),mean,na.rm=T)
rownames(p.survey)=sort(uspecies)
colnames(p.survey)<- levels(as.factor(survey))
round(p.survey,3)


