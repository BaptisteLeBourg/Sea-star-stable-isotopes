# This script focuses on the investigation of the impacts of the the trophic group,
# the environmental variables and their interactions on the stable isotope values 
# of sea stars in the Southern ocean with linear models and ANCOVAs.

################################################################################
# Packages
################################################################################

# Let's start by loading the necessary packages.
library(dplyr) # For dataset combination
library(car) # Contains the Anova function
library(agricolae) # To do the Scheffe tests.

################################################################################
# Data loading and preparation
################################################################################

# Let's load and prepare the stable isotope and environmental data.
DataSeaStar<-read.table("IsotopeData.csv",dec=".",sep=",",header=T)
DataEnv<-read.table("EnvironmentalData.csv",dec=".",sep=",",header=T)

# First we combine both datasets.
Data <- DataSeaStar %>% left_join(DataEnv, by = c("StationID","ExpeditionID","Date"))

# We need to log-transform depth (log 10) and chlorophyll concentration (natural
# log) data prior to doing the analysis. We also need to logit-transform sea ice
# concentration.
Data$logDepth<-log10(Data$Depth)
Data$logChl<-log(Data$chl_prev_month)
Data$logIce<-logit(Data$seaice_prev_month)

# We also need to specify that the trophic group is a factor.
Data$Trophic_group<-as.factor(Data$Trophic_group)

# For models with sea ice concentration and sea ice duration as explanatory 
# variables, we need a dataframe where individuals from Subantarctic regions are
# removed.  These individuals were collected in areas where sea ice has not been
# present since more than 1000 days.
DataAnt <- Data[-which(Data$days_since_melt > 1000),]

# For the models including trophic group, log-transformed depth, logit-transformed
# sea ice concentration, sea ice season and log-transformed chlorophyll as explanatory
# variables, it was necessary to create a dataset removing the data from the 
# predators of encrusting prey as only four individuals were available for this 
# trophic group in this model.
DataAntChl<-DataAnt[!is.na(DataAnt$chl_prev_month),]
which(DataAntChl$Trophic_group=="Pelagic-Carnivore2")
DataAntChl<-DataAntChl[-c(16,378,388,498),]

################################################################################
# Linear models
################################################################################

# Decreasing δ13C values in particulate organic matter from Subantarctic to 
# Antarctic areas and high δ13C values in the sea ice microbial community in 
# Antarctic may induce non-linear variations of δ13C values with sea ice 
# concentration or ice season duration if data from Subantarctic stations are 
# included in the ANCOVAs as no sea ice is present in Subantarctic areas. 
# Furthermore, chlorophyll concentration data were not equally available for all
# sampling stations. Consequently, the models and ANCOVAs were performed four times.

# First model: includes both Subantarctic and Antarctic sea stars. Includes trophic 
# group and log-transformed depth and their interactions as explanatory variables.
# A post-hoc Scheffé test was also performed for this model to further assess the
# effect of the trophic group on δ13C and δ15N values.
#####  Model for δ13C values.
LM13C<-lm(d13C~Trophic_group+logDepth
          +Trophic_group:logDepth,data=Data,contrasts=list(Trophic_group="contr.poly"))
Anova(LM13C,type="III")
qqPlot(LM13C)

posthoc<-scheffe.test(LM13C,"Trophic_group",group=T)
posthoc

#####  Model for δ15N values.
LM15N<-lm(d15N~Trophic_group+logDepth
          +Trophic_group:logDepth,data=Data,contrasts=list(Trophic_group="contr.poly"))
Anova(LM15N,type="III")
qqPlot(LM15N)

posthoc<-scheffe.test(LM15N,"Trophic_group",group=T)
posthoc

# Second model: includes Subantarctic and Antarctic sea stars. Includes trophic 
# group, log-transformed depth and log-transformed chlorophyll concentration
# and their interactions as explanatory variables.
#####  Model for δ13C values.
LM13C<-lm(d13C~Trophic_group+logDepth+logChl
          +Trophic_group:logDepth
          +Trophic_group:logChl
          +logDepth:logChl,data=Data,contrasts=list(Trophic_group="contr.poly"))
Anova(LM13C,type="III")
qqPlot(LM13C)

#####  Model for δ15N values.
LM15N<-lm(d15N~Trophic_group+logDepth+logChl
          +Trophic_group:logDepth
          +Trophic_group:logChl
          +logDepth:logChl,data=Data,contrasts=list(Trophic_group="contr.poly"))
Anova(LM15N,type="III")
qqPlot(LM15N)

# Third model: includes Antarctic sea stars only. Includes trophic group, log-
# transformed depth, logit-transformed sea ice concentration and sea ice season
# and their interactions as explanatory variables.
#####  Model for δ13C values.
LM13C<-lm(d13C~logIce+logDepth+seaice_last_730
          +Trophic_group:logIce
          +Trophic_group:logDepth
          +Trophic_group:seaice_last_730
          +logDepth:logIce
          +logDepth:seaice_last_730
          +logIce:seaice_last_730,data=DataAnt,contrasts=list(Trophic_group="contr.poly"))
Anova(LM13C,type="III")
qqPlot(LM13C)

#####  Model for δ15N values.
LM15N<-lm(d15N~Trophic_group+logIce+logDepth+seaice_last_730
          +Trophic_group:logIce
          +Trophic_group:logDepth
          +Trophic_group:seaice_last_730
          +logDepth:logIce
          +logDepth:seaice_last_730
          +logIce:seaice_last_730,data=DataAnt,contrasts=list(Trophic_group="contr.poly"))
Anova(LM15N,type="III")
qqPlot(LM15N)

# Fourth model: includes Antarctic sea star only and excludes the predators of 
# encrusting prey. Includes trophic group, log-transformed depth, logit-transformed
# sea ice concentration, sea ice season and log-transformed chlorophyll and their 
# interactions as explanatory variables.
#####  Model for δ13C values.
LM13C<-lm(d13C~Trophic_group+logDepth+logIce+seaice_last_730+logChl
          +Trophic_group:logDepth
          +Trophic_group:logIce
          +Trophic_group:seaice_last_730
          +Trophic_group:logChl
          +logIce:logDepth
          +logDepth:seaice_last_730
          +logDepth:logChl
          +logIce:seaice_last_730
          +logIce:logChl
          +seaice_last_730:logChl,data=DataAntChl
          ,contrasts=list(Trophic_group="contr.poly"))
Anova(LM13C,type="III")
qqPlot(LM13C)

#####  Model for δ15N values.
LM15N<-lm(d15N~Trophic_group+logDepth+logIce+seaice_last_730+logChl
          +Trophic_group:logDepth
          +Trophic_group:logIce
          +Trophic_group:seaice_last_730
          +Trophic_group:logChl
          +logIce:logDepth
          +logDepth:seaice_last_730
          +logDepth:logChl
          +logIce:seaice_last_730
          +logIce:logChl
          +seaice_last_730:logChl,data=DataAntChl
          ,contrasts=list(Trophic_group="contr.poly"))
Anova(LM15N,type="III")
qqPlot(LM15N)

################################################################################
# Effects of the interactions between trophic groups and environmental variables
################################################################################

# The influence of the interactions between the trophic group and environmental 
# variables on stable isotope values were further assessed with Pearson correlation 
# tests, with P-values being adjusted with the Bonferroni method 

# interactions between trophic groups and depth…
ngroups <- length(unique(Data$Trophic_group))
spx <- split(Data$logDepth,Data$Trophic_group)

##### …for δ13C values.
spy <- split(Data$d13C,Data$Trophic_group)
correlDepth <- data.frame(matrix(ncol = 3, nrow = ngroups))
for (i in c(1:ngroups)){
  o<-cor.test(spx[[i]],spy[[i]])
  correlDepth[i,1]<-levels(Data$Trophic_group)[i]
  correlDepth[i,2]<-o$estimate
  correlDepth[i,3]<-o$p.value
  correlDepth[i,4]<-p.adjust(correlDepth[i,3], method = "bonferroni",n=8)
}
colnames(correlDepth)<-c("Trophic group","r","P","Adjusted P")
correlDepth[c(1,2,5,6,8,3,4,7),]

##### …for δ15N values.
spy <- split(Data$d15N,Data$Trophic_group)
correlDepth <- data.frame(matrix(ncol = 3, nrow = ngroups))
for (i in c(1:ngroups)){
  o<-cor.test(spx[[i]],spy[[i]])
  correlDepth[i,1]<-levels(Data$Trophic_group)[i]
  correlDepth[i,2]<-o$estimate
  correlDepth[i,3]<-o$p.value
  correlDepth[i,4]<-p.adjust(correlDepth[i,3], method = "bonferroni",n=8)
}
colnames(correlDepth)<-c("Trophic group","r","P","Adjusted P")
correlDepth[c(1,2,5,6,8,3,4,7),]


# interactions between trophic groups and sea ice concentration
ngroups <- length(unique(DataAnt$Trophic_group))
spx <- split(DataAnt$logIce,DataAnt$Trophic_group)

##### …for δ13C values.
spy <- split(DataAnt$d13C,DataAnt$Trophic_group)
correlIce1 <- data.frame(matrix(ncol = 3, nrow = ngroups))
for (i in c(1:ngroups)){
  o<-cor.test(spx[[i]],spy[[i]])
  correlIce1[i,1]<-levels(Data$Trophic_group)[i]
  correlIce1[i,2]<-o$estimate
  correlIce1[i,3]<-o$p.value
  correlIce1[i,4]<-p.adjust(correlIce1[i,3], method = "bonferroni",n=8)
}
colnames(correlIce1)<-c("Trophic group","r","P","Adjusted P")
correlIce1[c(1,2,5,6,8,3,4,7),]

##### …for δ15N values.
spy <- split(DataAnt$d15N,DataAnt$Trophic_group)
correlIce1 <- data.frame(matrix(ncol = 3, nrow = ngroups))
for (i in c(1:ngroups)){
  o<-cor.test(spx[[i]],spy[[i]])
  correlIce1[i,1]<-levels(Data$Trophic_group)[i]
  correlIce1[i,2]<-o$estimate
  correlIce1[i,3]<-o$p.value
  correlIce1[i,4]<-p.adjust(correlIce1[i,3], method = "bonferroni",n=8)
}
colnames(correlIce1)<-c("Trophic group","r","P","Adjusted P")
correlIce1[c(1,2,5,6,8,3,4,7),]

# interactions between trophic groups and sea ice duration
ngroups <- length(unique(DataAnt$Trophic_group))
spx <- split(DataAnt$seaice_last_730,DataAnt$Trophic_group)

##### …for δ13C values.
spy <- split(DataAnt$d13C,DataAnt$Trophic_group)
correlIce2 <- data.frame(matrix(ncol = 3, nrow = ngroups))
for (i in c(1:ngroups)){
  o<-cor.test(spx[[i]],spy[[i]])
  correlIce2[i,1]<-levels(Data$Trophic_group)[i]
  correlIce2[i,2]<-o$estimate
  correlIce2[i,3]<-o$p.value
  correlIce2[i,4]<-p.adjust(correlIce2[i,3], method = "bonferroni",n=8)
}
colnames(correlIce2)<-c("Trophic group","r","P","Adjusted P")
correlIce2[c(1,2,5,6,8,3,4,7),]

##### …for δ15N values.
spy <- split(DataAnt$d15N,DataAnt$Trophic_group)
correlIce2 <- data.frame(matrix(ncol = 3, nrow = ngroups))
for (i in c(1:ngroups)){
  o<-cor.test(spx[[i]],spy[[i]])
  correlIce2[i,1]<-levels(Data$Trophic_group)[i]
  correlIce2[i,2]<-o$estimate
  correlIce2[i,3]<-o$p.value
  correlIce2[i,4]<-p.adjust(correlIce2[i,3], method = "bonferroni",n=8)
}
colnames(correlIce2)<-c("Trophic group","r","P","Adjusted P")
correlIce2[c(1,2,5,6,8,3,4,7),]

# interactions between trophic groups and chlorophyll concentration
ngroups <- length(unique(Data$Trophic_group))
spx <- split(Data$logChl,Data$Trophic_group)

##### …for δ13C values.
spy <- split(Data$d13C,Data$Trophic_group)
correlChl <- data.frame(matrix(ncol = 3, nrow = ngroups))
for (i in c(1:ngroups)){
  o<-cor.test(spx[[i]],spy[[i]])
  correlChl[i,1]<-levels(Data$Trophic_group)[i]
  correlChl[i,2]<-o$estimate
  correlChl[i,3]<-o$p.value
  correlChl[i,4]<-p.adjust(correlChl[i,3], method = "bonferroni",n=8)
}
colnames(correlChl)<-c("Trophic group","r","P","Adjusted P")
correlChl[c(1,2,5,6,8,3,4,7),]

##### …for δ15N values.
spy <- split(Data$d15N,Data$Trophic_group)
correlChl <- data.frame(matrix(ncol = 3, nrow = ngroups))
for (i in c(1:ngroups)){
  o<-cor.test(spx[[i]],spy[[i]])
  correlChl[i,1]<-levels(Data$Trophic_group)[i]
  correlChl[i,2]<-o$estimate
  correlChl[i,3]<-o$p.value
  correlChl[i,4]<-p.adjust(correlChl[i,3], method = "bonferroni",n=8)
}
colnames(correlChl)<-c("Trophic group","r","P","Adjusted P")
correlChl[c(1,2,5,6,8,3,4,7),]

# End of script.