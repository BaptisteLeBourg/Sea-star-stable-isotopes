# This script focuses on the investigation of the isotopic niches of sea stars in 
# the Southern Ocean, depending of the trophic groups and of the environmental 
# variables.

################################################################################
# Packages
################################################################################

# Let's start by loading the necessary packages.
library(dplyr) # For dataset combination.
library(plyr)
library(SIBER) # Contains the functions necessary for the investigation of isotopic niches.
library(agricolae) # Contains the function orderPvalue.
library(Rcmdr) # Contains the logit function
library(plotrix) # Contains the ablineclip function

################################################################################
# Functions
################################################################################

# The first function is groupMetricsML2. This function is a modification of the 
# GroupMetricMl function in SIBER. It provides most of the standard ellipse metrics
# (Standard Ellipse Area SEA and its corresponding small sample size corrected 
# version SEAc, eccentricity, angle in radians), but also Pseudo-standard deviation 
# for the x (PSDx) and y axes (PSDy), and their corresponding small sample size 
# corrected version (PSDxc and PSDyc).
groupMetricsML2<-function (siber) 
{
  tmp.names <- NULL
  for (i in 1:siber$n.communities) {
    tmp.names <- c(tmp.names, paste(siber$all.communities[i], 
                                    siber$group.names[[i]], sep = "."))
  }
  out <- matrix(NA, nrow = 8, ncol = sum(siber$n.groups[2, 
  ]), dimnames = list(c("SEA", "SEAc","Angle","Eccentricity","PSDx","PSDy","PSDxc","PSDyc"), tmp.names))
  cnt <- 1
  for (i in 1:siber$n.communities) {
    for (j in 1:siber$n.groups[2, i]) {
      tmp.SEA <- sigmaSEA(siber$ML.cov[[i]][, , j])
      out["SEA", cnt] <- tmp.SEA$SEA
      n <- siber$sample.sizes[i, paste(siber$group.names[[i]][j])]
      out["SEAc", cnt] <- tmp.SEA$SEA * (n - 1)/(n - 2)
      out["Eccentricity", cnt] <- tmp.SEA$eccentricity
      out["Angle", cnt] <- tmp.SEA$theta
      out["PSDx", cnt] <- cos(abs(tmp.SEA$theta))*2*tmp.SEA$a
      out["PSDy", cnt] <- sin(abs(tmp.SEA$theta))*2*tmp.SEA$a
      out["PSDxc", cnt] <- cos(abs(tmp.SEA$theta))*2*tmp.SEA$a*(sqrt((n-1)/(n-2)))
      out["PSDyc", cnt] <- sin(abs(tmp.SEA$theta))*2*tmp.SEA$a*(sqrt((n-1)/(n-2)))
      idx <- siber$raw.data[[i]]$group == siber$group.names[[i]][j]
      cnt <- cnt + 1
    }
  }
  return(out)
}

# The second function is posteriorPSDx. This function is a modification of the 
# posteriorSEA function in SIBER. This function loops over each posterior draw 
# of a single group's Bayesian bivariate ellipse and calculates the Pseudo-
# standard deviation on the x axis (PSDx) for each draw, thereby generating a 
# distribution of PSDx estimates. Like posteriorSEA, it is not intended for direct
# calling outside of siberPSDx.
posteriorPSDx <- function (post) {
  Nobs <- nrow(post)
  cosTHETA.B <- numeric(Nobs) 
  for (i in 1:Nobs) {
    estS <- post[i, 1:4]
    dim(estS) <- c(2, 2)
    cosTHETA.B[i] <- cos(abs(sigmaSEA(estS)$theta))*2*sigmaSEA(estS)$a
  } 
  return(cosTHETA.B)
}

# The third function is siberPSDx. It is a modification of the siberEllipse 
# function in SIBER. This function loops over each group within each community 
# and calculates the posterior distribution describing the corresponding 
# Pseudo-Standard Deviation on the x axis.
siberPSDx <- function (corrected.posteriors) {
  cosTHETA.B <- matrix(NA, 
                       nrow = nrow(corrected.posteriors[[1]]),
                       ncol = length(corrected.posteriors))
  for (i in 1:length(corrected.posteriors)){
    tmp <- posteriorPSDx(corrected.posteriors[[i]])
    cosTHETA.B[, i] <- tmp
  }
  return(cosTHETA.B)
}

# The fourth function is posteriorPSDy. it is a modification of the posteriorSEA 
# function in SIBER. This function loops over each posterior draw of a single 
# group's Bayesian bivariate ellipse and calculates the Pseudo-standard deviation
# on the y axis (PSDy) for each draw, thereby generating a distribution of PSDx 
# estimates. Like posteriorSEA, it is not intended for direct calling outside of
# siberPSDy.
posteriorPSDy <- function (post) {
  Nobs <- nrow(post)
  sinTHETA.B <- numeric(Nobs) 
  for (i in 1:Nobs) {
    estS <- post[i, 1:4]
    dim(estS) <- c(2, 2)
    sinTHETA.B[i] <- sin(abs(sigmaSEA(estS)$theta))*2*sigmaSEA(estS)$a
  } 
  return(sinTHETA.B)
}

# The fourth function is siberPSDy. It is a modification of the siberEllipse 
# function in SIBER. This function loops over each group within each community 
# and calculates the posterior distribution describing the corresponding Pseudo-
# Standard Deviation on the y axis.
siberPSDy <- function (corrected.posteriors) {
  sinTHETA.B <- matrix(NA, 
                       nrow = nrow(corrected.posteriors[[1]]),
                       ncol = length(corrected.posteriors))
  for (i in 1:length(corrected.posteriors)){
    tmp <- posteriorPSDy(corrected.posteriors[[i]])
    sinTHETA.B[, i] <- tmp
  }
  return(sinTHETA.B)
}

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

# For plots displaying relationship of stable isotope values with sea ice 
# concentration and sea ice duration, we need a dataframe where individuals from
# Subantarctic regions are removed. These individuals were collected in areas 
# where sea ice has not been present since more than 1000 days.
DataAnt <- Data[-which(Data$days_since_melt > 1000),]

################################################################################
# Mean correction of stable isotope for each station
################################################################################
# The ellipses will be generated at the scale of the whole Southern Ocean. Yet,
# stable isotope values in sea stars may differ between sampling stations because
# of spatial variations of stable isotope values of the primary food sources at 
# the baseline of the food webs, or of the potential influence of other 
# environmental parameters than the one of interest on stable isotope values. 
# To avoid this issue, we mean-corrected stable isotope values for each station.

##### Mean-correction for δ13C values.
Data<-ddply(Data,c("StationID"),transform,md13Cst = mean(d13C)) # mean of δ13C values for each station.
Data$md13Ctot = mean(Data$d13C) # mean of δ13C values for the whole dataset.
Data$d13Ccorr<-Data$d13C-(Data$md13Cst-Data$md13Ctot) # mean-correction of δ13C values.

##### Mean-correction for δ15N values.
Data<-ddply(Data,c("StationID"),transform,md15Nst = mean(d15N,na.rm=T)) # mean of δ15N values for each station.
Data$md15Ntot = mean(Data$d15N,na.rm=T) # mean of δ15N values for the whole dataset.
Data$d15Ncorr<-Data$d15N-(Data$md15Nst-Data$md15Ntot) # mean-correction of δ15N values.

################################################################################
# Isotopic niches for each trophic group
################################################################################

# We first select the columns necessary to generate the standard ellipses, i.e.,
# the mean-corrected δ13C and δ15N values, the species, and the trophic group as
# the community.
nicheSIBER<-Data[,c(23,26,3,11)]

# We remove suspension feeders from this dataset because only three taxa were 
# available for this trophic group.
nicheSIBER<-droplevels(nicheSIBER[-which(nicheSIBER$Trophic_group=="Pelagic-Suspension-feeder"),])

# We give the required names to the columns.
colnames(nicheSIBER) = c("iso1", "iso2", "group","community")

# For each trophic group, one stable isotope value has to be the mean stable 
# isotope value for one taxa, so we compute this mean for both δ13C and δ15N 
# values.
meand13C<-aggregate(iso1~group+community, data=nicheSIBER, FUN="mean")
meand15N<-aggregate(iso2~group+community, data=nicheSIBER, FUN="mean")
nicheSIBER<-left_join(meand13C,meand15N)

# We can now consider all taxa as a single species to ensure that the standard
# ellipses will be generated at the scale of the trophic groups only.
nicheSIBER$group<-rep("AllSpecies",nrow(nicheSIBER))

# Re-ordering of the columns in the right order and creation of the SIBER object.
nicheSIBER<-nicheSIBER[,c(3,4,1,2)]
siber.Test <- createSiberObject(nicheSIBER)

# Computation of PSDxc and PSDyc for each trophic group with our custom function
# groupMetricsML2.
group.ML <- groupMetricsML2(siber.Test)
print(group.ML)

# We can now conduct the Bayesian estimations of PSDx and PSDy.

# Options for running jags
parms <- list()
parms$n.iter <- 10^5   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains
parms$save.output = TRUE
parms$save.dir = tempdir()

# Definition of the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# Fitting of the ellipses.
ellipses.posterior <- siberMVN(siber.Test, parms, priors)

# Bayesian estimation of PSDx with our custom function siberPSDx.
PSDx.B <- siberPSDx(ellipses.posterior)
colnames(PSDx.B)<-colnames(group.ML)

# Estimation of the credible intervals.
PSDx.B.credibles <- lapply(
  as.data.frame(PSDx.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = 0.95)

# Estimation of the modes.
PSDx.B.modes <- lapply(
  as.data.frame(PSDx.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = 0.95, all.modes=T)

# Comparison of PSDx for each pair of ellipse in a matrix.
Q <- matrix(1, ncol = ncol(PSDx.B), nrow = ncol(PSDx.B))
for (i in 1:(ncol(PSDx.B) - 1)) {
  for (j in (i + 1):ncol(PSDx.B)) {
    Q[i, j] <- if(sum( PSDx.B[,i] < PSDx.B[,j] ) / nrow(PSDx.B) < 0.5) {
      sum( PSDx.B[,i] < PSDx.B[,j] ) / nrow(PSDx.B)}
    else{sum( PSDx.B[,i] > PSDx.B[,j] ) / nrow(PSDx.B)}
    Q[j, i] <- if(sum( PSDx.B[,i] < PSDx.B[,j] ) / nrow(PSDx.B) < 0.5) {
      sum( PSDx.B[,i] < PSDx.B[,j] ) / nrow(PSDx.B)}
    else{sum( PSDx.B[,i] > PSDx.B[,j] ) / nrow(PSDx.B)}
  }
}

# Sorting of ellipses in groups with similar PSDx.
mode <- data.frame(matrix(unlist(PSDx.B.modes), nrow=length(PSDx.B.modes), byrow=T))
colnames(mode)<-"mode"
mode$group<-colnames(PSDx.B)
groupsPSDx<- orderPvalue(mode[, 2], mode[, 1], 0.05, Q)
groupsPSDx

# We apply the same process for the Bayesian estimation of PSDy with our custom 
# function siberPSDy.
PSDy.B <- siberPSDy(ellipses.posterior)
colnames(PSDy.B)<-colnames(group.ML)

# Estimation of the credible intervals.
PSDy.B.credibles <- lapply(
  as.data.frame(PSDy.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = 0.95)

# Estimation of the modes.
PSDy.B.modes <- lapply(
  as.data.frame(PSDy.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = 0.95, all.modes=T)

# Comparison of PSDy for each pair of ellipse in a matrix.
Q <- matrix(1, ncol = ncol(PSDy.B), nrow = ncol(PSDy.B))
for (i in 1:(ncol(PSDy.B) - 1)) {
  for (j in (i + 1):ncol(PSDy.B)) {
    Q[i, j] <- if(sum( PSDy.B[,i] < PSDy.B[,j] ) / nrow(PSDy.B) < 0.5) {
      sum( PSDy.B[,i] < PSDy.B[,j] ) / nrow(PSDy.B)}
    else{sum( PSDy.B[,i] > PSDy.B[,j] ) / nrow(PSDy.B)}
    Q[j, i] <- if(sum( PSDy.B[,i] < PSDy.B[,j] ) / nrow(PSDy.B) < 0.5) {
      sum( PSDy.B[,i] < PSDy.B[,j] ) / nrow(PSDy.B)}
    else{sum( PSDy.B[,i] > PSDy.B[,j] ) / nrow(PSDy.B)}
  }
}

# Sorting of ellipses in groups with similar PSDy.
mode <- data.frame(matrix(unlist(PSDy.B.modes), nrow=length(PSDy.B.modes), byrow=T))
colnames(mode)<-"mode"
mode$group<-colnames(PSDy.B)
groupsPSDy<- orderPvalue(mode[, 2], mode[, 1], 0.05, Q)
groupsPSDy

# Summary of the results.

# Preparation of data for the plots and a table summarising the results of
# PSDx and PSDy estimation. We need to reorder the trophic groups in the order we
# want to see them on the table and the plots.
PSDxTrophicGroup<-PSDx.B[,c(1,2,5,6,3,4,7:ncol(PSDx.B)) ]
modePSDxTrophicGroup<-as.data.frame(groupsPSDx)
QPSDxTrophicGroup<-PSDx.B.credibles[c(1,2,5,6,3,4,7)]
modePSDxTrophicGroup$group<-rownames(modePSDxTrophicGroup)
modePSDxTrophicGroup<-modePSDxTrophicGroup[with(modePSDxTrophicGroup, order(group)),]
modePSDxTrophicGroup<-modePSDxTrophicGroup[c(1,2,5,6,3,4,7),]

PSDyTrophicGroup<-PSDy.B[,c(1,2,5,6,3,4,7:ncol(PSDy.B)) ]
QPSDyTrophicGroup<-PSDy.B.credibles[c(1,2,5,6,3,4,7)]
modePSDyTrophicGroup<-as.data.frame(groupsPSDy)
modePSDyTrophicGroup$group<-rownames(modePSDyTrophicGroup)
modePSDyTrophicGroup<-modePSDyTrophicGroup[with(modePSDyTrophicGroup, order(group)),]
modePSDyTrophicGroup<-modePSDyTrophicGroup[c(1,2,5,6,3,4,7),]

groupMLTrophicGroup<-group.ML[,c(1,2,5,6,3,4,7)]

# Creation of a table summarising the results of PSDx and PSDy estimation.
tableTrophicGroup<-data.frame(matrix(nrow=nrow(modePSDyTrophicGroup), ncol=3))
colnames(tableTrophicGroup)<-c("Trophic Group","PSDCB","PSDNB")

for (i in 1:nrow(tableTrophicGroup)){
  tableTrophicGroup[i,1]<-modePSDxTrophicGroup[i,3]
  tableTrophicGroup[i,2]<-paste(round(modePSDxTrophicGroup[i,1],1),round(as.data.frame(QPSDxTrophicGroup[i])[2,1],1),sep=" (")
  tableTrophicGroup[i,2]<-paste(tableTrophicGroup[i,2],round(as.data.frame(QPSDxTrophicGroup[i])[2,2],1),sep="-")
  tableTrophicGroup[i,2]<-paste(tableTrophicGroup[i,2],modePSDxTrophicGroup[i,2],sep=") ")
  tableTrophicGroup[i,3]<-paste(round(modePSDyTrophicGroup[i,1],1),round(as.data.frame(QPSDyTrophicGroup[i])[2,1],1),sep=" (")
  tableTrophicGroup[i,3]<-paste(tableTrophicGroup[i,3],round(as.data.frame(QPSDyTrophicGroup[i])[2,2],1),sep="-")
  tableTrophicGroup[i,3]<-paste(tableTrophicGroup[i,3],modePSDyTrophicGroup[i,2],sep=") ")
}

tableTrophicGroup

# Creation of plots summarizing the results.

# Preparation of a data frame for the plot displaying mean ± SD stable isotope 
# values for each trophic group.
mean13CTrophGr<-aggregate(d13C~Trophic_group,data=Data,FUN="mean")
sd13CTrophGr<-aggregate(d13C~Trophic_group,data=Data,FUN="sd")
colnames(sd13CTrophGr)<-c("Trophic_group" ,"sd13C")
mean15NTrophGr<-aggregate(d15N~Trophic_group,data=Data,FUN="mean")
sd15NTrophGr<-aggregate(d15N~Trophic_group,data=Data,FUN="sd")
colnames(sd15NTrophGr)<-c("Trophic_group" ,"sd15N")
sdTrophGr<-left_join(sd13CTrophGr,sd15NTrophGr)
meanTrophGr<-left_join(mean13CTrophGr,mean15NTrophGr)
meansdTrophGr<-left_join(meanTrophGr,sdTrophGr)

# Window settings.
dev.new(width=17, height=3*(5.8-0.4461538), unit="cm")
layout(matrix(c(
  6,2,
  6,2,
  1,2,
  1,2,
  1,2,
  1,3,
  1,3,
  1,3,
  5,3,
  5,3,
  5,4,
  5,4,
  5,4
), nrow=12, ncol=2, byrow = TRUE))
par(mar=c(1.5,12,0.5,0.3))

# Plot displaying mean ± SD stable isotope values for each trophic group.

plot(meansdTrophGr$d15N~meansdTrophGr$d13C,type="n",
     xlab="",ylab="",font.lab=2,cex.lab=2,cex=1,cex.axis=2,ylim=c(4,24),
     xlim=c(-26,-6),ax=F,font.main=2,cex.main=2)

axis(side=2, at=c(4,9,14,19,24), labels=c("4.0","9.0","14.0","19.0","24.0"),lwd=2,cex.axis=3.5,font=2, las=2)
axis(side=1, at=c(-26,-21,-16,-11,-6), labels=rep("",5),lwd=2,cex.axis=3.5,font=2)
axis(side=1, at=c(-26,-21,-16,-11,-6), labels=c("-26.0","-21.0","-16.0","-11.0","-6.0"),lwd=0,cex.axis=3.5,font=2,line=1)
box(lwd=2)

arrows(meansdTrophGr$d13C-meansdTrophGr$sd13C, meansdTrophGr$d15N, meansdTrophGr$d13C+ meansdTrophGr$sd13C, meansdTrophGr$d15N, length=0.05,lwd=2, angle=90, code=3,
       col=c("black","red","forestgreen","orange","magenta","cyan2","blue",
             "purple")[as.factor(meansdTrophGr$Trophic_group)])
arrows(meansdTrophGr$d13C, meansdTrophGr$d15N-meansdTrophGr$sd15N, meansdTrophGr$d13C, meansdTrophGr$d15N+meansdTrophGr$sd15N, length=0.05,lwd=2, angle=90, code=3,
       col=c("black","red","forestgreen","orange","magenta","cyan2","blue",
             "purple")[as.factor(meansdTrophGr$Trophic_group)])
points(meansdTrophGr$d15N~meansdTrophGr$d13C,pch=19,cex=3,cex.axis=2,
       font.main=2,cex.main=2,lwd=2,
       col=c("black","red","forestgreen","orange","magenta","cyan2","blue",
             "purple")[as.factor(meansdTrophGr$Trophic_group)])


text(-16,0.75,expression(bold("δ"^"13"*"C (‰)")), cex=3.5,xpd=NA,adj=0.5)
title(ylab=expression(bold("δ"^"15"*"N (‰)")), line=8, cex.lab=3.5)

legend(-27,23.9,c("Unkown",
                  "Predators of active prey",
                  "Predators of large sessile prey",
                  "Predators of encrusting prey",
                  "Suspension feeders",
                  "Sediment feeders",
                  "Omnivores",
                  "Pelagos-based omnivores"),
       pch=19,col=c("black","red","magenta","cyan2","purple","forestgreen","orange",
                    "blue"),cex=2.5,bty="n")

mtext("a",side=3,adj=0.01,line=-3,font=2,cex=2.75)

# Plot displaying siberPSDx modes and credible intervals and PSDxc values for 
# each trophic group.
siberDensityPlot(PSDxTrophicGroup,ylim=c(0,15),xticklabels = rep("",ncol(PSDxTrophicGroup)),
                 xlab="",ylab="",xaxs="i",cex.axis=3.5,font=2, las=2,ax=F,ct="null")

axis(side=2, at=c(0,5,10,15,20), labels=c("0.0","5.0","10.0","15.0","20.0"),las=1,cex.axis=3.5,font=2,lwd=2)
axis(side=1, at=c(1:nrow(modePSDxTrophicGroup)),labels=rep("",nrow(modePSDxTrophicGroup)),lwd=2)
box(lwd=2)

points(c(1:nrow(modePSDxTrophicGroup)),modePSDxTrophicGroup$means,cex=2.75,pch=19)
points(c(1:ncol(groupMLTrophicGroup)),groupMLTrophicGroup[5,],cex=3,pch=24,bg="white")

for (i in 1:nrow(modePSDxTrophicGroup)){
  text(i,as.data.frame(QPSDxTrophicGroup[i])[2,2]+0.5,modePSDxTrophicGroup[i,2],cex=2.5)
}

title(ylab=expression(bold("PSD"["CB"]~"(‰)")), line=8, cex.lab=3.5)

mtext("b",side=3,adj=0.01,line=-3,font=2,cex=2.75)

# Plot displaying siberPSDy modes and credible intervals and PSDyc values for 
# each trophic group.
siberDensityPlot(PSDyTrophicGroup,ylim=c(0,15),xticklabels = rep("",ncol(PSDyTrophicGroup)),
                 xlab="",ylab="",xaxs="i",cex.axis=3.5,font=2, las=2,ax=F,ct="null")
axis(side=2, at=c(0,5,10,15,20), labels=c("0.0","5.0","10.0","15.0","20.0"),las=1,cex.axis=3.5,font=2,lwd=2)
axis(side=1, at=c(1:nrow(modePSDyTrophicGroup)),labels=rep("",nrow(modePSDyTrophicGroup)),lwd=2)
box(lwd=2)

points(c(1:nrow(modePSDyTrophicGroup)),modePSDyTrophicGroup$means,cex=2.75,pch=19)
points(c(1:ncol(groupMLTrophicGroup)),groupMLTrophicGroup[6,],cex=3,pch=24,bg="white")

for (i in 1:nrow(modePSDyTrophicGroup)){
  text(i,as.data.frame(QPSDyTrophicGroup[i])[2,2]+0.5,modePSDyTrophicGroup[i,2],cex=2.5)
}

title(ylab=expression(bold("PSD"["NB"]~"(‰)")), line=8, cex.lab=3.5)

mtext("c",side=3,adj=0.01,line=-3,font=2,cex=2.75)

# x-axis tick labels and names for plots b and c.
labels=c("Unkown",
         "Predators of\nactive prey",
         "Predators of large\nsessile prey",
         "Predators of\nencrusting prey",
         "Sediment feeders",
         "Omnivores",
         "Pelagos-based\nomnivores")
par(mar=c(0,12,0.5,0.3))
plot(c(0,nrow(modePSDxTrophicGroup)),c(0,1),xlab="",ylab="",
     ax=F,col="white",xaxs="i",xlim=c(1-0.5,nrow(modePSDxTrophicGroup)+0.5))
text(x = 1:nrow(modePSDxTrophicGroup),
     y = 1,
     labels =labels,
     xpd = NA,
     srt = 60,
     font=2,
     cex = 2.5,adj=0.9) 
text(x=nrow(modePSDxTrophicGroup)/2+0.5,0,xpd=NA,"Trophic groups",adj=c(0.5,0),font=2,cex = 3.5)

################################################################################
# Isotopic niches for depth intervals
################################################################################

# As depth is a continuous variable and as standard ellipses are generated for 
# factor levels, we need to create depth intervals.
data<-transform(Data, Interval = cut(Depth, breaks = c(0, 99.99, 199.99, 499.99, 999.99, 1499.99, 1999.99, 2999.99, 4499.99,Inf),
                                     labels = c("0-100", "100-200", "200-500","500-1000","1000-1500","1500-2000","2000-3000","3000-4500","4500+")))

# We then select the columns necessary to generate the standard ellipses, i.e.,
# the mean-corrected δ13C and δ15N values, the species, and the chlorophyll
# concentration interval as the community. We also remove rows for which no 
# chlorophyll concentration data is available.
nicheSIBER<-data[,c(23,26,3,27)]

# We give the required names to the columns.
colnames(nicheSIBER) = c("iso1", "iso2", "group","community")

# For each depth interval, one stable isotope value has to be 
# the mean stable isotope value for one taxa, so we compute this mean for both 
# δ13C and δ15N values.
meand13C<-aggregate(iso1~group+community, data=nicheSIBER, FUN="mean")
meand15N<-aggregate(iso2~group+community, data=nicheSIBER, FUN="mean")
nicheSIBER<-left_join(meand13C,meand15N)

# We can now consider all taxa as a single species to ensure that the standard
# ellipses will be generated at the scale of the depth intervals only.
nicheSIBER$group<-rep("AllSpecies",nrow(nicheSIBER))

# Re-ordering of the columns in the right order and creation of the SIBER object.
nicheSIBER<-nicheSIBER[,c(3,4,1,2)]
siber.Test <- createSiberObject(nicheSIBER)

group.ML <- groupMetricsML2(siber.Test)
print(group.ML)

# We can now conduct the Bayesian estimations of PSDx and PSDy.

# Options for running jags
parms <- list()
parms$n.iter <- 10^5   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains
parms$save.output = TRUE
parms$save.dir = tempdir()

# Definition of the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# Fitting of the ellipses.
ellipses.posterior <- siberMVN(siber.Test, parms, priors)

# Bayesian estimation of PSDx with our custom function siberPSDx.
PSDx.B <- siberPSDx(ellipses.posterior)
colnames(PSDx.B)<-colnames(group.ML)

# Estimation of the credible intervals.
PSDx.B.credibles <- lapply(
  as.data.frame(PSDx.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = 0.95)

# Estimation of the modes.
PSDx.B.modes <- lapply(
  as.data.frame(PSDx.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = 0.95, all.modes=T)

# Comparison of PSDx for each pair of ellipse in a matrix.
Q <- matrix(1, ncol = ncol(PSDx.B), nrow = ncol(PSDx.B))
for (i in 1:(ncol(PSDx.B) - 1)) {
  for (j in (i + 1):ncol(PSDx.B)) {
    Q[i, j] <- if(sum( PSDx.B[,i] < PSDx.B[,j] ) / nrow(PSDx.B) < 0.5) {
      sum( PSDx.B[,i] < PSDx.B[,j] ) / nrow(PSDx.B)}
    else{sum( PSDx.B[,i] > PSDx.B[,j] ) / nrow(PSDx.B)}
    Q[j, i] <- if(sum( PSDx.B[,i] < PSDx.B[,j] ) / nrow(PSDx.B) < 0.5) {
      sum( PSDx.B[,i] < PSDx.B[,j] ) / nrow(PSDx.B)}
    else{sum( PSDx.B[,i] > PSDx.B[,j] ) / nrow(PSDx.B)}
  }
}

# Sorting of ellipses in groups with similar PSDx.
mode <- data.frame(matrix(unlist(PSDx.B.modes), nrow=length(PSDx.B.modes), byrow=T))
colnames(mode)<-"mode"
mode$group<-colnames(PSDx.B)
groupsPSDx<- orderPvalue(mode[, 2], mode[, 1], 0.05, Q)
groupsPSDx

# We apply the same process for the Bayesian estimation of PSDy with our custom 
# function siberPSDy.
PSDy.B <- siberPSDy(ellipses.posterior)
colnames(PSDy.B)<-colnames(group.ML)

# Estimation of the credible intervals.
PSDy.B.credibles <- lapply(
  as.data.frame(PSDy.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = 0.95)

# Estimation of the modes.
PSDy.B.modes <- lapply(
  as.data.frame(PSDy.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = 0.95, all.modes=T)

# Comparison of PSDy for each pair of ellipse in a matrix.
Q <- matrix(1, ncol = ncol(PSDy.B), nrow = ncol(PSDy.B))
for (i in 1:(ncol(PSDy.B) - 1)) {
  for (j in (i + 1):ncol(PSDy.B)) {
    Q[i, j] <- if(sum( PSDy.B[,i] < PSDy.B[,j] ) / nrow(PSDy.B) < 0.5) {
      sum( PSDy.B[,i] < PSDy.B[,j] ) / nrow(PSDy.B)}
    else{sum( PSDy.B[,i] > PSDy.B[,j] ) / nrow(PSDy.B)}
    Q[j, i] <- if(sum( PSDy.B[,i] < PSDy.B[,j] ) / nrow(PSDy.B) < 0.5) {
      sum( PSDy.B[,i] < PSDy.B[,j] ) / nrow(PSDy.B)}
    else{sum( PSDy.B[,i] > PSDy.B[,j] ) / nrow(PSDy.B)}
  }
}

# Sorting of ellipses in groups with similar PSDy.
mode <- data.frame(matrix(unlist(PSDy.B.modes), nrow=length(PSDy.B.modes), byrow=T))
colnames(mode)<-"mode"
mode$group<-colnames(PSDy.B)
groupsPSDy<- orderPvalue(mode[, 2], mode[, 1], 0.05, Q)
groupsPSDy

# Summary of the results.

# Preparation of data for the plots and a table summarising the results of PSDx 
# and PSDy estimation. We need to reorder the chlorophyll concentration 
# intervals in the order we want to see them on the table and the plots.
PSDxbathome<-PSDx.B
QPSDxbathome<-PSDx.B.credibles
modePSDxbathome<-as.data.frame(groupsPSDx)
modePSDxbathome$group<-rownames(modePSDxbathome)
modePSDxbathome<-modePSDxbathome[with(modePSDxbathome, order(group)),]
modePSDxbathome<-modePSDxbathome[c(1,2,5,9,3,4,6,7,8),]

PSDybathome<-PSDy.B
QPSDybathome<-PSDy.B.credibles
modePSDybathome<-as.data.frame(groupsPSDy)
modePSDybathome$group<-rownames(modePSDybathome)
modePSDybathome<-modePSDybathome[with(modePSDybathome, order(group)),]
modePSDybathome<-modePSDybathome[c(1,2,5,9,3,4,6,7,8),]

groupMLbathome<-group.ML

# Creation of a table summarising the results of PSDx and PSDy estimation.
tablebathome<-data.frame(matrix(nrow=nrow(modePSDybathome), ncol=3))
colnames(tablebathome)<-c("Bapthome","PSDCB","PSDNB")

for (i in 1:nrow(tablebathome)){
  tablebathome[i,1]<-modePSDxbathome[i,3]
  tablebathome[i,2]<-paste(round(modePSDxbathome[i,1],1),round(as.data.frame(QPSDxbathome[i])[2,1],1),sep=" (")
  tablebathome[i,2]<-paste(tablebathome[i,2],round(as.data.frame(QPSDxbathome[i])[2,2],1),sep="-")
  tablebathome[i,2]<-paste(tablebathome[i,2],modePSDxbathome[i,2],sep=") ")
  tablebathome[i,3]<-paste(round(modePSDybathome[i,1],1),round(as.data.frame(QPSDybathome[i])[2,1],1),sep=" (")
  tablebathome[i,3]<-paste(tablebathome[i,3],round(as.data.frame(QPSDybathome[i])[2,2],1),sep="-")
  tablebathome[i,3]<-paste(tablebathome[i,3],modePSDybathome[i,2],sep=") ")
}

tablebathome

# Creation of plots summarizing the results.
# Window settings.
dev.new(width=17, height=3*5.2, unit="cm")
layout(matrix(c(
  1,3,
  1,3,
  1,3,
  2,4,
  2,4,
  2,4,
  6,5
), nrow=7, ncol=2, byrow = TRUE))
par(mar=c(1.5,12.5,0.5,0.3))

# Plot displaying the relationship between δ13C values and depth.
plot(data$logDepth,data$d13C,pch=21,bg="white",col="black",
     xlab="",ylab="",font.lab=2,cex.lab=2,cex=1,cex.axis=2,ylim=c(-26,-6),xlim=c(log10(5),log10(5000)),
     ax=F,font.main=2,cex.main=2)

axis(side=2, at=c(-26,-21,-16,-11,-6), labels=c("-26.0","-21.0","-16.0","-11.0","-6.0"),lwd=2,cex.axis=3.5,font=2, las=2)
axis(side=1, at=c(log(5,10),log(50,10),log(500,10),log(5000,10)), labels=rep("",4),lwd=2,cex.axis=3.5,font=2)
box(lwd=2)

##### Correlation test between δ13C values and log-transformed depth. The result
##### is then displayed on the plot.
cor.test(data$d13C,data$logDepth)
legend("topright", c(expression(bold("r = -0.702, P < 0.001")),""),cex=2.75, bty="n",xjust=1) 
ablineclip(lm(data$d13C~data$logDepth),lwd=1.5,col="black",x1=min(data$logDepth),x2=max(data$logDepth))

title(ylab=expression(bold("δ"^"13"*"C (‰)")), line=8, cex.lab=3.5)
title(xlab=expression(bold("Depth (m, log"["10"]*"-scale)")), line=8, cex.lab=3.5)

mtext("a",side=3,adj=0.01,line=-3,font=2,cex=2.75)

# Plot displaying the relationship between δ15N values and depth.
plot(data$logDepth,data$d15N,pch=21,bg="white",col="black",
     xlab="",ylab="",font.lab=2,cex.lab=2,cex=1,cex.axis=2,ylim=c(4,24),xlim=c(log10(5),log10(5000)),
     ax=F,font.main=2,cex.main=2)

axis(side=2, at=c(4,9,14,19,24), labels=c("4.0","9.0","14.0","19.0","24.0"),lwd=2,cex.axis=3.5,font=2, las=2)
axis(side=1, at=c(log(5,10),log(50,10),log(500,10),log(5000,10)), labels=rep("",4),lwd=2,cex.axis=3.5,font=2)
axis(side=1, at=c(log(5,10),log(50,10),log(500,10),log(5000,10)), labels=c("5","50","500","5000"),lwd=0,cex.axis=3.5,font=2,line=1)
box(lwd=2)

##### Correlation test between δ15N values and log-transformed depth. The result
##### is then displayed on the plot.
cor.test(data$d15N,data$logDepth)
legend("topright", c(expression(bold("r = 0.075, P < 0.001")),""),cex=2.75, bty="n",xjust=1) 
ablineclip(lm(data$d15N~data$logDepth),lwd=2,col="black",x1=min(data$logDepth),x2=max(data$logDepth))

title(ylab=expression(bold("δ"^"15"*"N (‰)")), line=8, cex.lab=3.5)
text(mean(0.69897+1.5),0.6,expression(bold("Depth (m, log"["10"]*"-scale)")), cex=3.5,xpd=NA,adj=0.5)

mtext("b",side=3,adj=0.01,line=-3,font=2,cex=2.75)

# Plot displaying siberPSDx modes and credible intervals and PSDxc values for 
# each depth interval.
siberDensityPlot(PSDxbathome,ylim=c(0,15),xticklabels = rep("",ncol(PSDxbathome)),
                 xlab="",ylab="",xaxs="i",cex.axis=3.5,font=2, las=2,ax=F,ct="null")

axis(side=2, at=c(0,5,10,15,20), labels=c("0.0","5.0","10.0","15.0","20.0"),las=1,cex.axis=3.5,font=2,lwd=2)
axis(side=1, at=c(1:nrow(modePSDxbathome)),labels=rep("",nrow(modePSDxbathome)),lwd=2)
box(lwd=2)

points(c(1:nrow(modePSDxbathome)),modePSDxbathome$means,cex=2.75,pch=19)
points(c(1:ncol(groupMLbathome)),groupMLbathome[5,],cex=3,pch=24,bg="white")

for (i in 1:nrow(modePSDxbathome)){
  text(i,as.data.frame(QPSDxbathome[i])[2,2]+0.5,modePSDxbathome[i,2],cex=2.5)
}

title(ylab=expression(bold("PSD"["CB"]~"(‰)")), line=8, cex.lab=3.5)

mtext("c",side=3,adj=0.01,line=-3,font=2,cex=2.75)

# Plot displaying siberPSDy modes and credible intervals and PSDyc values for 
# each depth interval.
siberDensityPlot(PSDybathome,ylim=c(0,15),xticklabels = rep("",ncol(PSDybathome)),
                 xlab="",ylab="",xaxs="i",cex.axis=3.5,font=2, las=2,ax=F,ct="null")

axis(side=2, at=c(0,5,10,15,20), labels=c("0.0","5.0","10.0","15.0","20.0"),las=1,cex.axis=3.5,font=2,lwd=2)
axis(side=1, at=c(1:nrow(modePSDybathome)),labels=rep("",nrow(modePSDybathome)),lwd=2)
box(lwd=2)

points(c(1:nrow(modePSDybathome)),modePSDybathome$means,cex=2.75,pch=19)
points(c(1:ncol(groupMLbathome)),groupMLbathome[6,],cex=3,pch=24,bg="white")

for (i in 1:nrow(modePSDybathome)){
  text(i,as.data.frame(QPSDybathome[i])[2,2]+0.5,modePSDybathome[i,2],cex=2.5)
}

title(ylab=expression(bold("PSD"["NB"]~"(‰)")), line=8, cex.lab=3.5)

mtext("d",side=3,adj=0.01,line=-3,font=2,cex=2.75)

# x-axis tick labels and names for plots c and d.
labels=c("[0,100[",
         "[100,200[",
         "[200,500[",
         "[500,1000[",
         "[1000,1500[",
         "[1500,2000[",
         "[2000,3000[",
         "[3000,4500[",
         "[4500,+∞[")
par(mar=c(0,12,0.5,0.3))
plot(c(0,nrow(modePSDxbathome)),c(0,1),xlab="",ylab="",
     ax=F,col="white",xaxs="i",xlim=c(1-0.5,nrow(modePSDxbathome)+0.5))
text(x = 1:nrow(modePSDxbathome),
     y = 1.1,
     labels =labels,
     xpd = NA,
     srt = 60,
     font=2,
     cex = 3.5,adj=1) 
text(x=nrow(modePSDxbathome)/2+0.5,0,xpd=NA,"Bathome (m)",adj=c(0.5,0),font=2,cex = 3.5)

################################################################################
# Isotopic niches for sea ice concentration intervals
################################################################################

# As sea ice concentration is a continuous variable and as standard ellipses are 
# generated for factor levels, we need to create sea ice concentration intervals.
# We also need to specify that some regions where no sea ice is present are actually
# Subantarctic regions.
data<-transform(Data, Interval = cut(seaice_prev_month, breaks = c(-0.99, 9.99, 19.99, 29.99, 39.99, 49.99, 59.99, 69.99, 79.99, 89.99,Inf),
                                     labels = c("0-10", "10-20", "20-30","30-40","40-50","50-60","60-70","70-80","80-90","90-100")))
levels(data$Interval) <- c(levels(data$Interval), "0-10-sub") 
data$Interval[data$region=="Falklands"]  <- "0-10-sub"
data$Interval[data$region=="South Georgia"]  <- "0-10-sub"
data$Interval[data$region=="Kerguelen"]  <- "0-10-sub"
data$Interval[data$region=="Prince Edward Islands"]  <- "0-10-sub"
data$Interval[data$region=="Patagonia"]  <- "0-10-sub"
data$Interval[data$region=="SatlDeep"]  <- "0-10-sub"

# We then select the columns necessary to generate the standard ellipses, i.e.,
# the mean-corrected δ13C and δ15N values, the species, and the sea ice concentration
# interval as the community.
nicheSIBER<-data[,c(23,26,3,27)]

# We give the required names to the columns.
colnames(nicheSIBER) = c("iso1", "iso2", "group","community")

# For each sea ice concentration interval, one stable isotope value has to be 
# the mean stable isotope value for one taxa, so we compute this mean for both 
# δ13C and δ15N values.
meand13C<-aggregate(iso1~group+community, data=nicheSIBER, FUN="mean")
meand15N<-aggregate(iso2~group+community, data=nicheSIBER, FUN="mean")
nicheSIBER<-left_join(meand13C,meand15N)

# We can now consider all taxa as a single species to ensure that the standard
# ellipses will be generated at the scale of the sea ice concentration intervals
# only.
nicheSIBER$group<-rep("AllSpecies",nrow(nicheSIBER))

# Re-ordering of the columns in the right order and creation of the SIBER object.
nicheSIBER<-nicheSIBER[,c(3,4,1,2)]
siber.Test <- createSiberObject(nicheSIBER)

# Computation of PSDxc and PSDyc for each ice concentration interval with our 
# custom function groupMetricsML2.
group.ML <- groupMetricsML2(siber.Test)
print(group.ML)

# We can now conduct the Bayesian estimations of PSDx and PSDy.

# Options for running jags
parms <- list()
parms$n.iter <- 10^5   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains
parms$save.output = TRUE
parms$save.dir = tempdir()

# Definition of the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# Fitting of the ellipses.
ellipses.posterior <- siberMVN(siber.Test, parms, priors)

# Bayesian estimation of PSDx with our custom function siberPSDx.
PSDx.B <- siberPSDx(ellipses.posterior)
colnames(PSDx.B)<-colnames(group.ML)

# Estimation of the credible intervals.
PSDx.B.credibles <- lapply(
  as.data.frame(PSDx.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = 0.95)

# Estimation of the modes.
PSDx.B.modes <- lapply(
  as.data.frame(PSDx.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = 0.95, all.modes=T)

# Comparison of PSDx for each pair of ellipse in a matrix.
Q <- matrix(1, ncol = ncol(PSDx.B), nrow = ncol(PSDx.B))
for (i in 1:(ncol(PSDx.B) - 1)) {
  for (j in (i + 1):ncol(PSDx.B)) {
    Q[i, j] <- if(sum( PSDx.B[,i] < PSDx.B[,j] ) / nrow(PSDx.B) < 0.5) {
      sum( PSDx.B[,i] < PSDx.B[,j] ) / nrow(PSDx.B)}
    else{sum( PSDx.B[,i] > PSDx.B[,j] ) / nrow(PSDx.B)}
    Q[j, i] <- if(sum( PSDx.B[,i] < PSDx.B[,j] ) / nrow(PSDx.B) < 0.5) {
      sum( PSDx.B[,i] < PSDx.B[,j] ) / nrow(PSDx.B)}
    else{sum( PSDx.B[,i] > PSDx.B[,j] ) / nrow(PSDx.B)}
  }
}

# Sorting of ellipses in groups with similar PSDx.
mode <- data.frame(matrix(unlist(PSDx.B.modes), nrow=length(PSDx.B.modes), byrow=T))
colnames(mode)<-"mode"
mode$group<-colnames(PSDx.B)
groupsPSDx<- orderPvalue(mode[, 2], mode[, 1], 0.05, Q)
groupsPSDx

# We apply the same process for the Bayesian estimation of PSDy with our custom 
# function siberPSDy.
PSDy.B <- siberPSDy(ellipses.posterior)
colnames(PSDy.B)<-colnames(group.ML)

# Estimation of the credible intervals.
PSDy.B.credibles <- lapply(
  as.data.frame(PSDy.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = 0.95)

# Estimation of the modes.
PSDy.B.modes <- lapply(
  as.data.frame(PSDy.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = 0.95, all.modes=T)

# Comparison of PSDy for each pair of ellipse in a matrix.
Q <- matrix(1, ncol = ncol(PSDy.B), nrow = ncol(PSDy.B))
for (i in 1:(ncol(PSDy.B) - 1)) {
  for (j in (i + 1):ncol(PSDy.B)) {
    Q[i, j] <- if(sum( PSDy.B[,i] < PSDy.B[,j] ) / nrow(PSDy.B) < 0.5) {
      sum( PSDy.B[,i] < PSDy.B[,j] ) / nrow(PSDy.B)}
    else{sum( PSDy.B[,i] > PSDy.B[,j] ) / nrow(PSDy.B)}
    Q[j, i] <- if(sum( PSDy.B[,i] < PSDy.B[,j] ) / nrow(PSDy.B) < 0.5) {
      sum( PSDy.B[,i] < PSDy.B[,j] ) / nrow(PSDy.B)}
    else{sum( PSDy.B[,i] > PSDy.B[,j] ) / nrow(PSDy.B)}
  }
}

# Sorting of ellipses in groups with similar PSDy.
mode <- data.frame(matrix(unlist(PSDy.B.modes), nrow=length(PSDy.B.modes), byrow=T))
colnames(mode)<-"mode"
mode$group<-colnames(PSDy.B)
groupsPSDy<- orderPvalue(mode[, 2], mode[, 1], 0.05, Q)
groupsPSDy

# Summary of the results.

# Preparation of data for the plots and a table summarising the results of PSDx 
# and PSDy estimation. We need to reorder the sea ice concentration intervals in
# the order we want to see them on the table and the plots.
PSDxIceConc<-PSDx.B[,c(11,1:10)]
QPSDxIceConc<-PSDx.B.credibles[c(11,1:10)]
modePSDxIceConc<-as.data.frame(groupsPSDx)
modePSDxIceConc$group<-rownames(modePSDxIceConc)
modePSDxIceConc<-modePSDxIceConc[with(modePSDxIceConc, order(group)),]

PSDyIceConc<-PSDy.B[,c(11,1:10)]
QPSDyIceConc<-PSDy.B.credibles[c(11,1:10)]
modePSDyIceConc<-as.data.frame(groupsPSDy)
modePSDyIceConc$group<-rownames(modePSDyIceConc)
modePSDyIceConc<-modePSDyIceConc[with(modePSDyIceConc, order(group)),]

groupMLIceConc<-group.ML[,c(11,1:10)]

# Creation of a table summarising the results of PSDx and PSDy estimation.
tableIceConc<-data.frame(matrix(nrow=nrow(modePSDyIceConc), ncol=3))
colnames(tableIceConc)<-c("Sea Ice Concentration","PSDCB","PSDNB")

for (i in 1:nrow(tableIceConc)){
  tableIceConc[i,1]<-modePSDxIceConc[i,3]
  tableIceConc[i,2]<-paste(round(modePSDxIceConc[i,1],1),round(as.data.frame(QPSDxIceConc[i])[2,1],1),sep=" (")
  tableIceConc[i,2]<-paste(tableIceConc[i,2],round(as.data.frame(QPSDxIceConc[i])[2,2],1),sep="-")
  tableIceConc[i,2]<-paste(tableIceConc[i,2],modePSDxIceConc[i,2],sep=") ")
  tableIceConc[i,3]<-paste(round(modePSDyIceConc[i,1],1),round(as.data.frame(QPSDyIceConc[i])[2,1],1),sep=" (")
  tableIceConc[i,3]<-paste(tableIceConc[i,3],round(as.data.frame(QPSDyIceConc[i])[2,2],1),sep="-")
  tableIceConc[i,3]<-paste(tableIceConc[i,3],modePSDyIceConc[i,2],sep=") ")
}

tableIceConc

# Creation of plots summarizing the results.
# Window settings.
dev.new(width=17, height=3*5.2, unit="cm")
layout(matrix(c(
  1,3,
  1,3,
  1,3,
  2,4,
  2,4,
  2,4,
  6,5
), nrow=7, ncol=2, byrow = TRUE))
par(mar=c(1.5,12.5,0.5,0.3))

# Plot displaying the relationship between δ13C values and sea ice concentration.
plot(DataAnt$seaice_prev_month,DataAnt$d13C,pch=21,bg="white",col="black",
     xlab="",ylab="",font.lab=2,cex.lab=2,cex=1,cex.axis=2,ylim=c(-26,-6),
     xlim=c(-1,91),ax=F,font.main=2,cex.main=2)

axis(side=2, at=c(-26,-21,-16,-11,-6), labels=c("-26.0","-21.0","-16.0","-11.0","-6.0"),lwd=2,cex.axis=3.5,font=2, las=2)
axis(side=1, at=c(0,30,60,90), labels=rep("",4),lwd=2,cex.axis=3.5,font=2)
box(lwd=2)

##### Correlation test between δ13C values and logit-transformed sea ice concentration.
##### The result is then displayed on the plot.
cor.test(DataAnt$logIce,DataAnt$d13C)
legend("topright", c(expression(bold("r = 0.082, P < 0.001")),""),cex=2.75, bty="n",xjust=1) 

##### Since we want to show the relationship between δ15N values and sea ice 
##### concentration with a scale of raw sea ice concentration values instead of 
##### logit-transformed ones, we need to predict the relationship between δ13C 
##### values and sea ice concentration.
test.glm <- lm(d13C ~ logit(seaice_prev_month),data=DataAnt)
x_test <- seq(from=0,to=max(DataAnt$seaice_prev_month),length.out=1000)
y_test <- predict(test.glm,newdata = data.frame(seaice_prev_month = x_test),interval="confidence")
lines(x_test, y_test[,1],lwd=1.5,col="black")

title(ylab=expression(bold("δ"^"13"*"C (‰)")), line=8, cex.lab=3.5)
title(xlab=expression(bold("Sea ice concentration (%)")), line=8, cex.lab=3.5)

mtext("a",side=3,adj=0.01,line=-3,font=2,cex=2.75)

# Plot displaying the relationship between δ15N values and sea ice concentration.
plot(DataAnt$seaice_prev_month,DataAnt$d15N,pch=21,bg="white",col="black",
     xlab="",ylab="",font.lab=2,cex.lab=2,cex=1,cex.axis=2,ylim=c(4,24),xlim=c(-1,91),
     ax=F,font.main=2,cex.main=2)

axis(side=2, at=c(4,9,14,19,24), labels=c("4.0","9.0","14.0","19.0","24.0"),lwd=2,cex.axis=3.5,font=2, las=2)
axis(side=1, at=c(0,30,60,90), labels=rep("",4),lwd=2,cex.axis=3.5,font=2)
axis(side=1, at=c(0,30,60,90), labels=c("0","30","60","90"),lwd=0,cex.axis=3.5,font=2,line=1)
box(lwd=2)

##### Correlation test between δ15N values and logit-transformed sea ice concentration.
##### The result is then displayed on the plot.
cor.test(DataAnt$d15N,DataAnt$logIce)
legend("topright", c(expression(bold("r = 0.133, P < 0.001")),""),cex=2.75, bty="n",xjust=1) 

##### Since we want to show the relationship between δ15N values and sea ice 
##### concentration with a scale of raw sea ice concentration values instead of 
##### logit-transformed ones, we need to predict the relationship between δ13C 
##### values and sea ice concentration.
test.glm <- lm(d15N ~ logit(seaice_prev_month),data=DataAnt)
x_test <- seq(from=0,to=max(DataAnt$seaice_prev_month),length.out=1000)
y_test <- predict(test.glm,newdata = data.frame(seaice_prev_month = x_test),interval="confidence")
lines(x_test, y_test[,1],lwd=1.5,col="black")

title(ylab=expression(bold("δ"^"15"*"N (‰)")), line=8, cex.lab=3.5)
text(mean(45),0.6,expression(bold("Sea ice concentration (%)")), cex=3.5,xpd=NA,adj=0.5)

mtext("b",side=3,adj=0.01,line=-3,font=2,cex=2.75)

# Plot displaying siberPSDx modes and credible intervals and PSDxc values for 
# each sea ice concentration interval.
siberDensityPlot(PSDxIceConc,ylim=c(0,15),xticklabels = rep("",ncol(PSDxIceConc)),
                 xlab="",ylab="",xaxs="i",cex.axis=3.5,font=2, las=2,ax=F,ct="null")

axis(side=2, at=c(0,5,10,15,20), labels=c("0.0","5.0","10.0","15.0","20.0"),las=1,cex.axis=3.5,font=2,lwd=2)
axis(side=1, at=c(1:nrow(modePSDxIceConc)),labels=rep("",nrow(modePSDxIceConc)),lwd=2)
box(lwd=2)

points(c(1:nrow(modePSDxIceConc)),modePSDxIceConc$means,cex=2.75,pch=19)
points(c(1:ncol(groupMLIceConc)),groupMLIceConc[5,],cex=3,pch=24,bg="white")

for (i in 1:nrow(modePSDxIceConc)){
  text(i,as.data.frame(QPSDxIceConc[i])[2,2]+0.5,modePSDxIceConc[i,2],cex=2.5)
}

title(ylab=expression(bold("PSD"["CB"]~"(‰)")), line=8, cex.lab=3.5)

mtext("c",side=3,adj=0.01,line=-3,font=2,cex=2.75)

# Plot displaying siberPSDy modes and credible intervals and PSDyc values for 
# each sea ice concentration interval.
siberDensityPlot(PSDyIceConc,ylim=c(0,15),xticklabels = rep("",ncol(PSDyIceConc)),
                 xlab="",ylab="",xaxs="i",cex.axis=3.5,font=2, las=2,ax=F,ct="null")

axis(side=2, at=c(0,5,10,15,20), labels=c("0.0","5.0","10.0","15.0","20.0"),las=1,cex.axis=3.5,font=2,lwd=2)
axis(side=1, at=c(1:nrow(modePSDyIceConc)),labels=rep("",nrow(modePSDyIceConc)),lwd=2)
box(lwd=2)

points(c(1:nrow(modePSDyIceConc)),modePSDyIceConc$means,cex=2.75,pch=19)
points(c(1:ncol(groupMLIceConc)),groupMLIceConc[6,],cex=3,pch=24,bg="white")

for (i in 1:nrow(modePSDyIceConc)){
  text(i,as.data.frame(QPSDyIceConc[i])[2,2]+0.5,modePSDyIceConc[i,2],cex=2.5)
}

title(ylab=expression(bold("PSD"["NB"]~"(‰)")), line=8, cex.lab=3.5)

mtext("d",side=3,adj=0.01,line=-3,font=2,cex=2.75)

# x-axis tick labels and names for plots c and d.
labels=c("Subantarctic",
         "[0,10[",
         "[10,20[",
         "[20,30[",
         "[30,40[",
         "[40,50[",
         "[50,60[",
         "[60,70[",
         "[70,80[",
         "[80,90[",
         "[90,100]")
par(mar=c(0,12,0.5,0.3))
plot(c(0,nrow(modePSDxIceConc)),c(0,1),xlab="",ylab="",
     ax=F,col="white",xaxs="i",xlim=c(1-0.5,nrow(modePSDxIceConc)+0.5))
text(x = 1:nrow(modePSDxIceConc),
     y = 1.1,
     labels =labels,
     xpd = NA,
     srt = 60,
     font=2,
     cex = 3.5,adj=1) 
text(x=nrow(modePSDxIceConc)/2+0.5,0,xpd=NA,"Sea ice concentration (%)",adj=c(0.5,0),font=2,cex = 3.5)

################################################################################
# Isotopic niches for sea ice duration intervals
################################################################################

# As sea ice duration is a continuous variable and as standard ellipses are 
# generated for factor levels, we need to create sea ice duration intervals. We
# also need to specify that some regions where no sea ice is present are actually
# Subantarctic regions.
data<-transform(Data, Interval = cut(seaice_last_730, breaks = c(-0.99, 60, 121, 182, 243, 304,365,426,487,548,609,670, Inf),
                                     labels = c("0-60", "61-121", "122-182","183-243","244-304","305-365",
                                                "366-426","427-487","488-548","549-609","610-670","671-730")))
levels(data$Interval) <- c(levels(data$Interval), "subantarctic") 
data$Interval[data$region=="Falklands"]  <- "subantarctic"
data$Interval[data$region=="South Georgia"]  <- "subantarctic"
data$Interval[data$region=="Kerguelen"]  <- "subantarctic"
data$Interval[data$region=="Prince Edward Islands"]  <- "subantarctic"
data$Interval[data$region=="Patagonia"]  <- "subantarctic"
data$Interval[data$region=="SatlDeep"]  <- "subantarctic"


# We then select the columns necessary to generate the standard ellipses, i.e.,
# the mean-corrected δ13C and δ15N values, the species, and the sea ice duration
# interval as the community.
nicheSIBER<-data[,c(23,26,3,27)]

# We give the required names to the columns.
colnames(nicheSIBER) = c("iso1", "iso2", "group","community")

# For each sea ice duration interval, one stable isotope value has to be the 
# mean stable isotope value for one taxa, so we compute this mean for both δ13C 
# and δ15N values.
meand13C<-aggregate(iso1~group+community, data=nicheSIBER, FUN="mean")
meand15N<-aggregate(iso2~group+community, data=nicheSIBER, FUN="mean")
nicheSIBER<-left_join(meand13C,meand15N)

# We can now considered all taxa as a single species to ensure that the standard
# ellipses will be generated at the scale of the sea ice concentration intervals
# only.
nicheSIBER$group<-rep("AllSpecies",nrow(nicheSIBER))

# Re-ordering of the columns in the right order and creation of the SIBER object.
nicheSIBER<-nicheSIBER[,c(3,4,1,2)]
siber.Test <- createSiberObject(nicheSIBER)

# Computation of PSDxc and PSDyc for each ice duration interval with our custom 
# function groupMetricsML2.
group.ML <- groupMetricsML2(siber.Test)
print(group.ML)

# We can now conduct the Bayesian estimations of PSDx and PSDy.

# Options for running jags
parms <- list()
parms$n.iter <- 10^5   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains
parms$save.output = TRUE
parms$save.dir = tempdir()

# Definition of the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# Fitting of the ellipses.
ellipses.posterior <- siberMVN(siber.Test, parms, priors)

# Bayesian estimation of PSDx with our custom function siberPSDx.
PSDx.B <- siberPSDx(ellipses.posterior)
colnames(PSDx.B)<-colnames(group.ML)

# Estimation of the credible intervals.
PSDx.B.credibles <- lapply(
  as.data.frame(PSDx.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = 0.95)

# Estimation of the modes.
PSDx.B.modes <- lapply(
  as.data.frame(PSDx.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = 0.95, all.modes=T)

# Comparison of PSDx for each pair of ellipse in a matrix.
Q <- matrix(1, ncol = ncol(PSDx.B), nrow = ncol(PSDx.B))
for (i in 1:(ncol(PSDx.B) - 1)) {
  for (j in (i + 1):ncol(PSDx.B)) {
    Q[i, j] <- if(sum( PSDx.B[,i] < PSDx.B[,j] ) / nrow(PSDx.B) < 0.5) {
      sum( PSDx.B[,i] < PSDx.B[,j] ) / nrow(PSDx.B)}
    else{sum( PSDx.B[,i] > PSDx.B[,j] ) / nrow(PSDx.B)}
    Q[j, i] <- if(sum( PSDx.B[,i] < PSDx.B[,j] ) / nrow(PSDx.B) < 0.5) {
      sum( PSDx.B[,i] < PSDx.B[,j] ) / nrow(PSDx.B)}
    else{sum( PSDx.B[,i] > PSDx.B[,j] ) / nrow(PSDx.B)}
  }
}

# Sorting of ellipses in groups with similar PSDx.
mode <- data.frame(matrix(unlist(PSDx.B.modes), nrow=length(PSDx.B.modes), byrow=T))
colnames(mode)<-"mode"
mode$group<-colnames(PSDx.B)
groupsPSDx<- orderPvalue(mode[, 2], mode[, 1], 0.05, Q)
groupsPSDx

# We apply the same process for the Bayesian estimation of PSDy with our custom 
# function siberPSDy.
PSDy.B <- siberPSDy(ellipses.posterior)
colnames(PSDy.B)<-colnames(group.ML)

# Estimation of the credible intervals.
PSDy.B.credibles <- lapply(
  as.data.frame(PSDy.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = 0.95)

# Estimation of the modes.
PSDy.B.modes <- lapply(
  as.data.frame(PSDy.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = 0.95, all.modes=T)

# Comparison of PSDy for each pair of ellipse in a matrix.
Q <- matrix(1, ncol = ncol(PSDy.B), nrow = ncol(PSDy.B))
for (i in 1:(ncol(PSDy.B) - 1)) {
  for (j in (i + 1):ncol(PSDy.B)) {
    Q[i, j] <- if(sum( PSDy.B[,i] < PSDy.B[,j] ) / nrow(PSDy.B) < 0.5) {
      sum( PSDy.B[,i] < PSDy.B[,j] ) / nrow(PSDy.B)}
    else{sum( PSDy.B[,i] > PSDy.B[,j] ) / nrow(PSDy.B)}
    Q[j, i] <- if(sum( PSDy.B[,i] < PSDy.B[,j] ) / nrow(PSDy.B) < 0.5) {
      sum( PSDy.B[,i] < PSDy.B[,j] ) / nrow(PSDy.B)}
    else{sum( PSDy.B[,i] > PSDy.B[,j] ) / nrow(PSDy.B)}
  }
}

# Sorting of ellipses in groups with similar PSDy.
mode <- data.frame(matrix(unlist(PSDy.B.modes), nrow=length(PSDy.B.modes), byrow=T))
colnames(mode)<-"mode"
mode$group<-colnames(PSDy.B)
groupsPSDy<- orderPvalue(mode[, 2], mode[, 1], 0.05, Q)
groupsPSDy

# Summary of the results.

# Preparation of data for the plots and a table summarising the results of PSDx 
# and PSDy estimation. We need to reorder the sea ice duration intervals in the 
# order we want to see them on the table and the plots.
PSDxIceDur<-PSDx.B[,c(13, 1:12)]
QPSDxIceDur<-PSDx.B.credibles[c(13, 1:12)]
modePSDxIceDur<-as.data.frame(groupsPSDx)
modePSDxIceDur$group<-rownames(modePSDxIceDur)
modePSDxIceDur<-modePSDxIceDur[with(modePSDxIceDur, order(group)),]
modePSDxIceDur<-modePSDxIceDur[c(13,1,10,2:9,11,12),]

PSDyIceDur<-PSDy.B[,c(13, 1:12)]
QPSDyIceDur<-PSDy.B.credibles[c(13, 1:12)]
modePSDyIceDur<-as.data.frame(groupsPSDy)
modePSDyIceDur$group<-rownames(modePSDyIceDur)
modePSDyIceDur<-modePSDyIceDur[with(modePSDyIceDur, order(group)),]
modePSDyIceDur<-modePSDyIceDur[c(13,1,10,2:9,11,12),]

groupMLIceDur<-group.ML[,c(13,1:12)]

####### Table and figure ######
# Creation of a table summarising the results of PSDx and PSDy estimation.
tableIceDur<-data.frame(matrix(nrow=nrow(modePSDyIceDur), ncol=3))
colnames(tableIceDur)<-c("Sea Ice duration","PSDCB","PSDNB")

for (i in 1:nrow(tableIceDur)){
  tableIceDur[i,1]<-modePSDxIceDur[i,3]
  tableIceDur[i,2]<-paste(round(modePSDxIceDur[i,1],1),round(as.data.frame(QPSDxIceDur[i])[2,1],1),sep=" (")
  tableIceDur[i,2]<-paste(tableIceDur[i,2],round(as.data.frame(QPSDxIceDur[i])[2,2],1),sep="-")
  tableIceDur[i,2]<-paste(tableIceDur[i,2],modePSDxIceDur[i,2],sep=") ")
  tableIceDur[i,3]<-paste(round(modePSDyIceDur[i,1],1),round(as.data.frame(QPSDyIceDur[i])[2,1],1),sep=" (")
  tableIceDur[i,3]<-paste(tableIceDur[i,3],round(as.data.frame(QPSDyIceDur[i])[2,2],1),sep="-")
  tableIceDur[i,3]<-paste(tableIceDur[i,3],modePSDyIceDur[i,2],sep=") ")
}

tableIceDur

# Creation of plots summarizing the results.
# Window settings.
dev.new(width=17, height=3*5.2, unit="cm")
layout(matrix(c(
  1,3,
  1,3,
  1,3,
  2,4,
  2,4,
  2,4,
  6,5
), nrow=7, ncol=2, byrow = TRUE))
par(mar=c(1.5,12.5,0.5,0.3))

# Plot displaying the relationship between δ13C values and sea ice duration.
plot(DataAnt$seaice_last_730,DataAnt$d13C,pch=21,bg="white",col="black",
     xlab="",ylab="",font.lab=2,cex.lab=2,cex=1,cex.axis=2,ylim=c(-26,-6),
     xlim=c(0,730),ax=F,font.main=2,cex.main=2)

axis(side=2, at=c(-26,-21,-16,-11,-6), labels=c("-26.0","-21.0","-16.0","-11.0","-6.0"),lwd=2,cex.axis=3.5,font=2, las=2)
axis(side=1, at=c(0,146,292,438,584,730), labels=rep("",6),lwd=2,cex.axis=3.5,font=2)
box(lwd=2)

##### Correlation test between δ13C values and sea ice duration. The result is 
##### then displayed on the plot.
cor.test(DataAnt$seaice_last_730,DataAnt$d13C)
legend("topright", c(expression(bold("r = -0.248, P < 0.001")),""),cex=2.75, bty="n",xjust=1) 
ablineclip(lm(data$d13C~data$seaice_last_730),lwd=1.5,col="black",x1=min(DataAnt$seaice_last_730),x2=730)

title(ylab=expression(bold("δ"^"13"*"C (‰)")), line=8, cex.lab=3.5)
title(xlab=expression(bold("Sea ice duration (days.year"^"-1"*")")), line=8, cex.lab=3.5)

mtext("a",side=3,adj=0.01,line=-3,font=2,cex=2.75)

# Plot displaying the relationship between δ15N values and sea ice duration.
plot(DataAnt$seaice_last_730,DataAnt$d15N,pch=21,bg="white",col="black",
     xlab="",ylab="",font.lab=2,cex.lab=2,cex=1,cex.axis=2,ylim=c(4,24),xlim=c(0,730),
     ax=F,font.main=2,cex.main=2)

axis(side=2, at=c(4,9,14,19,24), labels=c("4.0","9.0","14.0","19.0","24.0"),lwd=2,cex.axis=3.5,font=2, las=2)
axis(side=1, at=c(0,146,292,438,584,730), labels=rep("",6),lwd=2,cex.axis=3.5,font=2)
axis(side=1, at=c(0,146,292,438,584,730), labels=c("0","146","292","438","584","730"),lwd=0,cex.axis=3.5,font=2,line=1)
box(lwd=2)

##### Correlation test between δ15N values and sea ice duration. The result is 
##### then displayed on the plot.
cor.test(DataAnt$d15N,DataAnt$seaice_last_730)
legend("topright", c(expression(bold("r = 0.344, P < 0.001")),""),cex=2.75, bty="n",xjust=1) 
ablineclip(lm(data$d15N~data$seaice_last_730),lwd=2,col="black",x1=min(DataAnt$seaice_last_730),x2=730)

title(ylab=expression(bold("δ"^"15"*"N (‰)")), line=8, cex.lab=3.5)
text(mean(365.5),0.6,expression(bold("Sea ice duration (days.year"^"-1"*")")), cex=3.5,xpd=NA,adj=0.5)

mtext("b",side=3,adj=0.01,line=-3,font=2,cex=2.75)

# Plot displaying siberPSDx modes and credible intervals and PSDxc values for 
# each sea ice concentration interval.
siberDensityPlot(PSDxIceDur,ylim=c(0,15),xticklabels = rep("",ncol(PSDxIceDur)),
                 xlab="",ylab="",xaxs="i",cex.axis=3.5,font=2, las=2,ax=F,ct="null")

axis(side=2, at=c(0,5,10,15,20), labels=c("0.0","5.0","10.0","15.0","20.0"),las=1,cex.axis=3.5,font=2,lwd=2)
axis(side=1, at=c(1:nrow(modePSDxIceDur)),labels=rep("",nrow(modePSDxIceDur)),lwd=2)
box(lwd=2)

points(c(1:nrow(modePSDxIceDur)),modePSDxIceDur$means,cex=2.75,pch=19)
points(c(1:ncol(groupMLIceDur)),groupMLIceDur[5,],cex=3,pch=24,bg="white")

for (i in 1:nrow(modePSDxIceDur)){
  text(i,as.data.frame(QPSDxIceDur[i])[2,2]+0.5,modePSDxIceDur[i,2],cex=2.5)
}

title(ylab=expression(bold("PSD"["CB"]~"(‰)")), line=8, cex.lab=3.5)

mtext("c",side=3,adj=0.01,line=-3,font=2,cex=2.75)

# Plot displaying siberPSDy modes and credible intervals and PSDyc values for 
# each sea ice duration interval.
siberDensityPlot(PSDyIceDur,ylim=c(0,15),xticklabels = rep("",ncol(PSDyIceDur)),
                 xlab="",ylab="",xaxs="i",cex.axis=3.5,font=2, las=2,ax=F,ct="null")

axis(side=2, at=c(0,5,10,15,20), labels=c("0.0","5.0","10.0","15.0","20.0"),las=1,cex.axis=3.5,font=2,lwd=2)
axis(side=1, at=c(1:nrow(modePSDyIceDur)),labels=rep("",nrow(modePSDyIceDur)),lwd=2)
box(lwd=2)

points(c(1:nrow(modePSDyIceDur)),modePSDyIceDur$means,cex=2.75,pch=19)
points(c(1:ncol(groupMLIceDur)),groupMLIceDur[6,],cex=3,pch=24,bg="white")

for (i in 1:nrow(modePSDyIceDur)){
  text(i,as.data.frame(QPSDyIceDur[i])[2,2]+0.5,modePSDyIceDur[i,2],cex=2.5)
}

title(ylab=expression(bold("PSD"["NB"]~"(‰)")), line=8, cex.lab=3.5)

mtext("d",side=3,adj=0.01,line=-3,font=2,cex=2.75)

# x-axis tick labels and names for plots c and d.
labels=c("Subantarctic",
         "[0,60]",
         "[61,121]",
         "[122,182]",
         "[183,243]",
         "[244,304]",
         "[305,365]",
         "[366,426]",
         "[427,487]",
         "[488,548]",
         "[549,609]",
         "[610,670]",
         "[671,730]")
par(mar=c(0,12,0.5,0.3))
plot(c(0,nrow(modePSDxIceDur)),c(0,1),xlab="",ylab="",
     ax=F,col="white",xaxs="i",xlim=c(1-0.5,nrow(modePSDxIceDur)+0.5))
text(x = 1:nrow(modePSDxIceDur),
     y = 1.1,
     labels =labels,
     xpd = NA,
     srt = 60,
     font=2,
     cex = 3.5,adj=1) 
text(x=nrow(modePSDxIceDur)/2+0.5,0,xpd=NA,expression(bold("Sea ice duration (days.year"^"-1"*")")),adj=c(0.5,0),font=2,cex = 3.5)

################################################################################
# Isotopic niches for chlorophyll concentration intervals
################################################################################

# As chlorophyll concentration is a continuous variable and as standard ellipses
# are generated for factor levels, we need to create chlorophyll concentration 
# intervals.
data<-transform(Data, Interval = cut(chl_prev_month, breaks = c(0.01, 0.1, 1, Inf),
                                     labels = c("0.01-0.1", "0.1-1", "1+")))

# We then select the columns necessary to generate the standard ellipses, i.e.,
# the mean-corrected δ13C and δ15N values, the species, and the chlorophyll
# concentration interval as the community. We also remove rows for which no 
# chlorophyll concentration data is available.
nicheSIBER<-data[,c(23,26,3,27)]
nicheSIBER<-droplevels(nicheSIBER[!is.na(nicheSIBER$Interval),])

# We give the required names to the columns.
colnames(nicheSIBER) = c("iso1", "iso2", "group","community")

# For each chlorophyll concentration interval, one stable isotope value has to be 
# the mean stable isotope value for one taxa, so we compute this mean for both 
# δ13C and δ15N values.
meand13C<-aggregate(iso1~group+community, data=nicheSIBER, FUN="mean")
meand15N<-aggregate(iso2~group+community, data=nicheSIBER, FUN="mean")
nicheSIBER<-left_join(meand13C,meand15N)

# We can now consider all taxa as a single species to ensure that the standard
# ellipses will be generated at the scale of the chlorophyll concentration 
# intervals only.
nicheSIBER$group<-rep("AllSpecies",nrow(nicheSIBER))

# Re-ordering of the columns in the right order and creation of the SIBER object.
nicheSIBER<-nicheSIBER[,c(3,4,1,2)]
siber.Test <- createSiberObject(nicheSIBER)

# Computation of PSDxc and PSDyc for chlorophyll concentration interval with our 
# custom function groupMetricsML2.
group.ML <- groupMetricsML2(siber.Test)
print(group.ML)

# We can now conduct the Bayesian estimations of PSDx and PSDy.

# Options for running jags
parms <- list()
parms$n.iter <- 10^5   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains
parms$save.output = TRUE
parms$save.dir = tempdir()

# Definition of the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# Fitting of the ellipses.
ellipses.posterior <- siberMVN(siber.Test, parms, priors)

# Bayesian estimation of PSDx with our custom function siberPSDx.
PSDx.B <- siberPSDx(ellipses.posterior)
colnames(PSDx.B)<-colnames(group.ML)

# Estimation of the credible intervals.
PSDx.B.credibles <- lapply(
  as.data.frame(PSDx.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = 0.95)

# Estimation of the modes.
PSDx.B.modes <- lapply(
  as.data.frame(PSDx.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = 0.95, all.modes=T)

# Comparison of PSDx for each pair of ellipse in a matrix.
Q <- matrix(1, ncol = ncol(PSDx.B), nrow = ncol(PSDx.B))
for (i in 1:(ncol(PSDx.B) - 1)) {
  for (j in (i + 1):ncol(PSDx.B)) {
    Q[i, j] <- if(sum( PSDx.B[,i] < PSDx.B[,j] ) / nrow(PSDx.B) < 0.5) {
      sum( PSDx.B[,i] < PSDx.B[,j] ) / nrow(PSDx.B)}
    else{sum( PSDx.B[,i] > PSDx.B[,j] ) / nrow(PSDx.B)}
    Q[j, i] <- if(sum( PSDx.B[,i] < PSDx.B[,j] ) / nrow(PSDx.B) < 0.5) {
      sum( PSDx.B[,i] < PSDx.B[,j] ) / nrow(PSDx.B)}
    else{sum( PSDx.B[,i] > PSDx.B[,j] ) / nrow(PSDx.B)}
  }
}

# Sorting of ellipses in groups with similar PSDx.
mode <- data.frame(matrix(unlist(PSDx.B.modes), nrow=length(PSDx.B.modes), byrow=T))
colnames(mode)<-"mode"
mode$group<-colnames(PSDx.B)
groupsPSDx<- orderPvalue(mode[, 2], mode[, 1], 0.05, Q)
groupsPSDx

# We apply the same process for the Bayesian estimation of PSDy with our custom 
# function siberPSDy.
PSDy.B <- siberPSDy(ellipses.posterior)
colnames(PSDy.B)<-colnames(group.ML)

# Estimation of the credible intervals.
PSDy.B.credibles <- lapply(
  as.data.frame(PSDy.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = 0.95)

# Estimation of the modes.
PSDy.B.modes <- lapply(
  as.data.frame(PSDy.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = 0.95, all.modes=T)

# Comparison of PSDy for each pair of ellipse in a matrix.
Q <- matrix(1, ncol = ncol(PSDy.B), nrow = ncol(PSDy.B))
for (i in 1:(ncol(PSDy.B) - 1)) {
  for (j in (i + 1):ncol(PSDy.B)) {
    Q[i, j] <- if(sum( PSDy.B[,i] < PSDy.B[,j] ) / nrow(PSDy.B) < 0.5) {
      sum( PSDy.B[,i] < PSDy.B[,j] ) / nrow(PSDy.B)}
    else{sum( PSDy.B[,i] > PSDy.B[,j] ) / nrow(PSDy.B)}
    Q[j, i] <- if(sum( PSDy.B[,i] < PSDy.B[,j] ) / nrow(PSDy.B) < 0.5) {
      sum( PSDy.B[,i] < PSDy.B[,j] ) / nrow(PSDy.B)}
    else{sum( PSDy.B[,i] > PSDy.B[,j] ) / nrow(PSDy.B)}
  }
}

# Sorting of ellipses in groups with similar PSDy.
mode <- data.frame(matrix(unlist(PSDy.B.modes), nrow=length(PSDy.B.modes), byrow=T))
colnames(mode)<-"mode"
mode$group<-colnames(PSDy.B)
groupsPSDy<- orderPvalue(mode[, 2], mode[, 1], 0.05, Q)
groupsPSDy

# Summary of the results.

# Preparation of data for the plots and a table summarising the results of PSDx 
# and PSDy estimation. We need to reorder the chlorophyll concentration 
# intervals in the order we want to see them on the table and the plots.
PSDxChlConc<-PSDx.B
QPSDxChlConc<-PSDx.B.credibles
modePSDxChlConc<-as.data.frame(groupsPSDx)
modePSDxChlConc$group<-rownames(modePSDxChlConc)
modePSDxChlConc<-modePSDxChlConc[with(modePSDxChlConc, order(group)),]

PSDyChlConc<-PSDy.B
QPSDyChlConc<-PSDy.B.credibles
modePSDyChlConc<-as.data.frame(groupsPSDy)
modePSDyChlConc$group<-rownames(modePSDyChlConc)
modePSDyChlConc<-modePSDyChlConc[with(modePSDyChlConc, order(group)),]

groupMLChlConc<-group.ML

# Creation of a table summarising the results of PSDx and PSDy estimation.
tableChlConc<-data.frame(matrix(nrow=nrow(modePSDyChlConc), ncol=3))
colnames(tableChlConc)<-c("Sea Ice Concentration","PSDCB","PSDNB")

for (i in 1:nrow(tableChlConc)){
  tableChlConc[i,1]<-modePSDxChlConc[i,3]
  tableChlConc[i,2]<-paste(round(modePSDxChlConc[i,1],1),round(as.data.frame(QPSDxChlConc[i])[2,1],1),sep=" (")
  tableChlConc[i,2]<-paste(tableChlConc[i,2],round(as.data.frame(QPSDxChlConc[i])[2,2],1),sep="-")
  tableChlConc[i,2]<-paste(tableChlConc[i,2],modePSDxChlConc[i,2],sep=") ")
  tableChlConc[i,3]<-paste(round(modePSDyChlConc[i,1],1),round(as.data.frame(QPSDyChlConc[i])[2,1],1),sep=" (")
  tableChlConc[i,3]<-paste(tableChlConc[i,3],round(as.data.frame(QPSDyChlConc[i])[2,2],1),sep="-")
  tableChlConc[i,3]<-paste(tableChlConc[i,3],modePSDyChlConc[i,2],sep=") ")
}

tableChlConc

# Creation of plots summarizing the results.
# Window settings.
dev.new(width=17, height=3*(5.8-0.4461538), unit="cm")
layout(matrix(c(
  1,3,
  1,3,
  1,3,
  1,3,
  1,3,
  2,4,
  2,4,
  2,4,
  2,4,
  2,4,
  6,5,
  6,5
), nrow=12, ncol=2, byrow = TRUE))
par(mar=c(1.5,12,0.5,0.3))

# Plot displaying the relationship between δ13C values and chlorophyll concentration.
plot(Data$logChl,Data$d13C,pch=21,bg="white",col="black",
     xlab="",ylab="",font.lab=2,cex.lab=2,cex=1,cex.axis=2,ylim=c(-26,-6),
     xlim=c(log(0.02),log(10)),ax=F,font.main=2,cex.main=2)

axis(side=2, at=c(-26,-21,-16,-11,-6), labels=c("-26.0","-21.0","-16.0","-11.0","-6.0"),lwd=2,cex.axis=3.5,font=2, las=2)
axis(side=1, at=c(log(0.02),log(0.10),log(1),log(10)), labels=rep("",4),lwd=2,cex.axis=3.5,font=2)
box(lwd=2)

##### Correlation test between δ13C values and log-transformed chlorophyll
##### concentration. The result is then displayed on the plot.
cor.test(data$logChl,data$d13C)
legend("topright", c(expression(bold("r = 0.333, P < 0.001")),""),cex=2.75, bty="n",xjust=1) 
ablineclip(lm(data$d13C~data$logChl),lwd=1.5,col="black",x1=min(data$logChl,na.rm=T),x2=max(data$logChl,na.rm=T))

title(ylab=expression(bold("δ"^"13"*"C (‰)")), line=8, cex.lab=3.5)

mtext("a",side=3,adj=0.01,line=-3,font=2,cex=2.75)

# Plot displaying the relationship between δ15N values and chlorophyll concentration.
plot(data$logChl,data$d15N,pch=21,bg="white",col="black",
     xlab="",ylab="",font.lab=2,cex.lab=2,cex=1,cex.axis=2,ylim=c(4,24),xlim=c(log(0.02),log(10)),
     ax=F,font.main=2,cex.main=2)

axis(side=2, at=c(4,9,14,19,24), labels=c("4.0","9.0","14.0","19.0","24.0"),lwd=2,cex.axis=3.5,font=2, las=2)
axis(side=1, at=c(log(0.02),log(0.1),log(1),log(10)), labels=rep("",4),lwd=2,cex.axis=3.5,font=2)
axis(side=1, at=c(log(0.02),log(0.1),log(1),log(10)), labels=c("0.02","0.10","1.00","10.00"),lwd=0,cex.axis=3.5,font=2,line=1)
box(lwd=2)

##### Correlation test between δ15N values and log-transformed chlorophyll
##### concentration. The result is then displayed on the plot.
cor.test(data$d15N,data$logChl)
legend("topright", c(expression(bold("r = 0.069, P = 0.021")),""),cex=2.75, bty="n",xjust=1) 
ablineclip(lm(data$d15N~data$logChl),lwd=2,col="black",x1=min(data$logChl,na.rm=T),x2=max(data$logChl,na.rm=T))

title(ylab=expression(bold("δ"^"15"*"N (‰)")), line=8, cex.lab=3.5)
text(mean(-0.804719),0,expression(bold(atop("Surface chlorophyll"~bolditalic("a")~"concentration","(mg.m"^"-3"*", log-scale)"))), cex=3.5,xpd=NA,adj=0.5)

mtext("b",side=3,adj=0.01,line=-3,font=2,cex=2.75)

# Plot displaying siberPSDx modes and credible intervals and PSDxc values for 
# each chlorophyll concentration interval.
siberDensityPlot(PSDxChlConc,ylim=c(0,15),xticklabels = rep("",ncol(PSDxChlConc)),
                 xlab="",ylab="",xaxs="i",cex.axis=3.5,font=2, las=2,ax=F,ct="null")

axis(side=2, at=c(0,5,10,15,20), labels=c("0.0","5.0","10.0","15.0","20.0"),las=1,cex.axis=3.5,font=2,lwd=2)
axis(side=1, at=c(1:nrow(modePSDxChlConc)),labels=rep("",nrow(modePSDxChlConc)),lwd=2)
box(lwd=2)

points(c(1:nrow(modePSDxChlConc)),modePSDxChlConc$means,cex=2.75,pch=19)
points(c(1:ncol(groupMLChlConc)),groupMLChlConc[5,],cex=3,pch=24,bg="white")

for (i in 1:nrow(modePSDxChlConc)){
  text(i,as.data.frame(QPSDxChlConc[i])[2,2]+0.5,modePSDxChlConc[1,2],cex=2.5)
}

title(ylab=expression(bold("PSD"["CB"]~"(‰)")), line=8, cex.lab=3.5)

mtext("c",side=3,adj=0.01,line=-3,font=2,cex=2.75)

# Plot displaying siberPSDy modes and credible intervals and PSDyc values for 
# each chlorophyll concentration interval.
siberDensityPlot(PSDyChlConc,ylim=c(0,15),xticklabels = rep("",ncol(PSDyChlConc)),
                 xlab="",ylab="",xaxs="i",cex.axis=3.5,font=2, las=2,ax=F,ct="null")

axis(side=2, at=c(0,5,10,15,20), labels=c("0.0","5.0","10.0","15.0","20.0"),las=1,cex.axis=3.5,font=2,lwd=2)
axis(side=1, at=c(1:nrow(modePSDyChlConc)),labels=rep("",nrow(modePSDyChlConc)),lwd=2)
box(lwd=2)

points(c(1:nrow(modePSDyChlConc)),modePSDyChlConc$means,cex=2.75,pch=19)
points(c(1:ncol(groupMLChlConc)),groupMLChlConc[6,],cex=3,pch=24,bg="white")

for (i in 1:nrow(modePSDyChlConc)){
  text(i,as.data.frame(QPSDyChlConc[i])[2,2]+0.5,modePSDyChlConc[i,2],cex=2.5)
}

title(ylab=expression(bold("PSD"["NB"]~"(‰)")), line=8, cex.lab=3.5)

mtext("d",side=3,adj=0.01,line=-3,font=2,cex=2.75)

# x-axis tick labels and names for plots c and d.
labels=c("[0.01,0.10[",
         "[0.10,1.00[",
         "[1.00,+∞[")
par(mar=c(0,12,0.5,0.3))
plot(c(0,nrow(modePSDxChlConc)),c(0,1),xlab="",ylab="",
     ax=F,col="white",xaxs="i",xlim=c(1-0.5,nrow(modePSDxChlConc)+0.5))
text(x = 1:nrow(modePSDxChlConc),
     y = 1.1,
     labels =labels,
     xpd = NA,
     srt = 60,
     font=2,
     cex = 3.5,adj=1) 
text(x=nrow(modePSDxChlConc)/2+0.5,0,xpd=NA,expression(bold(atop(
  "Surface chlorophyll"~bolditalic("a")~"concentration","(mg.m"^"-3"*")")))
  ,adj=c(0.5,0),font=2,cex = 3.5)

#End of script.
