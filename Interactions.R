# This script focuses on the investigation of the impacts of the interactions 
# between environmental impacts on the stable isotope values of sea stars in the 
# Southern Ocean.

################################################################################
# Packages
################################################################################

# Let's start by loading the necessary packages.
library(dplyr) # For dataset combination
library(Rcmdr) # Contains the logit function
library(viridis) # Contains the "turbo" colour palette used for the plots.

################################################################################
# Functions
################################################################################

# We need to code some functions that will be used in the subsequent analysis.

# The first one is inv.logit, to convert back logit-transformed values to their
# original values. It is necessary to assess and show the influence of sea ice
# concentration on stable isotope data.
inv.logit <- function(f,a) {
  a <- (1-2*a)
  (a*(1+exp(f))+(exp(f)-1))/(2*a*(1+exp(f)))
}

# The second one is predictgrid. Given a model, this function predicts the values
# of the z variable from the x and y variables. Defaults to range of x and y 
# variables, and a 16x16 grid.
predictgrid <- function(model, xvar, yvar, zvar, res = 16, type = NULL) {
  xrange <- range(mod$model[,2])
  yrange <- range(mod$model[,3])
  newdata <- expand.grid(x = seq(xrange[1], xrange[2], length.out = res),
                         y = seq(yrange[1], yrange[2], length.out = res))
  names(newdata) <- c(xvar, yvar)
  newdata[[zvar]] <- predict(model, newdata = newdata, type = type)
  newdata
}

# The third one is df2mat. It converts long-style data frame with x, y,
# and z variables into a list with x and y as row/column values, and z as a matrix.
df2mat <- function(p, xvar = NULL, yvar = NULL, zvar = NULL) {
  if (is.null(xvar)) xvar <- names(p)[1]
  if (is.null(yvar)) yvar <- names(p)[2]
  if (is.null(zvar)) zvar <- names(p)[3]
  
  x <- unique(p[[xvar]])
  y <- unique(p[[yvar]])
  z <- matrix(p[[zvar]], nrow = length(y), ncol = length(x))
  
  m <- list(x, y, z)
  names(m) <- c(xvar, yvar, zvar)
  m
}

# The fourth one is a modification of the filled.contour function. It has been 
# created to obtain more aesthetically pleasing and consistent plots (by modifying
# default style of the axes) and to remove the border in the legend in prevision
# of the plot assemblage.
filled.contour<-function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
                          length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
                          ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
                          levels = pretty(zlim, nlevels), nlevels = 20, color.palette = function(n) hcl.colors(n, 
                          "YlOrRd", rev = TRUE), col = color.palette(length(levels) - 
                          1), plot.title, plot.axes, key.title, key.axes, asp = NA, 
                          xaxs = "r", yaxs = "r", las = 1, axes = TRUE, frame.plot = axes, 
                          ...) 
{
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
  layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
  par(las = las)
  mar <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  par(mar = mar)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "r", 
              yaxs = "r")
  rect(0, levels[-length(levels)], 1, levels[-1L], col = col, border = NA)
  if (missing(key.axes)) {
    if (axes) 
      axis(4)
  }
  else key.axes
  box()
  if (!missing(key.title)) 
    key.title
  mar <- mar.orig
  mar[4L] <- 1
  par(mar = mar)
  plot.new()
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  .filled.contour(x, y, z, levels, col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}

################################################################################
# Data loading and preparation
################################################################################

# Let's load and prepare the stable isotope and environmental data…
DataSeaStar<-read.table("IsotopeData.csv",dec=".",sep=",",header=T)
DataEnv<-read.table("EnvironmentalData.csv",dec=".",sep=",",header=T)

# First we combine both datasets.
Data <- DataSeaStar %>% left_join(DataEnv, by = c("StationID","ExpeditionID","Date"))

# We need to log-transform depth (log 10) and chlorophyll concentration (natural
# log) data prior to doing the analysis. However, sea ice concentration will have 
# to be logit transformed within the models.
Data$logDepth<-log10(Data$Depth)
Data$logChl<-log(Data$chl_prev_month)


# For models with sea ice concentration and sea ice duration as explanatory 
# variables, we need a dataframe where individuals from Subantarctic regions are
# removed.  These individuals were collected in areas where sea ice has not been
# present since more than 1000 days.
DataAnt <- Data[-which(Data$days_since_melt > 1000),]

################################################################################
# Data analysis and plots
################################################################################

# The process of data analysis is the same for each combination of environmental
# variable: we start by creating a linear model assessing the impact of the two 
# environmental variables and their interactions. We then predict the stable 
# isotope values with the predicgrid function. We convert the result into a list 
# with the function df2mat and then use the result to make the plot.

# for models using sea ice concentration as an explanatory variable, we need to 
# logit-transform sea ice concentration values in the models and convert them 
# back to raw values before doing the prediction.

# The plots were saved as tiff to be put togeether in open office later. Axis 
# tick labels were also added manually with open office. Consequenly, no tick 
# labels were included on the plots.

# Let's start with the influence of the interaction between depth and sea ice 
# concentration on δ13C values. This first example is annotated to re-explain the 
# process.

##### Creation of the model. Notice that sea ice concentration is logit-transformed
##### in the model.
mod <- lm(d13C ~ logDepth + logit(seaice_prev_month) + 
            +logDepth:logit(seaice_prev_month), data = DataAnt)

##### Convert back logit-transformed sea ice concentration to raw values.
mod$model[,3]<-inv.logit(mod$model[,3],a=0.025)*100

##### Prediction of δ13C values.
mpgrid_df <- predictgrid(mod, "logDepth", "seaice_prev_month", "d13C")

##### Converion of the results into a list.
mpgrid_list <- df2mat(mpgrid_df)

##### We checked the range of predicted stable isotope values for each model in 
##### as we want to display the same range of predicted δ13C and δ15N values on
##### all plots.
range(mpgrid_list$d13C)

##### Plot creation.
dev.new(width=17, height=16, unit="cm")
filled.contour(mpgrid_list$logDepth, mpgrid_list$seaice_prev_month, mpgrid_list$d13C,
               zlim=c(-27,-10),xlim=c(log10(5),log10(5500)),ylim=c(0,90),color.palette = turbo,nlevels = 170,
               plot.axes = {
                 axis(1,at=c(log(5,10),log(50,10),log(500,10),log(5000,10)),labels=rep("",4),lwd=5,tck=-0.03)
                 axis(2,at=c(0,30,60,90),labels=rep("",4),lwd=5,tck=-0.03)
                 box(lwd=5)
                 points(DataAnt$logDepth, DataAnt$seaice_prev_month,pch=21,bg="white",lwd=2.5,cex=5)
                 contour(mpgrid_list$logDepth, mpgrid_list$seaice_prev_month,mpgrid_list$d13C, add = TRUE,vfont=c("sans serif", "bold"),labcex=3.5, lwd = 10,levels=seq(-27,-10))
               })

# now we can do the same thing for other interactions and both δ13C and δ15N values.

# Influence of the interaction between depth and sea ice concentration on δ15N values.
mod <- lm(d15N ~ logDepth + logit(seaice_prev_month) + 
            +logDepth:logit(seaice_prev_month), data = DataAnt)

mod$model[,3]<-inv.logit(mod$model[,3],a=0.025)*100

mpgrid_df <- predictgrid(mod, "logDepth", "seaice_prev_month", "d15N")
mpgrid_list <- df2mat(mpgrid_df)

range(mpgrid_list$d15N)

dev.new(width=17, height=16, unit="cm")
filled.contour(mpgrid_list$logDepth, mpgrid_list$seaice_prev_month, mpgrid_list$d15N, frame = F,
               zlim=c(6,18),xlim=c(log10(5),log10(5500)),ylim=c(0,90),color.palette = turbo,nlevels = 120,
               plot.axes = {
                 axis(1,at=c(log(5,10),log(50,10),log(500,10),log(5000,10)),labels=rep("",4),lwd=5,tck=-0.03)
                 axis(2,at=c(0,30,60,90),labels=rep("",4),lwd=5,tck=-0.03)
                 box(lwd=5)
                 points(DataAnt$logDepth, DataAnt$seaice_prev_month,pch=21,bg="white",lwd=2.5,cex=5)
                 contour(mpgrid_list$logDepth, mpgrid_list$seaice_prev_month,mpgrid_list$d15N, add = TRUE,vfont=c("sans serif", "bold"),labcex=3.5, lwd = 10,levels=seq(6,18))
               })

# Influence of the interaction between depth and sea ice duration on δ13C values.
mod <- lm(d13C ~ logDepth + seaice_last_730 + 
            +logDepth:seaice_last_730, data = DataAnt)

mpgrid_df <- predictgrid(mod, "logDepth", "seaice_last_730", "d13C")
mpgrid_list <- df2mat(mpgrid_df)

range(mpgrid_list$d13C)

dev.new(width=17, height=16, unit="cm")
filled.contour(mpgrid_list$logDepth, mpgrid_list$seaice_last_730, mpgrid_list$d13C, frame = F,
               zlim=c(-27,-10),xlim=c(log10(5),log10(5500)),ylim=c(0,730),color.palette = turbo,nlevels = 170,
               plot.axes = {
                 axis(1,at=c(log(5,10),log(50,10),log(500,10),log(5000,10)),labels=rep("",4),lwd=5,tck=-0.03)
                 axis(2,at=c(0,146,292,438,584,730),labels=rep("",6),lwd=5,tck=-0.03)
                 box(lwd=5)
                 points(DataAnt$logDepth, DataAnt$seaice_last_730,pch=21,bg="white",lwd=2.5,cex=5)
                 contour(mpgrid_list$logDepth, mpgrid_list$seaice_last_730,mpgrid_list$d13C, add = TRUE,vfont=c("sans serif", "bold"),labcex=3.5, lwd = 10,levels = seq(-27,-10))
               })

# Influence of the interaction between depth and sea ice duration on δ15N values.
mod <- lm(d15N ~ logDepth + seaice_last_730 + 
            +logDepth:seaice_last_730, data = DataAnt)


mpgrid_df <- predictgrid(mod, "logDepth", "seaice_last_730", "d15N")
mpgrid_list <- df2mat(mpgrid_df)

range(mpgrid_list$d15N)

dev.new(width=17, height=16, unit="cm")
filled.contour(mpgrid_list$logDepth, mpgrid_list$seaice_last_730, mpgrid_list$d15N, frame = F,
               zlim=c(6,18),xlim=c(log10(5),log10(5500)),ylim=c(0,730),color.palette = turbo,nlevels = 120,
               plot.axes = {
                 axis(1,at=c(log(5,10),log(50,10),log(500,10),log(5000,10)),labels=rep("",4),lwd=5,tck=-0.03)
                 axis(2,at=c(0,146,292,438,584,730),labels=rep("",6),lwd=5,tck=-0.03)
                 box(lwd=5)
                 points(DataAnt$logDepth, DataAnt$seaice_last_730,pch=21,bg="white",lwd=2.5,cex=5)
                 contour(mpgrid_list$logDepth, mpgrid_list$seaice_last_730,mpgrid_list$d15N, add = TRUE,vfont=c("sans serif", "bold"),labcex=3.5, lwd = 10,levels = seq(9,13))
               })

# Influence of the interaction between depth and chlorophyll concentration on δ13C 
# values.
mod <- lm(d13C ~ logDepth + logChl + 
            +logDepth:logChl, data = Data)

mpgrid_df <- predictgrid(mod, "logDepth", "logChl", "d13C")
mpgrid_list <- df2mat(mpgrid_df)

range(mpgrid_list$d13C)

dev.new(width=17, height=16, unit="cm")
filled.contour(mpgrid_list$logDepth, mpgrid_list$logChl, mpgrid_list$d13C,
               zlim=c(-27,-10),xlim=c(log10(5),log10(5500)),ylim=c(log(0.02),log(10)),color.palette = turbo,nlevels = 170,ax=F,
               plot.axes = {
                 axis(1,at=c(log(5,10),log(50,10),log(500,10),log(5000,10)),labels=rep("",4),lwd=5,tck=-0.03)
                 axis(2,at=c(log(0.02),log(0.1),log(1),log(10)), labels=rep("",4),lwd=5,tck=-0.03)
                 box(lwd=5)
                 points(Data$logDepth, Data$logChl,pch=21,bg="white",lwd=2.5,cex=5)
                 contour(mpgrid_list$logDepth, mpgrid_list$logChl, mpgrid_list$d13C, add = TRUE,vfont=c("sans serif", "bold"),labcex=3.5, lwd = 10,levels=seq(-27,-10))
               })

# Influence of the interaction between depth and chlorophyll concentration on δ15N 
# values.
mod <- lm(d15N ~ logDepth + logChl + 
            +logDepth:logChl, data = Data)

mpgrid_df <- predictgrid(mod, "logDepth", "logChl", "d15N")
mpgrid_list <- df2mat(mpgrid_df)

range(mpgrid_list$d15N)

dev.new(width=17, height=16, unit="cm")
filled.contour(mpgrid_list$logDepth, mpgrid_list$logChl, mpgrid_list$d15N,
               zlim=c(6,18),xlim=c(log10(5),log10(5500)),ylim=c(log(0.02),log(10)),color.palette = turbo,nlevels = 120,ax=F,
               plot.axes = {
                 axis(1,at=c(log(5,10),log(50,10),log(500,10),log(5000,10)),labels=rep("",4),lwd=5,tck=-0.03)
                 axis(2,at=c(log(0.02),log(0.1),log(1),log(10)), labels=rep("",4),lwd=5,tck=-0.03)
                 box(lwd=5)
                 points(Data$logDepth, Data$logChl,pch=21,bg="white",lwd=2.5,cex=5)
                 contour(mpgrid_list$logDepth, mpgrid_list$logChl, mpgrid_list$d15N, add = TRUE,vfont=c("sans serif", "bold"),labcex=3.5, lwd = 10,levels=seq(6,18))
               })

# Influence of the interaction between sea ice concentration and sea ice duration
# on δ13C values.
mod <- lm(d13C ~ logit(seaice_prev_month) + seaice_last_730 + 
            +logit(seaice_prev_month):seaice_last_730, data = DataAnt)


mod$model[,2]<-inv.logit(mod$model[,2],a=0.025)*100

mpgrid_df <- predictgrid(mod, "seaice_prev_month", "seaice_last_730", "d13C")
mpgrid_list <- df2mat(mpgrid_df)

range(mpgrid_list$d13C)

dev.new(width=17, height=16, unit="cm")
filled.contour(mpgrid_list$seaice_prev_month, mpgrid_list$seaice_last_730, mpgrid_list$d13C,
               zlim=c(-27,-10),xlim=c(0,90),ylim=c(0,730),color.palette = turbo,nlevels = 170,ax=F,
               plot.axes = {
                 axis(1,at=c(0,30,60,90),labels=rep("",4),lwd=5,tck=-0.03)
                 axis(2,at=c(0,146,292,438,584,730), labels=rep("",6),lwd=5,tck=-0.03)
                 box(lwd=5)
                 points(DataAnt$seaice_prev_month, DataAnt$seaice_last_730,pch=21,bg="white",lwd=2.5,cex=5)
                 contour(mpgrid_list$seaice_prev_month, mpgrid_list$seaice_last_730, mpgrid_list$d13C, add = TRUE,vfont=c("sans serif", "bold"),labcex=3.5, lwd = 10,levels=seq(-27,-10))
               })

# Influence of the interaction between sea ice concentration and sea ice duration
# on δ15N values.
mod <- lm(d15N ~ logit(seaice_prev_month) + seaice_last_730 + 
            +logit(seaice_prev_month):seaice_last_730, data = DataAnt)

mod$model[,2]<-inv.logit(mod$model[,2],a=0.025)*100

mpgrid_df <- predictgrid(mod, "seaice_prev_month", "seaice_last_730", "d13C")
mpgrid_list <- df2mat(mpgrid_df)

range(mpgrid_list$d13C)

dev.new(width=17, height=16, unit="cm")
filled.contour(mpgrid_list$seaice_prev_month, mpgrid_list$seaice_last_730, mpgrid_list$d13C,
               zlim=c(6,18),xlim=c(0,90),ylim=c(0,730),color.palette = turbo,nlevels = 120,ax=F,
               plot.axes = {
                 axis(1,at=c(0,30,60,90),labels=rep("",4),lwd=5,tck=-0.03)
                 axis(2,at=c(0,146,292,438,584,730), labels=rep("",6),lwd=5,tck=-0.03)
                 box(lwd=5)
                 points(DataAnt$seaice_prev_month, DataAnt$seaice_last_730,pch=21,bg="white",lwd=2.5,cex=5)
                 contour(mpgrid_list$seaice_prev_month, mpgrid_list$seaice_last_730, mpgrid_list$d13C, add = TRUE,vfont=c("sans serif", "bold"),labcex=3.5, lwd = 10,levels=seq(6,18))
               })

# Influence of the interaction between sea ice concentration and chlorophyll 
# concentration on δ13C values.
mod <- lm(d13C ~ logit(seaice_prev_month) + logChl + 
            +logit(seaice_prev_month):logChl, data = DataAnt)

mod$model[,2]<-inv.logit(mod$model[,2],a=0.025)*100

mpgrid_df <- predictgrid(mod, "seaice_prev_month", "logChl", "d13C")
mpgrid_list <- df2mat(mpgrid_df)

range(mpgrid_list$d13C)
range(mpgrid_list$seaice_prev_month)

dev.new(width=17, height=16, unit="cm")
filled.contour(mpgrid_list$seaice_prev_month, mpgrid_list$logChl, mpgrid_list$d13C,
               zlim=c(-27,-10),xlim=c(0,90),ylim=c(log(0.02),log(10)),color.palette = turbo,nlevels = 170,ax=F,
               plot.axes = {
                 axis(1,at=c(0,30,60,90),labels=rep("",4),lwd=5,tck=-0.03)
                 axis(2,at=c(log(0.02),log(0.1),log(1),log(10)), labels=rep("",4),lwd=5,tck=-0.03)
                 box(lwd=5)
                 points(DataAnt$seaice_prev_month, DataAnt$logChl,pch=21,bg="white",lwd=2.5,cex=5)
                 contour(mpgrid_list$seaice_prev_month, mpgrid_list$logChl, mpgrid_list$d13C, add = TRUE,vfont=c("sans serif", "bold"),labcex=3.5, lwd = 10,levels=seq(-27,-10))
               })

# Influence of the interaction between sea ice concentration and chlorophyll 
# concentration on δ15N values.
mod <- lm(d15N ~ logit(seaice_prev_month) + logChl + 
            +logit(seaice_prev_month):logChl, data = DataAnt)

mod$model[,2]<-inv.logit(mod$model[,2],a=0.025)*100

mpgrid_df <- predictgrid(mod, "seaice_prev_month", "logChl", "d15N")
mpgrid_list <- df2mat(mpgrid_df)

range(mpgrid_list$d15N)

dev.new(width=17, height=16, unit="cm")
filled.contour(mpgrid_list$seaice_prev_month, mpgrid_list$logChl, mpgrid_list$d15N,
               zlim=c(6,18),xlim=c(0,90),ylim=c(log(0.02),log(10)),color.palette = turbo,nlevels = 120,ax=F,
               plot.axes = {
                 axis(1,at=c(0,30,60,90),labels=rep("",4),lwd=5,tck=-0.03)
                 axis(2,at=c(log(0.02),log(0.1),log(1),log(10)), labels=rep("",4),lwd=5,tck=-0.03)
                 box(lwd=5)
                 points(DataAnt$seaice_prev_month, DataAnt$logChl,pch=21,bg="white",lwd=2.5,cex=5)
                 contour(mpgrid_list$seaice_prev_month, mpgrid_list$logChl, mpgrid_list$d15N, add = TRUE,vfont=c("sans serif", "bold"),labcex=3.5, lwd = 10,levels=seq(6,18))
               })

# Influence of the interaction between sea ice duration and chlorophyll concentration
# on δ13C values.
mod <- lm(d13C ~ seaice_last_730 + logChl + 
            +seaice_last_730:logChl, data = DataAnt)

mpgrid_df <- predictgrid(mod, "seaice_last_730", "logChl", "d13C")
mpgrid_list <- df2mat(mpgrid_df)

range(mpgrid_list$d13C)
dev.new(width=17, height=16, unit="cm")
filled.contour(mpgrid_list$seaice_last_730, mpgrid_list$logChl, mpgrid_list$d13C,
               zlim=c(-27,-10),xlim=c(0,730),ylim=c(log(0.02),log(10)),color.palette = turbo,nlevels = 170,ax=F,
               plot.axes = {
                 axis(1,at=c(0,146,292,438,584,730),labels=rep("",6),lwd=5,tck=-0.03)
                 axis(2,at=c(log(0.02),log(0.1),log(1),log(10)), labels=rep("",4),lwd=5,tck=-0.03)
                 box(lwd=5)
                 points(DataAnt$seaice_last_730, DataAnt$logChl,pch=21,bg="white",lwd=2.5,cex=5)
                 contour(mpgrid_list$seaice_last_730, mpgrid_list$logChl, mpgrid_list$d13C, add = TRUE,vfont=c("sans serif", "bold"),labcex=3.5, lwd = 10,levels=seq(-27,-10))
               })

# Influence of the interaction between sea ice duration and chlorophyll concentration
# on δ15N values.
mod <- lm(d15N ~ seaice_last_730 + logChl + 
            +seaice_last_730:logChl, data = DataAnt)

mpgrid_df <- predictgrid(mod, "seaice_last_730", "logChl", "d15N")
mpgrid_list <- df2mat(mpgrid_df)

range(mpgrid_list$d15N)

dev.new(width=17, height=16, unit="cm")
filled.contour(mpgrid_list$seaice_last_730, mpgrid_list$logChl, mpgrid_list$d15N,
               zlim=c(6,18),xlim=c(0,730),ylim=c(log(0.02),log(10)),color.palette = turbo,nlevels = 120,ax=F,
               plot.axes = {
                 axis(1,at=c(0,146,292,438,584,730),labels=rep("",6),lwd=5,tck=-0.03)
                 axis(2,at=c(log(0.02),log(0.1),log(1),log(10)), labels=rep("",4),lwd=5,tck=-0.03)
                 box(lwd=5)
                 points(DataAnt$seaice_last_730, DataAnt$logChl,pch=21,bg="white",lwd=2.5,cex=5)
                 contour(mpgrid_list$seaice_last_730, mpgrid_list$logChl, mpgrid_list$d15N, add = TRUE,vfont=c("sans serif", "bold"),labcex=3.5, lwd = 10,levels=seq(6,18))
               })

# End of script.
