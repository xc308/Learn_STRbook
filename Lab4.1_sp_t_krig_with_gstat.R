#********************************#
# Lab 4.1 Sp-t Kriging with gstat
#********************************#

# process of carryign out sp-t universal kriging
# using semivariogram with 

library("sp")
library("spacetime")
library("ggplot2")
library("dplyr")
library("gstat")
library("RColorBrewer")
library("STRbook")
library("tidyr")


# ST-kriging of the max-temp data set in july 1993
# need to fit a parametric function to the empirical 
# semivariogram vv 

data("STObj3", package = "STRbook")
STOjb4 <- STObj3[, "1993-07-01::1993-07-31"]
str(STOjb4)

vv <- variogram(object = z ~ 1 + lat,
          data = STOjb4,
          width = 80,
          cutoff = 1000,
          tlags = 0.01:6.01)

# need to fit a parametric function to the vv
# a number of covariance function models are availabe
# in gstat
# see
vignette("spatio-temporal-kriging")

# so the semivariogram corresponds to sp-t seperable covariance model
# to construct the sp-t semivariogram using gstat
  # function: vgmST
  # stModel = "seperable"
  # vgm: construct individual semivariograms for sp, t
      # partial sill, model type, range, the nugget
  # sill: defines the joint sp-t sill
  # nubmers used are initial values

  # fit.StVariogram: fit parameter semivariog to vv 

sepVgm <- vgmST(stModel = "separable",
      space = vgm(10, "Exp", 400, nugget = 0.1),
      time = vgm(10, "Exp", 1, nugget = 0.1),
      sill = 20)

sepVgm <- fit.StVariogram(vv, sepVgm)


# 2nd model we fit assume identical spatial and temporal
# covariance functions except for spatio-temporal anisotropy, 
# allows to use a spatio-temporal metric covariance model 


metricVgm <- vgmST(stModel = "metric",
      joint = vgm(100, "Exp", 400, nugget = 0.1),
      sill = 10,
      stAni = 100) # sp field is evolving on scales of the order of hundreds of km 

metricVgm <- fit.StVariogram(vv, metricVgm)

# compare the fits of the two semivariograms 
# by checking the mean squared error of the fits

metricMSE <- attr(metricVgm, "optim")$value
# 2.09
sepMSE <- attr(sepVgm, "optim")$value
# 1.42

# so separable model is better fit to the empirical semivariogram

plot(vv, list(sepVgm, metricVgm), main = "Semi-variance")
#str(sepVgm)


#-----------------------------------------------#
# use fitted ST covariance models for prediction
#-----------------------------------------------#

# using ST kriging

# 1st, create prediction grid
  # consider location betw 100 W to 80W, 32N to 46N
  # when coverting to SpatialPoints, ensure CRS of the
  # prediction is the same as that of the obs

spat_pred_grid <- expand.grid(lon = seq(-100, -80, length = 20),
            lat = seq(32, 46, length = 20)) %>%
  SpatialPoints(proj4string = CRS(proj4string(STObj3)))

gridded(spat_pred_grid) <- TRUE


# for temporal grid, consider 6 equally space day in Jul 1993

temp_pred_grid <- as.Date("1993-07-01") + seq(3, 28, length = 6)

# now combine sp-t grid to construct an STF obj for prediction
DE_pred <- STF(sp = spat_pred_grid, 
    time = temp_pred_grid)



# since there are missing obs in STObj4
# 1st need to turn STObj4 into STIDF (faster than STSDF)
# to show the capability of ST kriging to predict across time
# omitted data on 14 Jul 1994
STObj5 <- as(STOjb4[, -14], "STIDF")
str(STObj5)

STObj5 <- subset(STObj5, !is.na(STObj5$z))
#subset(sth, condition)

## now prediction on prediction grid using STObj5
pred_kriged <- krigeST(z ~ 1 + lat, 
       data = STObj5, # data set w/o 14 Jul
       newdata = DE_pred,  # prediction grid
       modelList = sepVgm, # semivariogram
       computeVar = TRUE) # compute var


#-------------------------------------#
# Plot predictions and prediction s.e
#-------------------------------------#

color_pal <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(16))

stplot(pred_kriged, 
       main = "Predictions (deg F)",
       layout = c(3, 2),
       col.regions = color_pal)

str(pred_kriged)
pred_kriged$se <- sqrt(pred_kriged$var1.var)

# plot
stplot(pred_kriged[, ,"se"],
       main = "Prediction std. errors (deg F)",
       layout = c(3, 2),
       col.regions = color_pal)


#--------#
# Comment
#--------#

# sp-t kriging is quick and easy to implement for 
# small data sets
# but it starts to become prohibitive as data sets grow in size
# unless approximation is used. 

# krigeST uses nmax to determine the max # of obs to do prediction
# the predictor is no longer optimal, but close enough to the optimal predictor

















