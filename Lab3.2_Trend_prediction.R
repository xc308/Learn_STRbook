#**************************#
# Lab 3.2 Trend Prediction
#**************************#

# leaps: contain function for stepwise regression
# lmtest: contains a suite of tests to carry out on fitted lm
# nlme: for fitting nonlinear mixed effects models
#   in the presence of correlated errors
# ape: for testing sp or sp-t independence with Moran's I statistic
# FRK: function to constructing the basis functions


install.packages("leaps")
install.packages("lmtest")
install.packages("nlme")
install.packages("ape")
install.packages("broom")
install.packages("FRK")
install.packages("purrr")
install.packages("lattice")
install.packages("dplyr")
install.packages("gstat")
install.packages("sp")
install.packages("spacetime")
install.packages("tidyr")
install.packages("STRbook")



library("leaps")
library("lmtest")
library("nlme")
library("ape")
library("broom")
library("FRK")
library("purrr")

library("lattice")
library("ggplot2")
library("RColorBrewer")

library("dplyr")
library("gstat")
library("sp")
library("spacetime")
library("STRbook")
library("tidyr")


#------------------#
# Basis functions
#------------------#

data("NOAA_df_1990", package = "STRbook")
Tmax <- filter(NOAA_df_1990,
       proc == "Tmax" &
        month == 7 &
         year == 1993)

# set of basis functions can be constructed using auto_basis in FRK
  # data: spatial obj
  # nres: # of resolutions to construct
  # type: type of basis function
  

# example: a single resolution of Gau radial basis functio

G <- auto_basis(data = Tmax[, c("lon", "lat")] %>%
             SpatialPoints(),
           nres = 1, 
           type = "Gaussian")

# These basis functions evaluated at data locations
# are then the covariates we seek for fitting the data

# the functions are evaluated at any arbitrary location using
# the function eval_basis
# which requires the locations as a matrix obj
# and it returns the evaluations as an obj of class Matrix
# which can be converted to matrx

str(Tmax) #  data.frame':	4122 obs. of  10 variables:


S <- eval_basis(basis = G, 
           s = Tmax[, c("lon", "lat")] %>%
             as.matrix()) %>%
  as.matrix()

str(S)  # num [1:4122, 1:12] 
colnames(S) <- paste0("B", 1:ncol(S))

head(S)


#------------------------#
# fitting the linear model
#------------------------#

# 1st remove all variables don't want 
Tmax2 <- cbind(Tmax, S) %>%
  select(-year, -month, -proc, -julian, -date)

# also need to remove 14 July to see how predictions on this day are affected. 
Tmax2_no_14 <- filter(Tmax2, !(day == 14))

# fit linear model using lon, lat, day, their interaction
# along with 12 basis function

# formula: z ~ (lon + lat + day)^2 + .

Tmax_Jul_lm <- lm(z ~ (lon + lat + day)^2 + .,
   data = select(Tmax2_no_14, -id))


# view result in summary
Tmax_Jul_lm %>% summary()

# lon:day     -17.615  < 2e-16 ***
# lat:day     -10.147  < 2e-16 ***



























