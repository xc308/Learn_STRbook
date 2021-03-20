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



#==================#
# Correlated Errors
#==================#

# there is clearly correlation in the residuals
# indicating that the fixed effects are not able to 
# explain the sp-t varaibility in the data


# if the sp_t covariance fun of these errors are known
# can use generalized least sq to fit the model

# if the cov function was a Gau, iso, with range of 0.5

Tmax_J_gls <- gls(z ~ (lon + lat + day)^2 + ., 
    data = select(Tmax2_no_14, -id),
    correlation = corGaus(value = 0.5,
                          form = ~ lon + lat + day,
                          fixed = TRUE))



#===================#
# Stepwise Selection
#===================#

# Stepwise is used to find a parsimounious model
# from a large seletion of explanatory varaible
# s.t. each varialbe is included or excluded in a step
# to minimizes the AIC of the fitted model

# simplest: start from intercept model
# each step is adding/removing one variable
# full model: its scope

# the following for loop retrieves the 
# fitted model for each step of the stepwise AIC forward selection

Tmax_J_lm4 <- list()
Tmax_J_lm4[[i + 1]] <- for(i in 0:4) { # 4 steps after intercept model
  step(lm(z ~ 1, data = select(Tmax2_no_14, -id)),
       scope = z ~ (lon + lat + day) ^ 2 + .,
       direction = "forward",
       steps = i)
  
}

head(Tmax2_no_14, 2)



#===================#
# Multicollinearity
#===================#


# fairly common in sp-t modelling to have multicollinarity
# both in sp and t

# multicollinearity will result in a high correlation
# betw the estimators of the regression coefficients

# the 13th basis function is the 5th basis function corrupted by noise
Tmax2_no_14_2 <- Tmax2_no_14 %>% 
  mutate(B13 = B5 + 0.01 * rnorm(nrow(Tmax2_no_14)))


Tmax_J_lm3 <- lm(z ~ (lon + lat + day) ^ 2 + .,
   data = Tmax2_no_14_2 %>%
     select(-id))

summary(Tmax_J_lm3)
# both 5th and 13th basis functions are no loner
# significant at 1% level. 

# The introduction will not adversely affect 
# the predicitons and prediction std errors.

# but will resutl in a high positive or neg correlation
# btw the estimator of the regr coefficient

vcov(Tmax_J_lm3)[c("B5", "B13"), c("B5", "B13")] %>%
  cov2cor()

head(a)



#====================#
# Analyzing Residuals
#====================#

# Having fitted a sp-t model, good practice to check residuals
# if they are still sp-t correlated, the our model won't have
# captured adequatedly the sp-t variability in the data

# extract the residuals from our linear model
# using function residuals

Tmax_no_14$residuals <- residuals(Tmax_Jul_lm)

# plot the residuals of the last eight days
# show these residuals are spatially correlated


g <- ggplot(filter(Tmax_no_14, day %in% 24:31)) +
  geom_point(aes(lon, lat, colour = residuals)) +
  facet_wrap(~ day, ncol = 4) + 
  col_scale(name = "degF") + 
  geom_point(data = filter(Tmax_no_14, 
                           day %in% 24:31 &
                             id %in% c(3810, 3889)),
             aes(lon, lat), colour = "black",
             pch = 2, size = 2.5) +
  theme_bw()

print(g)


#---------------#
# Moran's I test
#---------------#

# to test for spatial dependence in the residual on each day
# take each day, compute the distances to form the weight matrix
# carry out Moran's I test 

P <- list()
days <- c(1:13, 15:31)
for (i in seq_along(days)) {
  Tmax_day <- filter(Tmax_no_14,
         day == days[i]) 
  
  station.dist <- Tmax_day %>% 
    select(lon, lat) %>%
    dist() %>% 
    as.matrix()
  
  station.dist.inv <- 1 / station.dist
  
  diag(station.dist.inv) <- 0
  
  P[[i]] <- Moran.I(Tmax_day$residuals, 
          station.dist.inv) %>%
    do.call("cbind", .)
}


str(P)
# List of 30
# $ : num [1, 1:4]


# P is a list of single-row data frames
# bind each of these row
do.call("rbind", P) %>% head()


#--------------------------#
# Extending Moran'I to sp-t
#--------------------------#

# usual problem of how to scale time to make
# a Euclidean distance across sp-t and have a 
# realistic interpreation

# one way: fit a dependence model allows for scaling in time
# and scale time by an estimate of the scalling factor
# prior to compute the Euclidean distance
# which uses an anisotropic covariance fun

# for now, did with IDW, don't scale time and 
# compute dis on the sp-t domain

station.dists <- Tmax_no_14 %>%
  select(lon, lat, day) %>%
  dist() %>%
  as.matrix()

station.dists.inv <- 1/station.dists
diag(station.dists.inv) <- 0
Moran.I(Tmax_no_14$residuals, station.dists.inv)$p.value
# 0
# very small, strongly suggest that there's sp-t dependence


# if data are regularly spaced in time
# one may also look at the temporal residuals
# at some locations and test for temporal correlation
# in these residuals using D-W test

TS1 <- filter(Tmax_no_14, id == 3810)$residuals
TS2 <- filter(Tmax_no_14, id == 3889)$residuals


par(mfrow = c(1, 1))
plot(TS1, 
     xlab = "day of Jul 1993",
     ylab = "residuals (degF)",
     type = 'o',
     ylim = c(-8, 7))

lines(TS2, 
      xlab = "day of Jul 1993",
      ylab = "residuals (degF)",
      type = "o", 
      col = "red")

# strong temporal correlation in the residuals
# residuals close to each other in time tend to be
# more simila than residuals further apart

# also, resiuals are correlated betw stations. 
# given a certain time point

acf(TS1) # strong lag 1 correlation
acf(TS2)


#-----------#
# Group data
#-----------#

nested_Tmax_no_14 <- group_by(Tmax_no_14, lon, lat) %>% nest()

nest()

mutate # adds new varaible and preserve exisisting one



#============#
# Prediction
#============#

#----------------#
# prediction grid
#----------------#

pred_grid <- expand.grid(lon = seq(-100, -80, length = 20),
            lat = seq(32, 46, length = 20),
            day = seq(4, 29, length = 6))

str(pred_grid) # 2400 obs. of  3 variables

# would require all covariate values at all prediction locations
# so the 12 basis functions need to be evaluated on this grid
# do this by calling eval_basis and converting the result to matrix
# then attach to our prediction grid

Spred <- eval_basis(basis = G, 
           s = pred_grid[, c("lon", "lat")] %>% # predict location
             as.matrix()) %>%
  as.matrix()

colnames(Spred) <- paste0("B", 1:ncol(Spred))

str(Spred) # num [1:2400, 1:12]

# attach Spred to grid
pred_grid <- cbind(pred_grid, Spred)


linreg_pred <- predict(Tmax_Jul_lm,
        newdata = pred_grid,
        interval = "prediction")

head(linreg_pred)
# the pred interval = pred +/- 1.96 * pred std error


#---------------------------------------#
# Assign pred and pred s.e. to pred grid
#---------------------------------------#

pred_grid$z_pred <- linreg_pred[, 1]
pred_grid$z_err <- (linreg_pred[, 3] - linreg_pred[, 2]) / (2 * 1.96)

head(pred_grid, 2)

g_pred <- ggplot(pred_grid) + 
  geom_tile(aes(lon, lat, fill = z_pred),
            colour = tom.color) + 
  facet_wrap(~ day, ncol = 3) + 
  theme_bw()

print(g_pred)




































































