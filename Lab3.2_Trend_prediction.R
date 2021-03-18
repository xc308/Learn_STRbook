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



























