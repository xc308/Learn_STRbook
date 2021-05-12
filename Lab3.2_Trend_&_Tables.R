#*************************#
# lab 3.2 Trend Prediction
#*************************#

# leap: for stepwise regression
# lmtest: contains a suite of tests on fitted linear model
# nlme: nonlinear mixed model
      # but we use it to fit linear model in the presence of correlated errors

# ape: test sp or sp-t independence with Moran's I
# FRK: constructing basis functions

# broom and purrr: carry out multiple test on groups of dataset

library("leaps")
library("lmtest")
library("nlme")

library("ape")
library("broom")
library("FRK")
library('purrr')

library("lattice")
library("ggplot2")
library("RColorBrewer")


# for data wrangling
library("dplyr")
library("gstat")
library("sp")
library("spacetime")
library("STRbook")
library("tidyr")

install.packages("stargazer")
library(stargazer)


#======#
# Data
#======#

data("NOAA_df_1990", package = "STRbook")
Tmax <- NOAA_df_1990 %>% filter(proc == "Tmax" &
                          month == 7 &
                          year == 1993)

head(Tmax)
str(Tmax)
# 'data.frame':	4122 obs. of  10 variables:


# construct basis function 
G <- auto_basis(data = Tmax[, c("lon", "lat")] %>% SpatialPoints(),
           nres = 1,
           type = "Gaussian")

show_basis(G)

# evaluate the basis at locations will then be
# the covariates seek for fitting the data
  # require locations to be matrix
  # return obj of class Matrix
  # which can be easily converted to a matrix

S <- eval_basis(basis = G, 
           s = Tmax[, c("lon", "lat")] %>% as.matrix()) %>%
  as.matrix()

str(S)

colnames(S) <- paste0("B", 1:ncol(S))
head(S)


#===============#
# fit the model
#===============#

# use . to denote all variables in the df as covariates
# so append basis to the df, and rmv unwanted variables

Tmax2 <- cbind(Tmax, S) %>%
  select(-year, -month, -proc, -julian, -date)

str(Tmax2)

# also rmv 14 July 1993 to see how predictions on this day is affected 
# given no obs on this day

Tmax_no_14 <- filter(Tmax2, !(day == 14))


# (lon + lat + day)^2 : each of the var and their interactions

Tmax_July_lm <- lm(z ~ (lon + lat + day)^2 + .,
   data = select(Tmax_no_14, -id))

Tmax_July_lm %>% summary()


#------------------#
# Correlated Errors
#------------------#

# if there are correlation in the residuals
# then fixed effects are not able to explain 
# the sp-t variablity in the data

# and if we knew the sp covariance fucntion of these errors
# we could use generalized least squares to fit the model

  # cov is Gaussian, isotropic, range 0.5

install.packages("nlme")
library(nlme)

Tmax_July_gls <- gls(z ~ (lon + lat + day)^2 + .,
    data = select(Tmax_no_14, -id),
    correlation = corGaus(value = 0.5,
                          form = ~ lon + lat + day,
                          fixed = T))

Tmax_July_lm %>% summary()


#==================#
# Stepwise Selction
#==================#

# used t find a parsimounious model

# step()
  # start from intercept model 
  # scope: full model
  # direction: forward
  # each step introduce a new variable that minimizes tha AIC of the fitted model
  # which penalize complexity

Tmax_July_lm4 <- list()
for (i in 0:4){
  Tmax_July_lm4[[i + 1]] <- step(lm(z ~ 1, data = select(Tmax_no_14, -id)),
       scope = z ~ (lon + lat + day)^2 + .,
       direction = 'forward',
       steps = i)
}


stargazer(Tmax_July_lm4[[1]], Tmax_July_lm4[[2]], Tmax_July_lm4[[3]], Tmax_July_lm4[[4]], Tmax_July_lm4[[5]],
          title = "Stepwise Results",
          align = T, 
          omit.stat = c("LL", "ser", "f"),
          no.space = T,
          covariate.labels = c("Lat", "day", "Lon", "Lat:day", "Intercept"),
          order = c("Intercept", "Lat", "day", "Lat:day"))




summary_lm4 <- list()
for(i in 1:5) {
  summary_lm4[[i]] <- Tmax_July_lm4[[i]] %>% summary()
}




str(Tmax_July_lm4)
#-----------#
# regsubsets
#-----------#

# select covariates whose effect are significant 
# when minimizing the residual sum of squares at each step

regfit.full <- regsubsets(z ~ 1 + (lon + lat + day)^2 + .,
           data = select(Tmax_no_14, -id),
           method = "forward",
           nvmax = 4)


regfit.summary <- summary(regfit.full)
regfit.summary



#=================#
# Multicollinarity
#=================#

str(Tmax_no_14)
# 'data.frame':	3989 obs. of  17 variables:
# $ B5 : num  0.0062 0.0062 0.0062 0.0062 0.0062


Tmax_no_14_2 <- Tmax_no_14 %>% mutate(B13 = B5 + 0.01 * rnorm(nrow(Tmax_no_14)))

Tmax_July_lm3 <- lm(z ~ (lon + lat + day)^2 + .,
   data = Tmax_no_14_2 %>% select(-id))


install.packages("stargazer")
library(stargazer)


stargazer(Tmax_July_lm, Tmax_July_lm3, title = "Results", align = T,
          no.space = T)





A <-matrix(c(2, 4, 5, 6), nrow = 2)

t(A) %*% A

crossprod(A)









