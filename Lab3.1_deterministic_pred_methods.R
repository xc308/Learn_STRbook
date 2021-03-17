#****************************************#
# Lab 3.1 Deterministic Prediction Methods
#****************************************#

#===========================#
# Inverse Distance Weighting
#===========================#

# IDW: one of the simplest determininstic sp-t
# interpolation methods

# can be implemented easily using idw function in gstat

install.packages("dplyr")
install.packages("fields")
installed.packages("ggplot2")
install.packages("gstat")
install.packages("RColorBrewer")
install.packages("sp")
install.packages("spacetime")
install.packages("STRbook")


library(dplyr)
library(fields)
library(ggplot2)
library(gstat)
library(RColorBrewer)
library(sp)
library(spacetime)
library(STRbook)

data("NOAA_df_1990", package = "STRbook")
str(NOAA_df_1990)

Tmax <- filter(NOAA_df_1990, 
       proc == "Tmax" &
       month == 7 & year == 1993)


#------------------------------------------------------#
# Construct 3-d sp-t prediction grid using expand.grid
#------------------------------------------------------#

pre_grid <- expand.grid(lon = seq(-100, -80, length = 20),
            lat = seq(32, 46, length = 20),
            day = seq(4, 29, length = 6))


# function idw arg:
  # formula: identifies the variable to interpolate
  # loactions: spatial temporal varaibles
  # data: data in the df
  # newdata: sp-t grid locations at which to interpolate
  # idp: alpha, the larger, the less smoothing
  # which is set using CV, here set alpha = 5


# remove 14 July 1993
Tmax_14 <- filter(Tmax, day == 14)
str(Tmax_14)
is.na(Tmax_14$z )


Tmax_no_14 <- filter(Tmax, !(day == 14))

Tmax_Jul_idw <- idw(formula = z ~ 1, #dep.var
    locations = ~ lon + lat + day, # inputs
    data = Tmax_no_14,
    newdata = pre_grid,
    idp = 5)

str(Tmax_Jul_idw)
# var1.pred: the IDW interpolation over the prediction grid
# this df can be plotted using ggplot2

ggplot(Tmax_Jul_idw) +
  geom_tile(aes(x = lon, y = lat, fill = var1.pred)) + 
  fill_scale(name = "degF") + 
  xlab("Longitude (deg)") + 
  ylab("Latitude (deg)") + 
  facet_wrap(~ day, ncol = 3) + 
  coord_fixed(xlim = c(-100, -80),
              ylim = c(32, 46)) + 
  theme_bw()



# try without the missing 14 Jul
Tmax_Jul_idw_2 <- idw(formula = z ~ 1, #dep.var
                    locations = ~ lon + lat + day, # inputs
                    data = Tmax,
                    newdata = pre_grid,
                    idp = 5)

ggplot(Tmax_Jul_idw_2) + 
  geom_tile(aes(x = lon, y = lat, fill = var1.pred)) +
  fill_scale(name = "degF") + 
  facet_wrap(~ day, ncol = 3) + 
  coord_fixed(xlim = c(-100, -80),
              ylim = c(32, 46)) + 
  theme_bw()

# so could compare how the IDW is smoothed the missing 14 July


#=======================================#
# Implementing IDW from First Principles
#=======================================#

# IDW from scartch, code versatility for CV
# IDW is a kernel predictor and yields the kernel weight

# To construct kernel weights, 
  # 1st find the dist between all prediction locations and data locations
  # 2nd take reciprocals and raise them to power idp

# if wish to gernerate kernel weights for different obs and prediction sets
# and different bandwidth parameters, 
# we create Wt_IDW that generate the required kernel-weights matrix

pred_obs_dist_mat <- rdist(select(pre_grid, lon, lat, day), 
      select(Tmax_no_14, lon, lat, day))

Wt_IDW <- function(theta, dist_mat) 1/dist_mat ^ theta
Wtilde <- Wt_IDW(theta = 5, dist_mat = pred_obs_dist_mat)
# the (k, l) th element in Wtilde contains the 
# distance btw the kth prediction lication and the 
# lth obs location, raised to power5, then reciprocated. 


#----------------------------#
# Kernel weights normalizing 
#----------------------------#

# by sum of all kernel weights associated with each prediction location
# normalization every location using rowSums
Wtilde_rsums <- rowSums(Wtilde)
W <- Wtilde/Wtilde_rsums
# the desired weight matrix, known as influence matrix

# So the predicitons is just the influence matrix multiply by data
z_pred_IDW <- as.numeric(W %*% Tmax_no_14$z)


summary(Tmax_Jul_idw$var1.pred - z_pred_IDW)



#=================================#
# Generic Kernel Smoothing and CV
#=================================#

# implementing IDW from scratch is now we can 
# change the kernel function to whatever we want
# and compare predictions from different kernel function


# Gauss kernel
Wt_Gauss <- function(theta, dist_mat) exp(-dist_mat ^ 2 / theta)
 
Wtilde_G <- Wt_Gauss(theta = 0.5, dist_mat = pred_obs_dist_mat)

Wtilde_G_rsums <- rowSums(Wtilde_G)

W <- Wtilde_G / Wtilde_G_rsums

z_pred2 <- W %*% Tmax_no_14$z



#================================#
# Which predictions are the best?
#================================#

# in terms of squared prediction error?
# CV 
  # also allows us to choose bandwidth parameter
  # that are optimal for a given data set.

# To carry out CV
  # 1st: fit the model using subset of the data
  # 2nd: predict at the data location where omitted
  # 3rd: compute discrepancy, squared error btw predicted and observed

# denote the sum of the discrepancies for a particular parameter theta (bandwidth)
# as LOOCV score. 


# simple way: 
  # compute pairwise distances btw all obs locations
  # and associated kernel-weight matrix
  # select the appropriate rwos and cols from resulting matrix
  # to do prediction at a left-out obs
  # repeat this for every obs

obs_obs_dist_mat <- rdist(select(Tmax, lon, lat, day), 
      select(Tmax, lon, lat, day))



LOOCV_score <- function(Wt_fun, theta, dist_mat, Z) {
  Wtilde <- Wt_fun(theta, dist_mat)
  
  CV <- 0
  for(i in 1:length(Z)) {
    Wtilde2 <- Wtilde[i, -i]
    W2 <- Wtilde2 / sum(Wtilde2)
    z_pred <- W2 %*% Z[-i]
    
    CV[i] <- (z_pred - Z[i])^2
  }
  
  mean(CV)
}


#-------------------------------#
# Two different kernel functions
#-------------------------------#

LOOCV_score(Wt_fun = Wt_IDW,
            theta = 5,
            dist_mat = obs_obs_dist_mat,
            Z = Tmax$z)

# [1] 7.775333

LOOCV_score(Wt_fun = Wt_Gauss,
            theta = 0.5, 
            dist_mat = obs_obs_dist_mat,
            Z = Tmax$z)

# [1] 7.526056

# Gaussian kernel smoother performs marginally better

# But how do we know the chosen kernel bandwidth are suitable?
  # an objective to set the bandwidth par is to 
  # put them equal to those values that minimize the LOOCV scores
  # can be done by a set of plausible bandiwidths
  # 

theta_IDW <- seq(4, 6, length = 21)
theta_Gauss <- seq(0.1, 2.1, length = 21)

CV_IDW <- CV_Gauss <- 0
for (i in seq_along(theta_IDW)){
  CV_IDW[i] <- LOOCV_score(Wt_fun = Wt_IDW,
              theta = theta_IDW[i],
              dist_mat = obs_obs_dist_mat,
              Z = Tmax$z)
  
  
  CV_Gauss[i] <- LOOCV_score(Wt_fun = Wt_Gauss,
              theta = theta_Gauss[i], 
              dist_mat = obs_obs_dist_mat,
              Z = Tmax$z)
  
}


#------#
# plot
#------#

par(mfrow = c(1, 2))
plot(theta_IDW, CV_IDW,
     xlab = expression(alpha),
     ylab = expression(CV[(m)](alpha)),
     ylim = c(7.4, 8.5), type = "o")
abline(h = min(CV_IDW), col = "red")
abline(v = theta_IDW[which.min(CV_IDW)], col = "red")

theta_IDW[which.min(CV_IDW)] # [1] 5
min(CV_IDW) # [1] 7.775333



plot(theta_Gauss, CV_Gauss,
     xlab = expression(theta),
     ylab = expression(CV[(M)](theta)),
     ylim = c(7.4, 8.5),
     type = "o")

abline(h = min(CV_Gauss), col = "green")
abline(v = theta_Gauss[which.min(CV_Gauss)], col = "green")

theta_Gauss[which.min(CV_Gauss)] # [1] 0.6
min(CV_Gauss) # [1] 7.468624




























