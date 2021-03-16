#**********************************#
# Lab 2.3 Exploratory Data Analysis
#**********************************#

# focus on the construction and visulaization
# of empirical means and covariances
# the use of empirical othogonal functions 
# their associated prinicple component time series
# semi variogram analysis
# spatio-temporal canonical correaltion analysis

# 1st part of this lab, consider:
# daily maximum temp in the NOAA data set btw
# May 1993 to Sep 1993

install.packages("fields")
install.packages("CCA")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("gstat")
install.packages("sp")
install.packages("spacetime")


library(fields)
library(CCA)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gstat)
library(sp)
library(spacetime)
library(STRbook)

set.seed(01-02-2021)

#===================#
# Load the data set
#===================#

data("NOAA_df_1990", package = "STRbook")

# to keep the data size manageable, take a subset
Tmax <- filter(NOAA_df_1990, 
       proc == "Tmax" &
         month %in% 5:9 &
         year == 1993)


Tmax$t <- Tmax$julian - 728049
head(Tmax)

#~~~~~~~~~~~~~
length(unique(Tmax$lat)) # [1] 124
length(Tmax$lat) #  20306

head(Tmax$t)
length(unique(Tmax$t)) # 153 julian

t <- Tmax$t == 1
length(t[t != TRUE]) # 20173
length(t[t == TRUE]) # [1] 133 julian = 1

20306-20173 # [1] 133 t == 1

# order(Tmax$t)
#~~~~~~~~~~~~~~~~~



#=======================#
# Empirical Spatial Mean
#=======================#

# The empirical spatial mean is a spatial quantity
# that can be stored in a new df that contains
# spatial lications and the respective average max temp at each location

# group by long lat, then compute the average max temp at the each of the separate lon-lat coords

sp_av <- group_by(Tmax, lat, lon) %>% 
  summarise(mu_emp = mean(z))
# average over time for the same lat, lon group
sp_av # 133*3


count(Tmax, lat, lon)
# 133 differtn combinations of (lat,lon) groups
# 153 time points to average over within each coords group

#~~~~~~~~~~~~~~~~~~
#g <- group_by(Tmax, lat, lon) # 20306 * 12
# grouped by lat, lon

#print(g, n = 90, width = Inf)
#summarise(g, n = n())
#~~~~~~~~~~~~~~~~~~~



#------#
# plot
#------#

# plot the average max temp per station
# see how this varies according to the lon, lat

lat_means <- ggplot(sp_av) + 
  geom_point(aes(x = lat, y = mu_emp)) + 
  xlab("Latitude (deg)") + 
  ylab("Max Temp (degF)") + 
  theme_bw()


lon_means <- ggplot(sp_av) + 
  geom_point(aes(x = lon, y = mu_emp)) + 
  xlab("Longitutde (deg)") + 
  ylab("Max Temp (deg F)") + 
  theme_bw()


print(lat_means) # Tmax decreases with lat 
print(lon_means) # no pattern



#=========================#
# Empirical Temporal Means
#=========================#

Tmax_av <- group_by(Tmax, date) %>% 
  summarise(meanTmax = mean(z))

# Tmax_av a df containing the average max temp
# on each day (average across all the stations)


#------#
# plot
#------#

gTmaxav <- ggplot() +  # an empty framework
  geom_line(data = Tmax, # original data
            aes(x = date, y = z, group = id),
            colour = "blue", alpha = 0.04) +
  geom_line(data = Tmax_av, 
            aes(x = date, y = meanTmax)) + 
  xlab("Month") +
  ylab("Max temp (deg F)") + 
  theme_bw()

print(gTmaxav)


head(Tmax)
#======================#
# Empirical Covariances
#======================#

# before obtaining the empirical covariances
# important all trends are removed not just intercept

# one simply way: 
  # 1st fit a linear model (has spatial and/or temporal covariates) to the data
  # then plot the empirical covariances of the detrended data (i.e. residuals)


# linear model fitting use lm function,
# the residuals from lm can be incorporated into the original df Tmax

# as from time series plot, observe quadratic tendency of temp
# over the chosen time span

lm1 <- lm(z ~ lat + t + I(t ^ 2), data = Tmax)
# I() enables the power sing ^ to be an arithmetic operator

Tmax$residual <- residuals(lm1)


#-------------------------------------------#
# get the spatial locations of statioins
#-------------------------------------------#

spat_df <- filter(Tmax, t == 1) %>%  # chose 1 day snapshot, has 133 stations (different groups of (lat, lon)) 
  dplyr::select(lon, lat) %>%
  arrange(lon, lat) # ascending by lon/lat

head(spat_df)
str(spat_df) # data.frame':	133 obs. of  2 variables:

#----------
head(Tmax, 2)
t <- Tmax$t == 1
length(t[t != TRUE]) # 20173
20306-20173 # [1] 133

head(f, 2)
s <- dplyr::select(f, lon, lat)
A <- arrange(s, lon, lat)
head(A); tail(A) # ascend by lon
#----------


# number of stations (different groups of (lat, lon))

m <- nrow(sp_av) #  133


#============================#
# Empirical covariance matrix
#============================#

# use cov in R

# when few missing data, use = "complete.obs"

# when lots of missing data, need imputation or only a subset of stations

# to compute the empirical covaraince matrices, 
# 1st need to put the data into space-wide format using spread

X <- dplyr::select(Tmax, lon, lat, residual, t) %>%
  spread(key = t, value = residual) %>%
  dplyr::select(-lon, -lat) %>%
  t() # each row is t, each col is spatial residuals

dim(X) # [1] 153 time points 133 residuals at each locations


#~~~~~~~~~~
s <- dplyr::select(Tmax, lon, lat, residual, t)
head(s, 2)
sprd <- spread(s, key = t, value = residual )
head(sprd, 2)
#~~~~~~~~~~~


#-----------------#
# lag-0 covariance
#-----------------#

Lag0_cov <- cov(X, use = "complete.obs")

View(cov)
# "all.obs", "complete.obs", "pairwise.complete.obs"

dim(Lag0_cov) # 133 133

range(Lag0_cov) # [1] -13.85848  75.05377

#-----------------#
# lag-1 covariance
#-----------------#

# lat-1 cov is the cov between the residuals from X[-1, ]
# and X[-nrow(X), ]

lag1_cov <- cov(X[-1, ], X[-nrow(X), ], use = "complete.obs")
str(lag1_cov)

#----------------------#
# make sense of the cov
#----------------------#

# two dim space do not have any specific ordering
# order the station by long and then plot the permuted sp cov matrix

# works best when domain of interest is rectangular with long span much larger than lat

# here, a square domain, so to split the domain into lat strip
# then plot the spatial cov matrix associated with each strip

# here split the long into 4 strips
head(spat_df, 2)
spat_df$n <- 1:nrow(spat_df) # assign an index to each station
lim_lon <- range(spat_df$lon)
lon_strips <- seq(lim_lon[1], lim_lon[2], length.out = 5) # 4 long strips
spat_df$lon_strip <- cut(spat_df$lon, lon_strips, 
                         labels = FALSE, include.lowest = TRUE)

head(spat_df)
tail(spat_df)

cut

# now we know in which strip each station falls
# we subset the station data frame by strip then 
# sort the subsetted df by latitude

View(plot_cov_strips)



#================#
# plot_cov_strips
#================#

#---------------------#
# self-write function
#---------------------#

emp_cov_plt <- function(C, spat_df) {
  require(fields)
  
  for (i in seq_along(unique(spat_df$lon_strip))) {
    spat_strip <- filter(spat_df, lon_strip == i) %>% 
      arrange(lat)
    
    idx <- spat_strip$n  # the station idx within this strip
    jitter <- seq(0, 1e-4, length = length(idx))
    
    image.plot(x = spat_strip$lat + jitter, 
               y = spat_strip$lat + jitter,
               z = C[idx, idx], 
               xlab = "Latitude", ylab = "Latitude",
               zlim = c(-15, 85), 
               col = tim.colors(10), cex = 200)
  }
}

plotcov(Lag0_cov)

#--------------#
# View function
#--------------#

function (C, spat_df) 
{
  require(fields)
  for (i in seq_along(unique(spat_df$lon_strip))) {
    spat_strip <- spat_df %>% filter(lon_strip == i) %>% 
      arrange(lat)
    idx <- spat_strip$n
    jitter <- seq(0, 1e-04, length = length(idx))
    image.plot(spat_strip$lat + jitter, spat_strip$lat + 
                 jitter, C[idx, idx], xlab = "latitude", ylab = "latitude", 
               zlim = c(-15, 85), col = tim.colors(10), cex = 200)
  }
}



#------#
# plot
#------#

plot_cov_strips(Lag0_cov, spat_df = spat_df)

emp_cov_plt(Lag0_cov, spat_df = spat_df)

#------------#
# Understand
#------------#

# the empirical spatial covariance matrices reveal the presence 
# of spatial correlation in the residuals

# but 4 lag-0 plots seems to be qualitatively similar
# indicates no stronge depence on longitude

# but there's a dependence on latitude
# spatial covariance decreases with decreasing latitude
# so the covariance change as latitude changes
# such dependence is a type of spatial non-stationary
# such plots can be used to assess non-stationary 
# and the requireness of sp-temp models



x_try <- seq(0, 600, length.out = 10000)

cov


















