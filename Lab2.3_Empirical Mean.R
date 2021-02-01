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

length(unique(Tmax$lat)) # [1] 124
length(Tmax$lat) #  20306


#=======================#
# Empirical Spatial Mean
#=======================#

# The empirical spatial mean is a spatial quantity
# that can be stored in a new df that contains
# spatial lications and the respective average max temp at each location

# group by long lat, then compute the average max temp at the each of the separate lon-lat coords

sp_av <- group_by(Tmax, lat, lon) %>% summarise(mu_emp = mean(z))
# average over time for the same lat, lon group

g <- group_by(Tmax, lat, lon) # 20306
# grouped by lat, lon, 31 days


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



#======================#
# Empirical Covariances
#======================#








