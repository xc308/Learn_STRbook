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
