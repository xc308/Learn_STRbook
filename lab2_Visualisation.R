#***********************#
# Lab 2.2: Visualization
#***********************#

# visulaize maximum temp in NOAA data set
# max record btw May 1993 and Sep 1993

install.packages("animation")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("gstat")
install.packages("maps")


library(animation)
library(dplyr)
library(ggplot2)
library(gstat)
library(maps)
library("STRbook")


set.seed(31-01-2021)

data("NOAA_df_1990", package = "STRbook")

Tmax <- filter(NOAA_df_1990, proc == "Tmax" &
         month %in% 5:9 &
         year == 1993)

head(Tmax, 2)

Tmax %>% select(lon, lat, date, julian, z) %>% head()


# the 1st record of julian is 728050, 
# correspond to the 1 may 1993
# to ease the following operations,
# create a new time varaible t = 1 when julian == 728050

Tmax$t <- Tmax$julian - 728049


#================#
# Spatial Plots
#================#

#---------------------#
# Point-reference data
#---------------------#

# NOAA data are collected at stations that are fixed in space
# initial plots should give an overall spatial variation of the obs
# if there are many time points, only choose a selection of time points for visulaization
# 

Tmax_1 <- subset(Tmax, t %in% c(1, 15, 30))
# choose the 1st， 15th, 30th day in Tmax
head(Tmax_1)

NOAA_plot <- ggplot(Tmax_1) +               # plot points
  geom_point(aes(x = lon, y = lat, colour = z),
             size = 2) +       # enlarge points
  col_scale(name = "degF") +   # attach color scale
  xlab("Longitutde (deg)") + 
  ylab("Lattitude (deg)") +
  geom_path(data = map_data("state"),
            aes(x = long, y = lat, group = group)) +  # use group to identify which points corresponds to which county
  facet_grid(~ date) +
  coord_fixed(xlim = c(-105, -75),
              ylim = c(25, 50)) +    # zoom in
  theme_bw() +
  #coord_map(projection = "sinusoidal")
  coord_map(projection = "ortho")
 
print(NOAA_plot)


# also a good practice to put the spatial location
# of the data into perspective by plotting the boundaries 
# together with the data locations

# above, the US state boundaries are obtained from 
# the maps package through the command map_data("state")
# the boundaries are then overlayed on the plot using geom_path

# projections can be applied by adding another layer to the ggplot2 obj
# using coord_map
# +coord_map(projection = "sinusodal")

# or plot in 3-d projection = "ortho"


#-----------------------------#
# Plot of regular lattice data
#-----------------------------#

geom_tile()


#-----------------------------#
# Plot of irregular lattice data
#-----------------------------#

geom_polygon()



#===================#
# time-series plots
#===================#

# it's recommanded to plot the time seriess
# at all 139 staions
# but here just plot the time series at a set of stations
# selected at random

# 1st obtain a set of unique station identifiers
# choose 10 at random
# and then extract the data associated with these 10 stations from the data set

UIDs <- unique(Tmax$id) # 133
UIDs_sub <- sample(UIDs, 10)

filter # subset a df, retaining all rows satisfy all the conditions
Tmax_sub <- filter(Tmax, id %in% UIDs_sub)


# to visualize the time series at these sub stations
# use facets

# when given a long data frame
# one can first subdivide the df into groups
# and generate a plot for each group


# now plot the time sereise for each selected station
head(Tmax_sub, 2)

TmaxTS <- ggplot(Tmax_sub) + 
  geom_line(aes(x = t, y = z)) +
  facet_wrap(~ id, ncol = 5) + # facet by station
  xlab("Number of days (day)") + 
  ylab("Tmax (degF)") + 
  theme_bw() + 
  theme(panel.spacing = unit(1, "lines")) # facet spacing

# ~id is a formula used to denote the groups we are faceting
# x ~ y can be used to facet by two varaibles

print(TmaxTS)



#================#
# Hovmoller Plots
#================#

# a two-dim space time visualization, 
# where space is collapsed (projected or averaged)
# onto one dim; the 2nd dim is time

# the Hovmoller plot can be generated easily
# if the data are on space-time grid

# Consider a latitudinal H- plot
# 1st, generate a regular grid of 25 spatial points
# and 100 temporal points using expand.grid function
# with limit set to both the latitude and t limit

lim_lat <- range(Tmax$lat)
lim_t <- Tmax$t

lat_axis <- seq(lim_lat[1], lim_lat[2], length.out = 25)
t_axis <- seq(lim_t[1], lim_t[2], length.out = 100)

lat_t_grid <- expand.grid(lat_axis, t_axis)


# Next, associate each station's latitudal coords with the closest one on the grid
# by finding the distance from the station's latitudinal coords to each point of the grid
# finding which gridpoint is the closest, allocation that to it

# store the grided data in Tmax_grids

Tmax_grid <- Tmax
o <- outer(Tmax$lat, lat_axis, "-") 
# distance from the station's latitudinal coords to each point of grid
# 0 : 20306 by 25

dists <- abs(outer(Tmax$lat, lat_axis, "-"))
dim(dists) # [1] 20306    25

# find the closest grid point for each of the station
m <- apply(dists, 1, which.min) # a vector 20306
all(head(m, 153) == 14) # [1] TRUE
# mean the first 153 stations are all closest to the 14th grid
tail(m)
range(m) # [1]  1 25

min_g <- lat_axis[apply(dists, 1, which.min)]
lat_axis[14] # [1] 39.57222

Tmax_grid$lat <- lat_axis[apply(dists, 1, which.min)]
# the latitude of the grid point closest to each station 

# now have associated each station with closest grid point
# left is to group by latitude and time

Tmax_lat_Hov <- group_by(Tmax_grid, lat, t) %>% summarise(z = mean(z))
# group by lat, t, 
# summarize group by lat, average over time within each lat group

# Note, in this case, every lat-t band consists at least one data point
# so the H plot contains no missing points on the established grid

# if this is not the case, interp from akima package can be used to fill out the
# grid cell with no data


#------#
# plot 
#------#

Hovmoller_lat <- ggplot(Tmax_lat_Hov) +     # take the data
  geom_tile(aes(x = lat, y = t, fill = z)) +  # plot
  fill_scale(name = "deg F") +        # add the legend and it's title
  scale_y_reverse() +        # reverse time scale
  ylab("Days")      +
  xlab("Lattitude (degrees)")  + 
  theme_bw()


print(Hovmoller_lat)







































