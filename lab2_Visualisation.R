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



































