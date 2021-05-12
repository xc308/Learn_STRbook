#***************#
# Learn sp_t cv
#***************#

# ref: https://cran.r-project.org/web/packages/CAST/CAST.pdf

install.packages("CAST")
library(CAST)

View(CreateSpacetimeFolds)

View(createFolds)


install.packages("spatialrisk")
library(spatialrisk)
circle <- points_in_circle(Groningen, lon_center = 6.571561, lat_center = 53.21326, radius = 100)
circle

points_in_circle()





