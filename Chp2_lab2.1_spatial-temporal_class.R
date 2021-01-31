#********************************#
# working with Sp-Tp Data Classes
#********************************#

# convert data into obj of class STIDF, STFDF
# SP-temp locations which only need STI or STF ojb
# the obj of class STIDF, STFDF also contain data

# these classes are defined in package spacetime

# sometimes construct sp-temp objs using spatial objects
# defined in package sp

install.packages("sp")
install.packages("spacetime")
library(sp)
library(spacetime)


#***************************#
# Constructing an STIDF Obj
#***************************#

# the sp-temp obj for irregular data, STIDF
# can be constructed using two functions:
# stConstruct, STIDF

# the only thing need to do before call stConstruct
# is to define a formal time stamp from year, month, day field

# 1st construct a field with the date in year-month-day format
# using paste function along with "with" function


#===========================#
# Define a formal time stamp
#===========================#

NOAA_df_1990$date <- with(NOAA_df_1990, paste(year, month, day, sep = "-"))
head(NOAA_df_1990$date, 5)
# the date field is of type character. 
# now convert to Date Obj using as.Date

NOAA_df_1990$date <- as.Date(NOAA_df_1990$date)
class(NOAA_df_1990$date)
# [1] "Date"


#=====================================#
# Construct sp-temp obj of class STIDF
#=====================================#

#------------#
# stConstruct
#------------#

# provide the data frame in long format
# and indicate which are the spatial temporal coordinates

Tmax_long2 <- filter(NOAA_df_1990, proc == 'Tmax')

STObj <- stConstruct(Tmax_long2, 
            space = c("lon", "lat"),
            time = "date")

class(STObj) # [1] "STIDF"
# successfully generate an obj of class STIDF


#-------#
# STIDF
#-------#

# STIDF requires one to also specify the spatial
# part as an obj of class Spatial from the package sp

# in our case, the spatial component is an obj containing
# irregular spaced data, which is SpatialPoints obj in sp
# SpatialPoints function needs coords to be supplied

head(Tmax_long2, 3)
sp_part <- SpatialPoints(coords = Tmax_long2[, c("lon", "lat")])
temp_part <- Tmax_long2$date

STOjb2 <- STIDF(sp = sp_part, time = temp_part, 
      data = select(Tmax_long2, -date, -lon, -lat))

class(STOjb2) # [1] "STIDF"



#==========================#
# Constructing an STFDF Obj
#==========================#

# when spatial points are fixed in time
# only need to provide as many spatial coords as there are sp points
# i.e. coords for station locations

# also need to provide the regular time stamps
# one for each day btw 1Jan 1990 - 30 Dec 1993

# data can be space-wide or time-wide format


#-------------------------------------------#
# Get spatical and temporal part from origin
#-------------------------------------------#

head(locs, 2)
spat_part <- SpatialPoints(coords = locs[, c("lon", "lat")])
str(spat_part)

temp_part <- with(Times, paste(year, month, day, sep = "-"))
temp_part <- as.Date(temp_part)


#-------------------------------#
# gather the data in long format
#-------------------------------#

Tmax_long3 <- gather(Tmax, key = id, value = z, -julian, -year, -month, -day)

str(Tmax_long3)
# $ id    : chr  "3804"

Tmax_long3$id <- as.integer(Tmax_long3$id)


#------------------------------------------------------#
# make sure spatial index moving faster than time index
#------------------------------------------------------#

# so need to order time first then spatial index
Tmax_long3_ord <- arrange(Tmax_long3, julian, id)


#-------------------------------------------#
# make sure spatial index is the same order
# as that of the spatial component supplied
#-------------------------------------------#

all(unique(Tmax_long3_ord$id) == locs$id)
# [1] TRUE
# all()# are all the values of the supplied vec ture?


#----------------#
# construct STFDF
#----------------#

STObj3 <- STFDF(sp = spat_part,
      time = temp_part,
      data = Tmax_long3_ord)

class(STObj3)
# [1] "STFDF"


#---------#
# Add CRS
#---------#

proj4string(STObj3) <- CRS("+proj=longlat +ellps=WGS84")
proj4

STObj3$z[STObj3$z == -9999] <- NA

str(STObj3)






