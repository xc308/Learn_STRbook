install.packages("devtools")
install.packages("usethis")

library(devtools)
install_github("andrewzm/STRbook")

install.packages("dplyr")
install.packages("tidyr")

library(dplyr)
library(tidyr)
library(STRbook)

install.packages("spacetime")
library(spacetime)


#****************#
# Load sp-t data 
#****************#

# station info 328 * 3
locs <- read.table(system.file("extdata", "stationinfo.dat",
            package = "STRbook"),
            col.names = c("id", "lat", "lon"))


# daily time points 1461 * 4
Times <- read.table(system.file("extdata", "Times_1990.dat", 
            package = "STRbook"), 
           col.names = c("julian", "year", "month", "day"))


# max Temp 1461 * 328
Tmax_2 <- read.table(system.file("extdata", "Tmax_1990.dat", 
            package = "STRbook"))

head(Tmax_2)
names(Tmax) <- locs$id


# min temp 1461 * 328
Tmin <- read.table(system.file("extdata", "Tmin_1990.dat", 
            package = "STRbook"))
head(Tmin, 2)
names(Tmin) <- locs$id


# temp dew point
TDP <- read.table(system.file("extdata", "TDP_1990.dat", 
            package = "STRbook"))

names(TDP) <- locs$id
head(TDP, 2)

# preciption
Precip <- read.table(system.file("extdata", "Precip_1990.dat",
            package = "STRbook"))



#****************************#
# Transform into long format
#****************************#

# attach Tmax to Times since each row is an obs

Tmax <- cbind(Times, Tmax)
head(names(Tmax), 10)
# space wide
# need to transform into a long format
# key-value pair, keys are station IDs
# values are the MaxTemp stored in z

head(Tmax, 10)
str(Tmax)

Tmax_long <- gather(Tmax, id, z, -julian, -year, -month, -day)

head(Tmax_long)
str(Tmax_long)
Tmax_long$id <- as.integer(Tmax_long$id)

# filter out missing data 
# better to use inequality (e.g. less than)
# rather than equality as criterion

nrow(Tmax_long) # [1] 479208

Tmax_long <- filter(Tmax_long, !(z <= -9998)) # filter out unwanted
nrow(Tmax_long) # [1] 196253

# wish to include minimum temp and other 
# varaibles inside the df

# use mutate function
# add one column proc, indicating what process the measurement relates to

Tmax_long <- mutate(Tmax_long, proc = "Tmax")
head(Tmax_long, 3)


#==============#
# load the rest
#==============#
data("Tmin_long", package = "STRbook")
data("TDP_long", package = "STRbook")
data("Precip_long", package = "STRbook")

# now concatenating rowwise of all data 
NOAA_df_1990 <- rbind(Tmax_long, Tmin_long, TDP_long, Precip_long)


#========================#
# Advantages of long from
#========================#

# easy to make grouping and summarizing 
summ <- group_by(NOAA_df_1990, year, proc) %>%
  summarize(mean_proc = mean(z))

str(summ)

head(summ)


# want to find out the number of days when
# it didnot rain at each station in Jun of every year

NOAA_precip <- filter(NOAA_df_1990, proc == "Precip" & month == 6)
summ_no_precip <- group_by(NOAA_precip, year, id) %>% 
  summarise(days_no_precip = sum(z == 0))

head(summ_no_precip)

range(summ_no_precip$days_no_precip) # [1] 12 28
median(summ_no_precip$days_no_precip) # [1] 20


#===================#
# commands in dplyr
#===================#

#-----------------------#
# arrange sorts by a col
#-----------------------#

# want the NOAA_df_1990 sorted 1st by time and then by ID
NOAA_df_sorted <- arrange(NOAA_df_1990, julian, id)
head(NOAA_df_sorted)
head(NOAA_df_1990)


#-------#
# select 
#-------#

# used to selcet or discard cols
df1 <- select(NOAA_df_1990, julian, z)
df2 <- select(NOAA_df_1990, -julian)


#-----------#
# left_join
#-----------#

# left_join is considerablity faster than merge
# need to supply the col fiedl name common to both df
NOAA_df_1990 <- left_join(NOAA_df_1990, locs, by = "id")


#--------#
# spread
#--------#

# spead is the revert function of gather
# also works by identifying the key-value pair in the df
# the values are widened into a  table
# keys are used to label the col

dim(Tmax_long) # [1] 196253  7
Tmax_long_sel <- select(Tmax_long, julian, id, z)


Tmax_wide <- spread(Tmax_long_sel, key = id, value = z)
dim(Tmax_wide) # [1] 1461  138
head(Tmax_wide) # the 1st column is julian

# to construct a standard matrix containing id and z only
M <- select(Tmax_wide, -julian) %>% as.matrix()
dim(M) # [1] 1461  137






















