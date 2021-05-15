#********************#
# Analyzing Residuals
#********************#

# if there're still sp-t correlation in the residual
# then our model will not have caputured adequately 
# the sp-t variablilty in the data. 

install.packages("fields")
library(fields)

install.packages("ape")
library(ape) # Moran's I test



Tmax_no_14$residuals <- residuals(Tmax_July_lm)


#===================#
# Plot the residuals
#===================#

# choose 8 days
g <- ggplot(filter(Tmax_no_14, day %in% 24:31)) + 
  geom_point(aes(lon, lat, colour = residuals)) + 
  facet_wrap(~ day, ncol = 4) + 
  col_scale(name = "defF") + 
  geom_point(data = filter(Tmax_no_14, day %in% 24:31 &
                             id %in% c(3810, 3889)),
             aes(lon, lat), 
             colour = "black",
             pch = 2, size = 2.5) + 
  theme_bw()

print(g)

# the residual time series exhibit considerable temporal correlation
  # the closer the residual in time, the more similar residuals are

# the spatial residual exhibit spatial 



#===============#
# Moran's I test
#===============#

#----------#
# day 1 ref
#----------#

Tmax_day1 <- filter(Tmax_no_14, day == 1)
station1.dist <- filter(Tmax_no_14, day == 1) %>%
  select(lon, lat) %>%
  dist() %>%
  as.matrix()

str(station1.dist) # 133: 133

stations.dists.inv <- 1/station1.dist
diag(stations.dists.inv) <- 0

M1 <- Moran.I(Tmax_day1$residuals, weight = stations.dists.inv)

str(M1)
# List of 4
#$ observed: num 0.272
#$ expected: num -0.00758
#$ sd      : num 0.0124
#$ p.value : num 0

M1 %>% do.call("cbind", .)
#       observed     expected         sd p.value
# [1,] 0.2716679 -0.007575758 0.01235583       0


#-----------------#
# day 1:13, 15:31
#-----------------#

str(Tmax_no_14)
head(Tmax_no_14)
3989/ 30

P <- list()
days <- c(1:13, 15:31)
for (i in seq_along(days)) {
  Tmax_day <- filter(Tmax_no_14, day == days[i])
  
  station.dists <- Tmax_day %>% select(lon, lat) %>% dist() %>% as.matrix()
  station.dists.inv <- 1/station.dists # weight matrix
  diag(station.dists.inv) <- 0
  
  P[[i]] <- Moran.I(Tmax_day$residuals, station.dists.inv) %>%
    do.call("cbind", .)
  
}

str(P) # list of 30 rows

do.call("rbind", P) %>% head() # df

P_df <- do.call("rbind", P)

range(P_df[, 4]) # [1] 0.000000e+00 8.040186e-06
# very small, rj Null hypo, strong sp dependence


#===================#
# Durbin-Watson Test
#===================#

# when data are regularly spaced in time
# look at the temporal residuals at some location
# and test for temporal correlation in these residuals 
# using Durbin-Watson test

# 1st filter out two locations to see time series

TS1 <- filter(Tmax_no_14, id == 3810)$residuals
TS2 <- filter(Tmax_no_14, id == 3889)$residuals

str(TS1)
range(TS1) # [1] -6.672968  6.387341


#---------------#
# plot residuals
#---------------#

par(mar = c(4, 4, 1, 1))

plot(TS1, xlab = "day of July 1993",
     ylab = "residuals (degF)",
     type = "o", ylim = c(-8, 7))

lines(TS2, xlab = "day of July 1993",
     ylab = "residuals (degF)",
     type = "o", col = "red")

# residuals close to each other in time tend to 
# be more similar than residuals further apart

# also note that at the same day, two close staions
# residuals are correlated (similar values)
# due to sp correlation in the residuals

# also look at the empirical autocorrelation function of residuals
acf(TS1)
acf(TS2)
# both lag-1 correlation in the residuals


#-------------------#
# Durbin-Watson Test
#-------------------#

# for residuals at every stations
# use tidyr, purrr, broom for tests and predications on groups of data

# group data: group_by
# df containing one row for each group: nest()
  # so each group is a df
# perform opreration on each of group: summarise()

# group Tmax_no_14 by lon, lat
nested_Tmax_no_14 <- group_by(Tmax_no_14, lon, lat) %>% nest() 
# each (lon, lat) group has a df of 30 days obs

head(nested_Tmax_no_14, 3)
str(nested_Tmax_no_14)
#     lat   lon data              
#   <dbl> <dbl> <list>            
#  1  39.3 -81.4 <tibble [30 × 16]>
#  2  35.7 -81.4 <tibble [30 × 16]>
#  3  35.6 -88.9 <tibble [30 × 16]>

# Note: data in the 3rd col is a list


# Now define a function that 
  # takes the data frame of a single group,
  # carries out the test dwtest
  # returns the results

# here we test autocorrelation in the residuals
# after removing a temporal constant trend by residuals ~ 1

dwtest_one_station <- function(data) {
  dwtest(residuals ~ 1, data = data)
}

dwtest_one_station(nested_Tmax_no_14$data[[1]]) # 1st (lon, lat)


## carry out the test at each station in the nested df, 
# use map() in purrr

library(purrr)
map(nested_Tmax_no_14$data, dwtest_one_station) %>% head()


## assign results to another col in the nested df: mutate()
## these results are of class htest, uneasy to analyze: tidy() from broom
  # to extract key info: statistic, p-value, method, hypo

dwtest_one_station_tidy <- nested_Tmax_no_14$data[[1]] %>%
  dwtest_one_station() %>%
  tidy()

dwtest_one_station_tidy 
# # A tibble: 1 x 4
#statistic p.value method  alternative  
#<dbl>   <dbl> <chr>   <chr>        
#  1     0.982 0.00122 Durbin… true autocor…


## To assign the test result to each station in nested df as added fields
  # using: unnest()

unnest(dwtest_one_station_tidy)


Tmax_DW_no_14 <- nested_Tmax_no_14 %>%
  mutate(dwtest = map(data, dwtest_one_station)) %>%
  mutate(test_df = map(dwtest, tidy)) %>%
  unnest(test_df)


head(Tmax_DW_no_14) 
# # A tibble: 6 x 8
# Groups:   lat, lon [6]
#lat   lon data  dwtest statistic p.value method
#<dbl> <dbl> <lis> <list>     <dbl>   <dbl> <chr> 

Tmax_DW_no_14 %>% select(-method, -alternative) %>% head(3)


## the proportion of p-values below (5% divided by the # of test)
mean(Tmax_DW_no_14$p.value < (0.05 / nrow(Tmax_DW_no_14))) * 100
# 21.8% of stations are tested with autocorrelation > 0 in the residuals







































































