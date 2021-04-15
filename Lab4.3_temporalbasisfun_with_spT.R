#******************************************#
# Lab 4.3 Temporal Basis Functions with SpT
#******************************************#

# model maximum temp in NOAA dataset using 
# temporal basis functions and sp rd fields

# Y(s;t) = x(s;t)'beta + Sum phi_i(t)a_i(s) + nu(s;t)
  # x: covariates
  # beta: regression coefficients
  # phi_i(t): temporal basis functions
  # a_i(s) : coeffients of the temporal basis funs, 
             # and modeled as mulitvariate spatial random fields
  # nu(s;t): spatially correlated, but temp independent rd process
# sp-t modelling using temporal basis function using SpatioTemporal 

library("dplyr")
library("ggplot2")
library("gstat")
library("RColorBrewer")
library("sp")
library("spacetime")
library("STRbook")
library("tidyr")

install.packages("SpatioTemporal")
library("SpatioTemporal")

# would need a dataframe with col: 
  # ID station as char
  # obs: data
  # date: date
  # can be created using transmute from dplyr

data("NOAA_df_1990", package = "STRbook")
NOAA_sub <- filter(NOAA_df_1990, 
       year == 1993 &
         month == 7 &
         proc == "Tmax")

NOAA_sub_for_STdata <- NOAA_sub %>% 
  transmute(ID = as.character(id),
            obs = z,
            date = date)
# transmute adds new variables/ delete existing one/ preserve existing ones

head(NOAA_sub) # 4122 * 10
head(NOAA_sub_for_STdata) # 4122 * 3

# as the covariates will be modeled the spatially varying effects
# also need to be supplied as df
  # only station coordinates will be considered as covariates
  # stations are extracted from the Tmax

covars <- select(NOAA_sub, id, lat, lon) %>%
  unique() %>%
  rename(ID = id)


# since the model assumes nu(s;t) is temporally uncorrelated
# all temporal variability needs to be captured via covariates
# or the basis function

  # so need to check if data exhibit temporal autocorelation before
  # adding any temporal basis functions


## the role of temporal basis function is to adequately 
# capture temporal modes of variation

# by adding temporal basis functions, lag-1 autocorrelation coeffiect
# is no longer significant

# one should add basis until temporal autocorrelation in the data
# at most stations is considerablely reduced. 

  # e.g. when n.basis = 2 two temporal basis fun for capturing temp variation
  # the prop of stations with Residuals exhibiting a sig lag-2 autocorrelation coeffi is 26% 
  # much lower than w/o basis function to explain the temp variablity 69%

# the basis function will be available in STdata$trend


## the spatial quantities {a_i(s)} are themselves modeled as spatial fields
  # once the {phi_i(t)} are decleared, empirical estimates of {a_i(s)} can be found


all_dates <- NOAA_sub$date %>% unique()
lookup <- data.frame(date = all_dates, 
           V1 = scale(as.numeric(all_dates)))

par(mai = c(rep(0.8, 4)))
plot(lookup)


# want to create a function that takes a Date obj as input
# returns the required covariate values

# function returns the covariates in a data frame
# at required dates
fnc <- function(dates) {
  left_join(data.frame(data = dates), lookup, by = "date") %>%
    select(-date)
}














































