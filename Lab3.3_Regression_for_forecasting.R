#******************************************#
# Lab 3.3 Regression Models for Forecasting
#******************************************#

#broom and purrr for fintting and predicting 
# with multiple models simulatneously
install.packages("broom")

library("broom")
library("dplyr")
library("ggplot2")
library("STRbook")
library("purrr")
install.packages("tidyr")
library("tidyr")
remove.packages("broom")


data("SST_df", package = "STRbook")
head(SST_df, 2)
# 1st: wrangle the SST data into a long-format df
# that's amenable to linear fiting and plotting
# SST are provided in 3 dfs, 
  # 1: describing the land mask
  # 2: containing SST value in wide format
  # 3: containing the coordinates

# 1st combine land mask data with the coords df
lonlatmask_df <- data.frame(cbind(SSTlonlat, SSTlandmask))
head(lonlatmask_df)
names(lonlatmask_df) <- c("lon", "lat", "mask")

# then form SST df in wide format by attaching SSTdata 
# to the coordinate mask data frame

head(SSTdata, 2)

SSTdata <- cbind(lonlatmask_df, SSTdata)
head(SSTdata, 3)

# put data into long format
SST_df <- gather(SSTdata, key = date, value = sst, -lon, -lat, -mask)

SST_df %>% head(3)
#   lon lat mask date        sst
# 1 124 -29    1   V1 -0.3628883
# 2 126 -29    1   V1 -0.2846108


## replace the date field wit two fields
# by mapping table that linds v1 to Jan 1970
# v2 to Feb 1970, 
# then merge using left_join

date_grid <- expand.grid(Month = c("Jan", "Feb", "Mar",
                      "Apr", "May", "Jun", 
                      "Jul", "Aug", "Sep",
                      "Oct", "Nov", "Dec"), 
                      Year = 1970:2002,
                      stringsAsFactors = FALSE)

date_grid$date <- paste0("V", 1:396) # 12 * 33 = 396
str(SST_df)
SST_df <- left_join(SST_df, date_grid) %>% 
  select(-date)


# for good measure, also add the date field
# in month-year format

SST_df$date <- paste(SST_df$Month, SST_df$Year)
SST_df %>% head(3)


# NOW the datasets is ready
g <- ggplot(filter(SST_df, Year == 1997 &
                Month %in% c("Apr", "Aug", "Jun", "Oct"))) + 
  geom_tile(aes(lon, lat, fill = pmin(sst, 4))) + 
  facet_wrap(~ date, dir = "v") + 
  fill_scale(limits = c(-4, 4),
             name = "degC") + 
  theme_bw() + coord_fixed()

print(g)

# Also need to add the SOI data to the SST df
# SOI: 
  # 1st col: year
  # 2nd - 13 col: SOI for each month within a year
  # last col: the mean SOI for that year

data("SOI", package = "STRbook")
head(SOI, 2)

SOI_df <- select(SOI, -Ann) %>%
  gather(key = Month, value = soi, -Year)

head(SOI_df)

SST_df <- left_join(SST_df, SOI_df,
          by = c("Month", "Year"))

head(SST_df)


#=============================#
# Fitting the Models Pixelwise
#=============================#

# fit linear time-series models to the SSTs
# in each pixedl using data up to Apr 1997

SST_pre_May <- filter(SST_df, Year <= 1997) %>%
  filter(Year == 1997 &
           Month %in% c("Jan", "Feb",
                        "Mar", "Apr"))

fit_one_pixel <- function(data) {
  mod <- lm(sst ~ 1 + soi, data = data)
}


pixel_lms <- SST_pre_May %>%
  filter(!is.na(sst)) %>%
  group_by(lon, lat) %>%
  nest() %>%
  mutate(model = map(data, fit_one_pixel)) %>%
  mutate(model_df = map(model, tidy))

# tidy turn an object into a tidy tible

pixel_lms %>% head(3)

# to extract the model parametr from the linar-fit
# df, use unnest

lm_pars <- pixel_lms %>%
  unnest(model_df)

lm_pars %>% head(5)


## Plot spatial maps of the intercept and regression coefficient
# associated with soi

# 1st merge this df with the coods df using left_join
lm_pars <- left_join(lonlatmask_df, lm_pars) # regression coefficient over land pixels in lonlat mask wil be NA


#------#
# Plot
#------#

g2 <- ggplot(filter(lm_pars, term == "(Intercept)" | mask == 1)) +
  geom_tile(aes(lon, lat, fill = estimate)) + 
  fill_scale() + 
  theme_bw() + coord_fixed()


g3 <- ggplot(filter(lm_pars, term == "soi" | mask == 1)) + 
  geom_tile(aes(lon, lat, fill = estimate)) + 
  fill_scale() + 
  theme_bw() + coord_fixed()


print(g2)
print(g3)



#==========================#
# Predicting SST Pixelwise
#==========================#

# use the linear models at the pixel level to predict
# SST in OCT 1997 using SOI index for that month

soi_pred <- filter(SOI_df, Month == "Oct" & Year == "1997") %>%
  select(soi)


# write a function takes lm and SOI at prediction date soi_pred
# run predict function for this date, returns a df containing the 
# prediction and prediction std error

predict_one_pixel <- function(lm, soi_pred) {
  predict(lm, 
          newdata = soi_pred,
          interval = "prediction") %>%
    data.frame() %>%
    mutate(se = (upr - lwr) / (2 * 1.96)) %>%
    select(fit, se)
}


# prediction proceeds at each pixel by calling
# predict_one_pixel on each row in nested df pixel_lms

SST_Oct_1997 <- pixel_lms %>%
  mutate(preds = map(model, predict_one_pixel,
                     soi_pred = soi_pred)) %>%
  unnest(preds)


SST_Oct_1997 %>% head(3)



































