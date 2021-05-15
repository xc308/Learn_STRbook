#*****************************************#
# Lab 3.3 Regression Models for Forecasting
#*****************************************#

# fix a simple linear regression to every pixel in SST daat
# use these model to predict SST for a month in which 
# we have no SST data

library(broom)
library(dplyr)
library(ggplot2)
library(STRbook)
library(purrr)
library(tidyr)

data("SST_df", package = "STRbook")


#=================#
# Tidy Up the data
#=================#

data("SSTlandmask", package = "STRbook")
data("SSTdata", package = "STRbook")
data("SSTlonlat", package = "STRbook")

## combine land mask data wht coords df
lonlatmask_df <- data.frame(cbind(SSTlonlat, SSTlandmask))
names(lonlatmask_df) <- c("lon", "lat", "mask")

SSTdata <- cbind(lonlatmask_df, SSTdata)
str(SSTdata) # 'data.frame':	2520 obs. of  402 variables:

# long format
SST_df <- gather(SSTdata, date, sst, -lon, -lat, -mask)
str(SST_df)

length(unique(SST_df$date))

# to rename the date label from V1, V2, to 
# Jan 1970, Feb 1970,...using a mapping table

date_grid <- expand.grid(Month = c("Jan", "Feb", "Mar", 
                      "Apr", "May", "Jun",
                      "Jul", "Aug", "Sep",
                      "Oct", "Nov", "Dec"),
            Year = 1970:2002,
            stringsAsFactors = F)
str(date_grid)


date_grid$date <- paste0("V", 396)
SST_df <- left_join(SST_df, date_grid) %>%
  select(-date)

str(SST_df)
head(SST_df)

## set SST data that are coincident with land locations to NA


data("SOI", package = "STRbook")
data("SST_df", package = "STRbook")

str(SOI)
SOI_df <- SOI %>% select(-Ann) %>%
  gather(Month, soi, -Year)

str(SOI_df)
# 'data.frame':	1824 obs. of  3 variables:

str(SST_df)


SST_df <- left_join(SST_df, SOI_df, by = c("Month", "Year"))
str(SST_df ) # 'data.frame':	1005480 obs. of  9 variables:



#=========================#
# Fitting models pixelwise
#=========================#

# fit linear time-series models to SSTs in each pixel
# using data from Jan 1970 - Apr 1997

SST_pre_May <- filter(SST_df, Year <= 1970) %>%
  filter(!(Year == 1997 & Month %in% 
             c("May", "Jun", "Jul", "Aug",
               "Sep", "Oct", "Nov", "Dec")))


## linear model fit data over time to each pixel
fit_one_pixel <- function(data) {
  mod <- lm(sst ~ 1 + soi.x, data = data)
}

res_one_pixel <- function(mod) {
  resd <- residuals(mod)
}

pixel_lms <- SST_pre_May %>%
  filter(!is.na(sst)) %>%
  group_by(lon,lat) %>%
  nest() %>%
  mutate(model = map(data, fit_one_pixel)) %>%
  mutate(resd = map(model, res_one_pixel)) %>%
  mutate(model_df = map(model, tidy))

pixel_lms %>% head(3)
# # A tibble: 3 x 5
# Groups:   lon, lat [3]
# lon   lat   data              model     model_df        
#<dbl> <dbl> <list>             <list>     <list>          
#  1   154   -29 <tibble [12 × 7]> <lm>   <tibble [2 × 5]>
#  2   156   -29 <tibble [12 × 7]> <lm>   <tibble [2 × 5]>


lm_parmeters <- pixel_lms %>% unnest(model_df)
lm_parmeters %>% head(3)
# for each pixel, we have an estimate of intercept and slop for covariate soi
## A tibble: 3 x 9
# Groups:   lon, lat [2]
#lon   lat data         model  term     estimate std.error statistic p.value
#<dbl> <dbl> <list>       <list> <chr>       <dbl>     <dbl>     <dbl>   <dbl>
#  1   154   -29 <tibble [12… <lm>   (Interc…   0.0639    0.0819     0.780 0.454  
#  2   154   -29 <tibble [12… <lm>   soi.x     -0.0327    0.0772    -0.424 0.681  
#  3   156   -29 <tibble [12… <lm>   (Interc…   0.275     0.0773     3.56  0.00516


#-----------#
# Residuals
#-----------#

2261 * 12
head(pixel_lms$resd)
length(pixel_lms$resd) # 2261



R <- pixel_lms %>%
  select(lon, lat, resd) %>%
  unnest(resd)

head(R, 3)


lm_resd_mean <- pixel_lms %>% mutate(mean_resd = map(resd, mean)) 
head(lm_resd_mean, 2)

# # A tibble: 2 x 7
# Groups:   lon, lat [2]
# lon   lat data    model  resd  model_df  mean_resd
#<dbl> <dbl> <list>  <list> <lis> <list>    <list>   
#  1   154   -29 <tibbl… <lm>   <dbl… <tibble … <dbl [1]>
#  2   156   -29 <tibbl… <lm>   <dbl… <tibble … <dbl [1]>

#lm_resd_mean$mean_resid[1]
#map(pixel_lms$resd[1], mean)


lm_resd_mean$mean_resid <- unlist(lm_resd_mean$mean_resd)
length(unlist(lm_resd_mean$mean_resd)) #  2261

head(lm_resd_mean)
# # A tibble: 6 x 8
# Groups:   lon, lat [6]
#lon   lat data          model  resd     model_df     mean_resd mean_resid
#<dbl> <dbl> <list>        <list> <list>   <list>       <list>         <dbl>
#  1   154   -29 <tibble [12 … <lm>   <dbl [1… <tibble [2 … <dbl [1]>   5.78e-18
#  2   156   -29 <tibble [12 … <lm>   <dbl [1… <tibble [2 … <dbl [1]>   1.56e-17




#================#
# Plot out coeffs
#================#

# merge this df with coords df using left_join
str(lonlatmask_df)
left_join(lonlatmask_df, lm_parmeters)
# the regression coeff over land are marked as NA
lm_pars <- left_join(lonlatmask_df, lm_parmeters)
str(lm_pars)



g2 <- ggplot(filter(lm_pars, term == "(intercept)" | mask == 1)) + 
  geom_tile(aes(lon, lat, fill = estimate)) + 
  fill_scale() + 
  theme_bw() + coord_fixed()

g3 <- ggplot(filter(lm_pars, term == "soi.x" | mask == 1)) +
  geom_tile(aes(lon, lat, fill = estimate)) + 
  fill_scale() + 
  theme_bw() + coord_fixed()


plot(g2)
plot(g3)


ggplot(lm_resd_mean) +
  geom_tile(aes(lon, lat, fill = mean_resid)) +
  fill_scale() + 
  theme_bw() +
  coord_fixed()


g4 <- ggplot(lm_resd_mean) + 
  geom_tile(aes(lon, lat, fill = mean_resd) + 
  fill_scale() + 
  theme_bw() + coord_fixed()

