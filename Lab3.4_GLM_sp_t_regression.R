#********************************************#
# Lab 3.4 Gerneralized Linear Sp-T Regression
#********************************************#

data("MOcarolinawren_long", package = "STRbook")
MOcarolinawren_long <- MOcarolinawren_long %>% filter(!is.na(cnt))

# use the same covariates to fit these data as we did to fit max temp
# 12 of these covariates were basis functions 
# constructed using auto_basis from FRK


install.packages("FRK")
library(FRK)

install.packages("sp")
library(sp)
G <- auto_basis(data = MOcarolinawren_long[, c("lon", "lat")] %>%
             SpatialPoints(),
           nres = 1,
           type = "Gaussian")


# evaluate the basis functions at the Carolina wren obs locations
S <- eval_basis(basis = G, 
           s = MOcarolinawren_long[, c("lon", "lat")] %>%
             as.matrix()) %>%
  as.matrix()
colnames(S) <- paste0("B", 1:ncol(S))


# attach the basis-function covariate information
# to the df containning the counts

head(MOcarolinawren_long, 3)
#   loc.ID t cnt      lat       lon year
# 1      1 1   4 36.82526 -89.24779 1994

Wren_df <- cbind(MOcarolinawren_long, S) %>%
  select(-loc.ID, -t)

Wren_df %>% head(3)


# fit 5 col of the 1st 3 records of constructed df Wren_df
Wren_df[1:3, 1:5]


#------#
# GLMS
#------#

# requires speicifying the exponential family modle: Poisson
# link funtion: log function, which is canonical link

Wren_GLM <- glm(cnt ~ (lon + lat + year) ^ 2 + ., 
    family = poisson("log"),
    data = Wren_df)

# since Poisson, the mean and varaiance is the same
# In cases where the variance in the data is greater
# than that suggested by this model
# the data are said to exhibit " over dispersion"

# An estimate of the dispersion is given by the ratio
# of the deviance t the df (the number of data points - the number of covariates)

## dispersion esimates
Wren_GLM$deviance / Wren_GLM$df.residual
# [1] 3.782881 > 1
# a sign of dispersion

# Another way to obtain estimate of the dispersion parameter
# to replace poisson with quasipoisson when calling glm
# then use summary
# the quasi-Poisson assumes the variance is propotionate 
# to the mean, then the constant of the proportionality is the over-dispersion parameter

# then need to use other models in the exponential family
# eg. negaitive-binomial distribution
# to account explicityly for the over-dispersion



#===========#
# Prediction
#===========#

# 1st generate sp-t prediction grid
# 80*80*21 grid in degrees * degrees * years
# covering the obs in sp and t

pred_grid <- expand.grid(lon = seq(min(MOcarolinawren_long$lon) - 0.2,
                      max(MOcarolinawren_long$lon) + 0.2,
                      length.out = 80),
            lat = seq(min(MOcarolinawren_long$lat) - 0.2,
                      max(MOcarolinawren_long$lat) + 0.2,
                      length.out = 80),
            year  = 1994:2014)


# Require all the covariate values at all the pred location
# Now evaluate the basis functions at the pred locations
S_pred <- eval_basis(basis = G,
           s = pred_grid[, c("lon", "lat")] %>%
             as.matrix()) %>%
  as.matrix()

colnames(S_pred) <- paste0("B", 1:ncol(S_pred))

# attach to grid
pred_grid <- cbind(pred_grid, S_pred)


# type = "link" to indicate that we predict
# the link function of the repsonse and not the response

wren_preds <- predict(Wren_GLM,
        newdata = pred_grid,
        type = "link",
        se.fit = TRUE)

# the predictions of the link function of the response
# are attached to our prediction grid for plotting

pred_grid <- pred_grid %>%
  mutate(log_cnt = wren_preds$fit,
         se = wren_preds$se.fit)


# when fitting GLMs, good practice to check the deviance residuals
# and inspect them for any residual correlation
# default GLM residuals are deviance residuals

Wren_df$residuals <- residuals(Wren_GLM)


#------#
# plot
#------#

g2 <- ggplot(Wren_df) + 
  geom_point(aes(lon, lat, colour = residuals)) +
  col_scale(name = "residuals") +
  facet_wrap(~ year, nrow = 3) +
  theme_bw()

plot(g2)

# noisy: indicating lack of spataial correlation

install.packages("ape")
library(ape)

P <- list()
years <- 1994:2014
for(i in seq_along(years)) {
  Wren_year <- filter(Wren_df,
         year == years[i])
  
  obs_dists <- Wren_year %>%
    select(lon, lat) %>%
    dist() %>%
    as.matrix()
  
  obs_dists_inv <- 1 / obs_dists
  diag(obs_dists_inv) <- 0
  
  P[[i]] <- Moran.I(Wren_year$residuals, 
          obs_dists_inv) %>%
    do.call("cbind", .)
  
}

do.call("rbind", P) %>% summary(digits = 2)


#-------------------------------------#
# Empirical semivarariogram of deviance
#-------------------------------------#

Wren_STIDF <- STIDF(sp = SpatialPoints(
                          Wren_df[, c("lon", "lat")],
                          proj4string = CRS("+proj=longlat")),
                    time = as.Date(Wren_df[, "year"] %>%
                                     as.character(),
                                   format = "%Y"),
                    data = Wren_df)


#--------------#
# Semivariogram
#--------------#

# consider time bins of width 1 year (of width 52.1429 weeks)
# Bins specified in units of weeks are required

install.packages("gstat")
library(gstat)


tlags <- seq(0.01, 52.1429 * 6 + 0.01, by = 52.1429)

vv <- variogram(object = residuals ~ 1, 
          data = Wren_STIDF,
          tlags = tlags,
          width  =  25, # spatial bin (25km)
          cutoff = 150, # use points < 150km apart
          tunit = "weeks")

plot(vv)
# given a certain location, not much variance in the time lat
# suggesting strong temporal corrlation in the residuals
# but huge difference in the distance
# indicating little spatial correlation

# a clear sign that a more sophisticated sp-t random effects model 
# should be considered for these data.







































