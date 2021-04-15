#****************************************#
# Lab 4.4 Non-Gaussina SpT GAMs with mgcv
#****************************************#

# aim to predict the expected counts at arbitray sp-t locations
# from the vector of observed counts Z

library("dplyr")
library("ggplot2")
library("mgcv")
library("STRbook")
library("tidyr")

data("MOcarolinawren_long", package = "STRbook")

## GAMs and GAMMs rely on constructing smooth functions
# of the covariates

# g(Y(s;t)) = b + f(s;t) + nu(s;t)
  # g(.) a link function
  # b intercept
  # f(s;t) is a random smooth function of sp, t
  # nu(s;t) sp t white-noise error process

# in mgcv, the rd function f(s;t) is generally decomposed using 
# a seperable spline basis

  # several basis function can be used to reconstruct f(s;t)
  # knot-based (B-spline)

# f(s;t) is decomposed as Sum phi_{1i}(s;t)a_{1i}
  # where {a_[1i]} are unknown rd effects need to be predicted
  # and {phi_{1i}(s;t)} are given below


# There are a number of basis functions can be chosen
# thin-plate regresssion splines are convenient
# easy to multiple covaraiates
  # thin-plate splines are isotropic and invariante to rotation
  # but not invariant to covariate scaling 
  # not recommand it for fitting curv over sp, t
  # since units in time are different from those in space

# To combine interacting covaraites with different units, sp, t
# mgcv uses a tensor-product structure
# whereby the basis functions smoothing the individual covariates
# are combined productwise

# f(s;t) = Sum_{i=1, r1} Sum{j=1, r2} phi_{1i}(s) phi_{2j}(t) a_{ij
        # = phi(s;t)' a


# the function te forms the product from the marginals
  # e.g. te(lon, lat, t)
    # bs: basis furntion class. 
      # eg. thin-plate over space, "tp", cubic regression over time, "cr"
    # k: the number of basis function
    # d: dimension of each spline 


f <- cnt ~ te(lon, lat, t,
         bs = c("tp", "cr"),
         k = c(50, 10), # 50 basis fun for space, 10 basis fun for t
         d = c(2, 1)) # space dim = 2, time dim = 1


# the number of knots could be set using cv
# the estimated df shall be considerably lower than total # of knots
# or, increase the # of knots



# As Carolina wren counts are over-dispersed
# to account fo rthis, use neg-binomial distr
# to model the response  a quasi-poission

cnts <- gam(f, family = nb(link = "log"),
    data = MOcarolinawren_long)


## to see the plots of the tensor plots
plot.gam(cnts)







# to predict unobserved locations using hierachical model
# 1st construct a sp-t grid upon which to predict
MOlon <- MOcarolinawren_long$lon
MOlat <- MOcarolinawren_long$lat

## construct grid
grid_locs <- expand.grid(lon = seq(min(MOlon) - 0.2,
                      max(MOlon) + 0.2,
                      length.out = 80),
            lat = seq(min(MOlat) - 0.2,
                      max(MOlat) + 0.2,
                      length.out = 80),
            t = 1:max(MOcarolinawren_long$t))


X <- predict(cnts, grid_locs, se.fit = TRUE)
# se.fit returns a list containing the predictions and pred std.error

grid_locs$pred <- X$fit
grid_locs$se <- X$se.fit


# plot predictions and overlay observations
ggplot() + 
  geom_raster(data = grid_locs, 
              aes(lon, lat, fill = pmin(pmax(pred, -1), 5))) + # plot predictions
  facet_wrap(~t, nrow = 3, ncol = 7) + 
  geom_point(data = filter(MOcarolinawren_long, !is.na(cnt)), # layover observation with black dots
             aes(lon, lat),
             colour = "black", size = 3) +
  geom_point(data = filter(MOcarolinawren_long, !is.na(cnt)), # lay on obs with values on log scale
             aes(lon, lat, colour = log(cnt)),
             size = 2) +
  fill_scale(name = expression(log(Y[t])), # prediction Y on log scale
             limits = c(-1, 5)) +
  col_scale(name = "log(cnt)", limits = c(-1, 5)) + # obs cnt on log 
  theme_bw()


## plot prediction std.errors
ggplot() + 
  geom_raster(data = grid_locs,
              aes(lon, lat, fill = pmin(se, 2.5))) + 
  facet_wrap(~t, nrow = 3, ncol = 7) + 
  fill_scale(name = expression(s.e.),
             palette = "BrBG",
             limits = c(0, 2.5)) + 
  theme_bw()




















