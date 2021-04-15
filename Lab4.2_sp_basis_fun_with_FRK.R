#**************************************#
# Lab 4.2 Sp-t Basis Functions with FRK
#**************************************#

# using sp-t basis functions

library("dplyr")
library("FRK")
library("ggplot2")
library("gstat")
library("RColorBrewer")
library("sp")
library("spacetime")
library("STRbook")
library("tidyr")

#==========#
# FRK intro
#==========#

# FRK implements a low-rank approach to sp, sp-t modeling
# known as fixed rank kriging

# FRK considers the rd-effects model
# and provides functionality to the user for choosing
# the basis function from the data

# A key difference bwt FRK and other packaeg is 
# FRK modelling and prediction are carried out
# on a fine, regular discretization of the sp-t
# domain

# the small grid cells are basic areal unit
# primiary utility is to account for problems of COS
# (varying measurement footprint)

str(STObj3)
data("STObj3", package = "STRbook")


STOjb4 <- STObj3[, "1993-07-01::1993-07-31"]
STObj5 <- as(STOjb4[, -14], "STIDF")
STObj5 <- subset(STObj5, !is.na(STObj5$z))


#----------#
# SP-T BAU
#----------#

# Constructed using auto_BAUs
  # cellsize = c(1, 0.75, 1) 
  # 1 deg long * 1 deg lat * 1 day to ensure the BAUs are similar to prediction grid

  # convex: an extension radius used in domain construction via INLA
  # help(inla.nonconvex.hull)


#~~~~~~~~~~~~~#
# Install INLA
#~~~~~~~~~~~~~#

install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

library(INLA)

## EXAMPLE:
n = 100; a = 1; b = 1; tau = 100
z = rnorm(n)
eta = a + b*z

scale = exp(rnorm(n))
prec = scale*tau
y = rnorm(n, mean = eta, sd = 1/sqrt(prec))

data = list(y=y, z=z)
formula = y ~ 1+z
result = inla(formula, family = "gaussian", data = data)

summary(result)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``

# gridded BAUs arranged within a non-convex hull 
BAUs <- auto_BAUs(manifold = STplane(), # ST field on the plane
          type = "grid", # grided, not hex
          data = STObj5,
          cellsize = c(1, 0.75, 1),
          convex = -0.12, # hull extension
          tunit = "days")


# BAUs are class STFDF as they are 3-d pixels arranged regularly in sp, t

str(BAUs)
#SpatialPoints(proj4string = CRS(proj4string(STObj3)))


proj4string(BAUs)


# To plot the spatial BAUs overlaid with data location
par(mai = c(0.1, 0.1, 0.1, 0.1))

plot(as(BAUs[, 1], "SpatialPixels")) # choose day 1
plot(SpatialPoints(STObj5), add = TRUE, 
     col = "red")


# hexagonal BAUs

BAUs_hex <- auto_BAUs(manifold = STplane(),
          type = "hex", 
          data = STObj5,
          cellsize = c(1, 0.75, 1),
          nonconvex_hull = FALSE,
          tunit = "days")


# plot
plot(as(BAUs_hex[, 1], "SpatialPolygons"))
plot(SpatialPoints(STObj5), add = TRUE, col = "red")


#-------------------------#
# Construct basis function
#-------------------------#

# In FRK, these are constructed by tensor product
# of sp basis functions with temporal basis functions
  # a set of rs spatial basis functions 
  # a set of rt temporal basis functions

# generic basis function FRK uses by default is bisquare function

# basis functions can either be regularly placed or irregularly placed
# and often mutliresolution

# choose two resolutions
  # yielding rs = 94 spaial basis funcitons 
  # and place them irregularly in the domain
  # which are determined by auto_basis

G_spatial <- auto_basis(manifold = plane(),
           data = as(STObj5, "Spatial"),
           nres = 2, 
           type = "bisquare",
           regular = 0) # irregular


## temporal basis functions
  # use the function local_basis to construct regular seq of rt = 20
  # bisquare fucntiosn btw day 1 and day 31 of the month
  # each of these bisquare basis functions is assiganed an opening two days
  # so the support of each bisquare function is 4 days


# define temporal grid
t_grid <- matrix(seq(1, 31, length = 20))

# construct temporal basis function
G_temporal <- local_basis(manifold = real_line(),
            type = "bisquare",
            loc = t_grid, # location of center of cetroids
            scale = rep(2, 20)) # aperture paramter


## construct rsrt = 1880 = 94 * 20
G <- TensorP(G_spatial, G_temporal)

## PLOT each basis
show_basis(G_spatial)
show_basis(G_temporal)

## Note: 
  # While the basis functions are of tensor-product
  # the resulting ST covariance function obtained from the sp-t
  # rd effects model is not seperalbe in space and time


# In FRK, the fine-scale variation term at BAU level
# is assumed to be Gau with covariance matrix prop to diag(sig_nu,i ^2),
# i = 1, 2, ..., n_y
# where sig_nu,i ^2 are pre-specified at the BAU level
  # the constant of prop is then esmiated by FRK
  # related to some geographically related quantity e.g. surface roughness
  # here, we set = 1 for all i

BAUs$fs <- 1 # fine-scale


# the fine-scale variance at the BAU level is confounded with 
# measure-ment error variance 
  # some cases, measurement error is known
  # if not known, just the value of semivariog at the origin, partial sill
  # so simply assume the nugget effect is the measurement error variance

  # thus, any residual nugget compent is assumed to be fine-scale varaince introducd
  # as a consequence of the low-rank approximation to the processes


#sepVgm$space
#model      psill    range
#1   Nug 0.04910985   0.0000
#2   Exp 0.95089015 465.7563

STObj5$std <- sqrt(0.049)


## use lat as covariante and response is z
f <- z ~ lat + 1


## now ready to call the main FRK,
  # which estimates all the unknown parameters in the models
  # incld: covariance matrix  of the basis-function coeffience Ca
  # and the fine-scale variance sigma_nu ^2



S <- FRK(f = z ~ lat + 1,     # formular
    data = list(STObj5), # list of STIDF
    basis = G,   # sp-t basis
    BAUs = BAUs, # gridded BAUs
    n_EM = 3,   # max no. of EM iterations
    tol = 0.01) # convergence tolerance


grid_BAUs <- predict(S)






































































