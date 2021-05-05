#**************************************#
# Lab 4.5 Non_Gau sp_t models with INLA
#**************************************#

# INLA is a Bayesian methods that provides approximate
# marginal posterior distribution over all states and parameters

# Predict expected counts at arbitrary sp-t locations
# from the vector of observed counts Z

library("INLA")
library("dplyr")
library("tidyr")
library("ggplot2")
library("STRbook")
data("MOcarolinawren_long", package = "STRbook")


head(MOcarolinawren_long, 2)
coords <- unique(MOcarolinawren_long[c("loc.ID", "lon", "lat")])
boundary <- inla.nonconvex.hull(as.matrix(coords[, 2:3]))
# a list of 4
str(boundary)


# the triangulation of the domain is using inla.mesh.2d
  # max.edgeg: two components (length 2 vector):
    # max edge length in interior domain
    # max edge length in the exterior of the domain
      # to reduce boundary effects

  # cutoff: min edge length
  # choose: max 0.8 for interior domain

MOmesh <- inla.mesh.2d(boundary = boundary,
             max.edge = c(0.8, 1.2),
             cutoff = 0.1)
# list of 8
str(MOmesh)


par(mai = c(rep(0.1, 4)))
plot(MOmesh, asp = 1, main = "")
lines(coords[c("lon", "lat")], col = "red", type = "p")


# In standard Gau, modelling effort lies in the
# establishing covariance matrix of a = (a1', a2', ...aT')'
# INLA, use the covariance matrix of a is sepberable of t and s
    # Sum_t(rho) krpdt Sum_s(tau, kappa, nu)
#s.t. the its inverse, i.e. precision matrix is sparse


# the matrix Sum_t(rho) is an AR(1) process
# so using single AR parameter rho
  # which dictates the correlation of a across time
  # the higher rho close to 1, the stronger temporal correlation


# the matrix Sum_s is parametrized using 3 paramters
# reflects the spatial covaraince required s.t. the 
# the reconstructed field is a solution to the stochastice
# SPDE
    #(kappa^2 - Laplacian)^(a/2)(tau Y(.)) = error(.)

  # where error(.) : white-noise process
  # tau: controls the variance
  
  # resulting field has Matern covarinace function
  # kappa : scaling parameter that translates to a 
      # "practical" spatial corrlation length 
    # i.e. spatial separation at which correlation is 0.1

  # l = sqrt(8 nu)/kappa, range parameter
  # a = nu + d/2 : smooth parameter
  # d: sp dimension

  # here, fix nu = 1, a = 2
  # this nu is notoriously difficult to estimate
  # and frequently set using CV


## SPDE can be constructed on the mesh using 
# the function inla.spde2.pcmatern

  # "pc" in pcmatern is penerlized complexity
  # used to refer to prior dist over the hyperparameters that
  # are both interprestable and have good theoretical properties
  
  # define prior dist over range parameter l st
    # P(l < 1) = 0.01
    # marginal std.deviation st P(sigma > 4) =0.01

spde <- inla.spde2.pcmatern(mesh = MOmesh,
                    alpha = 2, 
                    prior.range = c(1, 0.01),
                    prior.sigma = c(4, 0.01))

str(spde)


## 
# with mesh triangulation disretization
# a_{t,i} can be viewed as the weight of the ith basis fun
# at time t
# the obs matrix phi_t then maps the obs to finite-element space at time t

# if the obs lies exactly  on the vertex, then the associated row
# in Phi_t will be 0 everywhere except for a 1 in the col corresponding to the vertex

# or the row has 3 non-zero elements, with each representing
# the propotion being assigned to each vertex

# for point predictions or areal averages, all rows in Phi_t
# sum to 1



#========================#
# Create index usingINLA 
#========================#

# inla.spde.make.index
  # index name
  # the # of spatial points in the mesh
  # the # of time points

n_years <- length(unique(MOcarolinawren_long$t)) # 21
n_spatial <- MOmesh$n # 259 unique sp
s_index <- inla.spde.make.index(name = "spatial.field",
                     n.spde = n_spatial,
                     n.group = n_years)

str(s_index)
# List of 3
# $ spatial.field      : int [1:5439] 1 2 3 4 5 6 7 8 9 10 ...
# $ spatial.field.group: int [1:5439] 1 1 1 1 1 1 1 1 1 1 .

# spatial.field index: runs from 1 to n_spatial for n_year times
# spatial.field.group: runs from 1 to n_years, with each element replicated n_spatial times



#=============#
# Creating PHI
#=============#

# using inla.spde.make.A
  # mesh
  # loc: measurement locations
  # group: the measurement group year of obs
  # # of groups: n_years

coords.allyear <- MOcarolinawren_long[c("lon", "lat")] %>% as.matrix()

PHI <- inla.spde.make.A(mesh = MOmesh,
                 loc = coords.allyear,
                 group = MOcarolinawren_long$t,
                 n.group = n_years)

dim(PHI) # [1] 1575 5439 

nrow(MOcarolinawren_long) # [1] 1575 obs

length(s_index$spatial.field) # [1] 5439 index different sp within 1 time, and all the time

# so PHI is matrix of # obs * # indices of basis functions in sp and t


#===========================#
# Construct Latent Gau Model
#===========================#

# the latent Gau model is constructed in INLA throu a stack

# Stacks: allow one to define data, effects, obs matrices in groups
  # one accounting for measurement locations, 
  # the other accounting for prediction locations
  # which can be stack into a bigger stack
  # 

# A stack:
  # containging data and the covaraites at the measuremnet locations
  # is constructed by 
    # data: data
    # matrix PI:A
    # info on gamma
    # the label of stack is tagged
str(s_index)

n_data <- nrow(MOcarolinawren_long) # 1575

stack_est <- inla.stack(data = list(cnt = MOcarolinawren_long$cnt),
           A = list(PHI, 1),
           effects = list(s_index, list(Intercept = rep(1, n_data))),
           tag = "est")
str(stack_est)


#---------------------------------------------------#
# Stack containig matrices and vectors for prediction
#---------------------------------------------------#

# construct a stack containt the matrices and vectors
# defining the model at prediction locatiosn

# chosse the triangulation vertices as the prediction location
  # then PHI is simply the identity matrix
  # X: a vector of onces, no other weights
  # stroe the info on the prediction locations in df_pred

MOmesh$n #259 sp locations
head(MOmesh$loc, 2)

df_pred <- data.frame(lon = rep(MOmesh$loc[, 1], n_years),
           lat = rep(MOmesh$loc[, 2], n_years),
           t = rep(1:n_years, each = MOmesh$n))

n_pred <- nrow(df_pred) # 5439

PHI_pred <- Diagonal(n = n_pred)


#==================#
# Prediction stack
#==================#

# similar way to the estimation stack
# but set data values to NA to indicate predictions
# to be carried out at these locations

stack_pred <- inla.stack(data = list(cnt = NA),
           A = list(PHI_pred, 1), 
           effects = list(s_index, list(Intercept = rep(1, n_pred))),
           tag = "pred")


#==================================#
# Combine estimation and pred stack
#==================================#

stack <- inla.stack(stack_est, stack_pred)
str(stack)


#===============#
# Define fomula
#===============#

# a combination of a std R formula for the fix effects
# and an INLA formula for the sp-t residual component
  # which specify the name of the index as the 1st arg
  # the model: spde
  # the name of the grouping /time index: spatial.field.group
  # model constructed across the group: AR(1)

# choice for the prior on the AR(1) coefficient rho 
# is a pealized complexicty prior st P(rho > 0) = 0.9
# to reflect the prior belief that hightly doubt a neg teomp correlation


## PC penalty on rho
rho_hyper <- list(theta = list(prior = "pccor1",
                  param = c(0, 0.9)))



## fomula
formula <- cnt ~ -1 + Intercept + f(spatial.field, 
                         model = spde, 
                         group = spatial.field.group,
                         control.group = list(model = "ar1",
                                              hyper = rho_hyper))



#==============================#
# Fit the model and Predictions 
#==============================#

# needs data from stack extracted form inla.stack.data
# exponential family: neg-binomial

output <- inla(formula, 
     data = inla.stack.data(stack, spde = spde),
     family = "nbinomial",
     control.predictor = list(A = inla.stack.A(stack),
                              compute = TRUE))


## load directly from book
data("INLA_output", package = "STRbook")
  
  # INLA provides approximate marginal posterior distribution
  # for each a_t in a, and {beta, rho, tau, kappa}

  # returned objects, output, provides all results and summaries
  
  # from the posterior distributions over the precision par tau
  # and scale par kappa, can readily obtain marginal posterior
  # over the more interpretable varaiance par sigma^2 and 
  # practical range parameter l.

  # Posterior distri of some parameters can be seen 
    # AR(1) coefficient of the latent field rho is large
    # practical range par l is of order of 2 degress 
    # and the posterior dist of the marginal variance of the latent field is larage 2-4

  # all suggest strong spatial, temp dependencies in the data 


output.field <- inla.spde2.result(inla = output, 
                  name = "spatial.field",
                  spde = spde,
                  do.transform = TRUE)


#==============================#
# Plot posterior marginal distr
#==============================#

par(mai = c(rep(1, 4)))
par(mfrow = c(1, 1))

## plot P(beta0 | Z)
plot(output$marginals.fixed$Intercept,
     type = "l",
     xlab = expression(beta[0]),
     ylab = expression(p(beta[0]*"|"*Z)))


## plot P(rho | Z)
plot(output$marginals.hyperpar$`GroupRho for spatial.field`,
     type = "l",
     xlab = expression(rho),
     ylab = expression(p(rho*"|"*Z)))


## plot(sigma^2 | Z)
plot(output.field$marginals.variance.nominal[[1]],
     type = "l",
     xlab = expression(sigma^2),
     ylab = expression(p(sigma^2*"|"*Z)))


## plot(range | Z)
plot(output.field$marginals.range.nominal[[1]],
     type = "l",
     xlab = expression(l),
     ylab = expression(p(l*"|"*Z)))



#===================================#
# Prediction and std.e. of Prediction
#===================================#

# Prediction and associated std.error for log(Y(.))

# figures were generated by linerly interpolating the posterior mean
# and posterior std. deviationg of log(Y*) on a fine grid

# NOTE: 
  # how a high observed count at a certain location in one year
  # affects the predictions at the same location in neighbour years
  # even unobserved. 

# PLotting spatial fields from INLA output requires
  # each prediction and prediction std error of a_t for
  # each t needs to be projected spatially. 


  # 1st, extract the prediction and prediction std.error
# of a = (a1', a2',...., aT')

index_pred <- inla.stack.index(stack, 'pred')$data
lp_mean <- output$summary.fitted.values$mean[index_pred]
lp_sd <- output$summary.fitted.values$sd[index_pred]


#----------------#
# Create sp grid 
#----------------#

# upon which we map the predictions and std error
# construct 80*80 grid

grid_locs <- expand.grid(lon = seq(min(MOcarolinawren_long$lon) - 0.2,
                      max(MOcarolinawren_long$lon) - 0.2,
                      length.out = 80),
            lat = seq(min(MOcarolinawren_long$lat) - 0.2,
                      max(MOcarolinawren_long$lat) + 0.2,
                      length.out = 80))


## project
proj.grid <- inla.mesh.projector(MOmesh, 
                    xlim = c(min(MOcarolinawren_long$lon) - 0.2,
                             max(MOcarolinawren_long$lon) + 0.2),
                    ylim = c(min(MOcarolinawren_long$lat) - 0.2,
                             max(MOcarolinawren_long$lat) + 0.2),
                    dims = c(80, 80))


## Map each at onto sp grid. 
  # iterate through t, for each t = 1, ..., T
  # map both prediction and std. error of predictio of a_t on sp grid


  # t1: sp runs from 1:259
  # t2: sp runs from (259 + 1) : 2*259
  # t3: sp runs from (2*259 + 1) : 3*259
pred <- sd <- NULL
for(i in 1:n_years) {
  ii <- (i - 1) * MOmesh$n  +1
  jj <- i * MOmesh$n
  
  pred[[i]] <- cbind(grid_locs, 
        z = c(inla.mesh.project(proj.grid,
                                lp_mean[ii:jj])),
        t = i)
  
  sd[[i]] <- cbind(grid_locs, 
        z = c(inla.mesh.project(proj.grid,
                                lp_sd[ii:jj])),
        t = i)
  
}


#-------------------------#
# compile all data into df
#-------------------------#

# compile all the data which as in lists into 
# one df for plotting with ggplot2

# do.call(fun, arg) # arg of fun can be a list

pred <- do.call("rbind", pred) %>% filter(!is.na(z))
sd <- do.call("rbind", sd) %>% filter(!is.na(z))

str(pred)
range(pred$z)

range(sd$z)

#-------#
# Plot
#-------#

ggplot() + 
  geom_raster(data = pred, 
              aes(lon, lat, fill = z)) + 
  fill_scale(name = expression(alpha[t]), limits = c(-1, 5)) + 
  facet_wrap(~t, nrow = 3, ncol = 7) +
  theme_bw()


ggplot() + 
  geom_raster(data = sd, 
              aes(lon, lat, fill = z)) + 
  fill_scale(name = expression(sigma[t]), limits = c(0, 3)) + 
  facet_wrap(~t, nrow = 3, ncol = 7) + 
  theme_bw()










































































