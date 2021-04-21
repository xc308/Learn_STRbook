#***************************************#
# Lab 5.2 Sp-T inference using IDE model
#***************************************#

# use IDE package to fit sp-t IDE models & predict, forecast from sp_t data

# 3 cases:  
  # 1st two: simulated data where true model is known
  # 3rd: Sydney radar data


#=========#
# Package
#=========#

# IDE
# FRK: CONSTRUCT basis functions to model the spatially varying parameters of the kernel
# plyr: binding data frames with unequal column number

library("plyr")
library("dplyr")
install.packages("IDE")
library("IDE")

library("FRK")
library("ggplot2")
library("sp")
library("spacetime")
library("STRbook")

#=====================================================#
# simulation example with a spatially invariant kernel
#=====================================================#

# Package IDE contains a function simIDE 
  # that simulates the behavior of a typical dynamic system
  # governed by linear transport. Z = HY + e

  # can simulate from a user-defined IDE model or pre-defined one

  # In pre-defined model
    # T: # of time points
    # nobs: # of spatially fixed observations
    # k_spat_invariant = 1: spatially invariant kernel
    # k_spat_invariant = 0: 
    # includes a linear trend in s1 and s2

View(simIDE)


simIDE_self <- function(T = 9, obs = 100, k_spat_invariant = 1, IDEmodel = NULL) {
  
  ## suppress warning
  z <- val <- s1 <- s2 <- timeind <- NULL
  
  ## obs/spatial locations
  set.seed(21-4-2021)
  zlocs <- data.frame(s1 = runif(100), s2 = runif(100))
  
  ## spatial process decomposition
  Y_basis <- auto_basis(manifold = plane(), 
                        data = SpatialPoints(zlocs),
                        regular = 1,
                        nres = 2)
  
  r <- nbasis(Y_basis)
  
  ## invariant kernel decompositon
  G_const <- constant_basis()
  
  
  ## regression coefficients beta
  beta <- rep(0.2, 3)
  
  ## Other parameters sigmas
  sigma2_eps <- 0.01^2
  sigma2_eta <- 0.01^2
  
  ## spatial domain $ grid dots
  bbox<- matrix(c(0, 0, 1, 1), nrow = 2, byrow = F)
  s <- construct_grid(bbox = bbox, ngrid = 41)
  str(s$s_grid_df)
  
  ## rd coefficients matrix alpha
  alpha <- matrix(0, r, T)
  
  ## kernel basis 2 types
  if(k_spat_invariant) {
    K_basis <- list(G_const, G_const, G_const, G_const)
    k <- c(150, 0.002, -0.1, 0.1) # kernel parameters
    alpha[65, 1] <- 1 # 
  }else {
    G_var <- auto_basis(manifold = plane(), data = SpatialPoints(s$s_grid_df), 
                        regular = 1, nres = 1) # @df 9 obs 4 var
    nbk <- nbasis(G_var) #9
    K_basis <- list(G_const, G_const, G_var, G_var)
    k <- c(150, 0.002, 0.1 * rnorm(nbk), 0.1 * rnorm(nbk))
    alpha[sample(1:r, 10), 1] <- 1 
  }
  
  
  ## time map: temporal index & corresponding date
  T <- 9
  timemap <- data.frame(timeind = paste0("Y", 0:(T - 1)),
                        time = as.Date(0:(T - 1), origin = "2017-12-01"),
                        stringsAsFactors = F)
  
  
  ## construct Sigma matrices for process and data
  Sigma2_eps <- sigma2_eps * Diagonal(nobs * T)
  Sigma2_eta <- sigma2_eta * Diagonal(r)
  
  
  ## construct propogation matrix M
  Mfun <- construct_M(Y_basis = Y_basis, s)
  M <- Mfun(K_basis = K_basis, k)
  
  
  ## construct PHI decomposed spatial basis evaluated at spatial grid
  PHI <- eval_basis(Y_basis, s$s_grid_mat)
  str(PHI) # 1681 90
  
  
  ## Y0 = PHI %*% alpha0
  s$s_grid_df$Y0 <- (PHI %*% alpha[, 1]) %>% as.numeric()# 1681
  
  
  ## dynamical of alpha_(t) = M %*% alpha + t(PHI) %*% eta_t
  for(i in 1:(T - 1)) {
    alpha[, i + 1] <- (M %*% alpha[, i]) %>% as.numeric() + sqrt(sigma2_eta) * rnorm(r)
    s$s_grid_df[paste0("Y", i)] <-  (PHI %*% alpha[, i + 1]) %>% as.numeric()
    
  }

  
  ## Process in long format along with temporal data
  s_long <- gather(s$s_grid_df, key = timeind, value = val, -s1, -s2) %>%
    left_join(timemap, by = "timeind") %>% dplyr::select(-timeind)
  
  str(s_long$val) # num [1:1681] 0.307
  
  
  ## Covairates, intercept, s1, s2
  X_proc <- cbind(1, s_long[, c("s1", "s2")]) %>% as.matrix()
  str(X_proc) # num [1:1681, 1:3] 
  
  
  ## simulate data
    
    # fixed effects
  fixed_effects <- (X_proc %*% beta) %>% as.numeric()
  str(fixed_effects) # num [1:1681]
  
    # update process value with fixed effects
  
  str(s_long$val) # num [1:1681]
  s_long$val <- s_long$val + fixed_effects


    # df of obs locations rdly generated
  zlocs <- data.frame(s1 = runif(nobs), s2 = runif(nobs))
  
    # spat process evaluated at the observation locations PHI_obs_1
  PHI_obs_1 <- eval_basis(Y_basis, zlocs[, 1:2] %>% as.matrix())
  
    # PHI_obs_1 expand to include 9 days
  PHI_obs <- do.call("bdiag", lapply(1:T, function(x) PHI_obs_1))
  str(PHI_obs) # 450 810 sparse matrix sturcture 50 * 9 90 *9
  
    # X_obs expand to include 9 days
  X_obs <- cbind(1, do.call("rbind", lapply(1:T, function(x) zlocs))) %>% as.matrix
  str(X_obs) # num [1:450, 1:3]

    # Z = X_obs beta + PHI alpha + measurement error
  Z <- X_obs %*% beta + PHI_obs %*% c(alpha) + sqrt(sigma2_eps) * rnorm(nrow(PHI_obs))
  
  
    # create sp-t grid
  z_df <- data.frame(expand.grid(space = zlocs, time = data.frame(timemap$time)))

  
    # add one column of Z value above
  z_df$z <- Z %>% as.numeric()
  
  ## ggplot the observations facet by time
  g_obs <- ggplot(z_df) + geom_point(aes(s1, s2, color = z)) +
    facet_wrap(~ time) + 
    scale_color_distiller(palette = "Spectral") + 
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1))
  
  ## ggplot the process 
  g_truth <- ggplot(s_long) + geom_tile(aes(s1, s2, fill = val)) +
    facet_wrap(~ time) + 
    scale_color_distiller(pallete = "Spectral",
                          limits = c(min(z_df$z, s_long$val),
                                     max(z_df$z, s_long$val))) + 
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1))
  
  
}







