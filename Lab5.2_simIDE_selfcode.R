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
  
  
  ## spatial locations
  
  
  ## spatial process decomposition
  
  
  ## invariant kernel decompositon
  
  
  
  ## regression coefficients beta
  
  
  ## Other parameters sigmas
  
  
  ## spatial domain $ grid dots
  
  
  
  ## rd coefficients matrix alpha
  
  
  ## kernel basis 2 types
  
  
  ## time map: temporal index & corresponding date
  
  
  ## construct Sigma matrices for process and data
  
  
  ## construct propogation matrix M
  
  
  
  ## construct PHI decomposed spatial basis evaluated at spatial grid
  
  
  ## Y0 = PHI %*% alpha0
  
  
  ## dynamical of alpha_(t) = M %*% alpha + t(PHI) %*% eta_t
  
  
  ## Yt = PHI %*% alpha_t
  
  
  ## Process in long format along with temporal data
  
  
  ## Covairates, intercept, s1, s2
  
  
  
  ## simulate data
    
    # fixed effects
  
    # update process value with fixed effects
  
    # df of obs locations rdly generated
  
    # spat process evaluated at the observation locations PHI_obs_1
  
    # PHI_obs_1 expand to include 9 days
  
    # X_obs expand to include 9 days
  
    # Z = X_obs beta + PHI alpha + measurement error
  
    # create sp-t grid
  
    # add one column of Z value above
  

  
  ## ggplot the observations facet by time
  
  
  ## ggplot the process 
  
  
  ## 
  
  
}







