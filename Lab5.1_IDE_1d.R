#*****************************************#
# Lab 5.1 Implementing an IDE in 1-d space
#*****************************************#

# Implement a stochastic integro-diff eq (IDE)
# in 1-d sp and t

# Y_t(s) = int_Ds m(s,x; theta_p) Y_(t-1) (x) dx + eta_t(s),

  # Y_t(s) : sp-t process a time t
  # theta_p: parameters fixed, to be estimated in practice
  # eta_t(.): spatial process, indpt of Y_t(.), with known cov functions


#=========#
# Packages
#=========#

library("dplyr")
library("ggplot2")
library("STRbook")
set.seed(14-04-2021)


#=========================================#
# Constructing the process grid and kernel
#=========================================#

# 1. constructing a discretization of the 1-d sp domain 
# Ds= [0, 1]
# 2. use this discretization, containing cells of width delta_s
# for approx integration and visualizations
# 3. call this discretization spatial grid


## spatial grid
ds<- 0.01
s_grid <- seq(0, 1, by = ds) # vector
N <- length(s_grid) # 101
str(s_grid) # num 1:101


## temp grid  (0~200)
nT <- 201
t_grid <- 0:(nT - 1) 
# not using seq(0, 200) as start from 1 not 0


## Sp-t grid
st_grid <- expand.grid(s = s_grid, t = t_grid) # 20301
head(st_grid)
tail(st_grid)


## transition kernel m(s,x; theta_p)
  # a bivariate function on sp grid
  # theta_p = (theta_1, theta_2, theta_3)'
    # theta_1: aplitude
    # theta_2: scale(aperture, 2*var)
    # theta_3: offset of the kernel

  # m(s,x; theta_p) = theta_1 * exp(-(x - theta_3 -s) / theta_2)


m <- function(s, x, thetap) {
  gamma <- thetap[1]  
  width <- thetap[2]
  offset <- thetap[3]
  
  D <- outer(s + offset, x, "-") # outer of a number and 1D array gets 2D array
  
  gamma * exp(-D ^ 2 / width)
}


#----------------#
# Understanding D
#----------------#

o <- outer(0.5, s_grid, "-") # 1* 101 row matrix, 2D array
o_vec_num <- outer(s_grid, 0.5, "-") # 101*1 col matrix, 2D array

head(o) # 2D arry, [1, 1:101], row matrix
head(outer(s_grid, 0.5, "-")) # 2D array, [1:101, 1], col matrix


## comparision
oo <- outer(o, o) # 4D array
str(oo)
head(oo)

o_sq <- o^2 # 1* 101 rwo matrix, 2D array
o_sq_u <- as.numeric(o_sq) # vector, 1D array


#-------------------------------------------------#
# Understding difference btw expand.grid and outer
#-------------------------------------------------#

x <- seq(0, 1)
y <- seq(2, 3)
z <- c("M", "F")
expand.grid(x = x, y = y, z = z) 
# the 1st factor fastest, y to folw, then z slowest

outer(x, y, "*") # x for row, y for col xy'



#============================#
# Visualize Kernel at s = 0.5
#============================#

# see how the process at s = 0.5 depends on x
# 4 kernels
  # 1st, narrow and center at .5
  # 2nd, slighlty wider
  # 3rd, shift to the right of .5
  # 4th, shift to the left of .5


thetap <- list()
thetap[[1]] <- c(40, 0.0002, 0)
thetap[[2]] <- c(5.75, 0.01, 0)
thetap[[3]] <- c(8, 0.005, 0.1)
thetap[[4]] <- c(8, 0.005, -0.1)

a <- m(s = 0.5, x = s_grid, thetap = thetap[[1]])
str(a)

m_x_0.5 <- m(s = 0.5, x = s_grid, thetap = thetap[[1]]) 
str(m_x_0.5)

nu <- as.numeric(m_x_0.5)


m_x <- list()
for(i in seq_along(thetap)) {
  m_x[[i]] <- m(s = 0.5, x = s_grid, thetap = thetap[[i]]) %>%
    as.numeric()
}


# Ref: allocate Gaussian kernel value at each sp grid to sp grid df
df <- data.frame(x = s_grid, m = m_x_0.5) # 101*2


df <- list()
for(i in seq_along(m_x)){
  df[[i]] <- data.frame(x = s_grid, m = m_x[[i]])
}

str(df)


plt <- list()
for(i in seq_along(df)){
  plt[[i]] <- ggplot(df[[i]]) + 
    geom_line(aes(x, m)) + 
    theme_bw()
}

print(plt)


#================#
# Define eta_t(.)
#================#

# define it as a sp process with an exponential cov function
# range parameter 0.1, and variance 0.1

Sigma_eta <- 0.1 * exp(-abs(outer(s_grid, s_grid, "-")) / 0.1)
# 2D array, [1:101, 1:101] matrix 


#==============================#
# Simulate eta_t(s) over s_grid
#==============================#

# by generating a multivariate Gau vector 
# with mean zero and cov matrix Sigma_eta


# One method:
  # use mvrnorm from MASS package
# Other method:
  # use lower Cholesky factor of Sigma_eta
  # and mutliply this by a vector of numbers generated from univariate rd normal
  # eta_t = Lz = L * rnrom(nrow(Sigam_eta))
  # so cov(eta_t) = cov(Lz)


chol_sigma_eta <- chol(Sigma_eta) # upper chol factor (0s upper)
L <- t(chol(Sigma_eta)) # to get lower chol factor (0s lower)
                        # matrix 101*101

sim <- L %*% rnorm(nrow(Sigma_eta)) # 2D array, col matrix, [1:101,1]
#z <- rnorm(nrow(Sigam_eta)) # 1D array, vector, 1:101
  # simulation from standard normal 

# since sim has the same 1st, 2nd moments with MVN eta
# simulation from Lwr_chol %*% simu_N(0,I) is the same as 
# the realization of 
# eta ~ MVN(0, Sigma_eta)


## 2nd method using mvrnorm in MASS
library("MASS")

sim2 <- mvrnorm(n = 1, 
                mu = rep(0, nrow(Sigma_eta)), 
                Sigma = Sigma_eta)
# multivariate distr, so 1 sample has the same number
# of elements as the length of mu



#==================================#
# Plot this realization of eta_t(.)
#==================================#

plot(s_grid, sim, "l")

plot(s_grid, sim2, "l")



#====================#
# Simulating from IDE
#====================#

# HAVE everything in place to simulate from the IDE
# using for loop
# 4 simulations, one for each kernel, stored in a list of df

# for each simulation setting i
  # simulate time points j 
   # to obtain the process

a <- t_grid[-1] # 1:200

## spatial grid
ds<- 0.01
s_grid <- seq(0, 1, by = ds) # vector
N <- length(s_grid) # 101


library("data.table")

Y <- list()
Yt <- list()
for(i in 1:4) {
  M <- m(s_grid, s_grid, thetap[[i]]) # construct kernel
  
  Y[[i]] <- data.frame(s = s_grid, 
             t = 0,    # init. time point
             Y = 0)    # init. process value 0
  
  for(j in t_grid[-1]) {
    prev_Y <- filter(Y[[i]], t == j-1)$Y  # Y at t-1
    
    eta <- L %*% rnorm(N)   # simulate eta rd sp effects at current time
    new_Y <- (M %*% prev_Y * ds + eta) %>% as.numeric() # from col matrix to 1D array, vector
    
    Y[[i]] <- rbind(Y[[i]], data.frame(s = s_grid, t = j, Y = new_Y))
    # appeding data frames repeatedly
    
    #Yt[[j]] <- data.frame(s = s_grid, t = j, Y = new_Y)
    #Y <- rbindlist(Yt)
    
    }
}



#=======================#
# Visualize 4 processes
#=======================#

# Y[[i]], i = 1,..., 4 contains df in long format
# straight forward to visualize

# plot Hovmoller plot for IDE process for i = 1

ggplot(Y[[1]]) + 
  geom_tile(aes(s, t, fill = Y)) + 
  scale_y_reverse() + 
  fill_scale(name = "process Y") + 
  theme_bw()



#========================#
# Simulating Observations
#========================#

# want to simulate noisy observations from process Y
# Because the only way to test whether algo for inference
# are working as they should is to mimic both
# underlying true process and the measurement process

# working with simulated data is the 1st step in developing 
# reliable algorithms to be applied to real data


# to map the observations of process to the data 
# need an incidence matrix that picks out the process value
# that has been observed. 

  # each row of the incidence matrix is composed of rows of observations of data
  # for each row, 0s everywhere except for the entry corresponding to
  # process value has been observed 

  # if locations of observation are changing overtime
  # incidence matrix changes overtime as well

# Assume at each time point
  # observe process at 50 locations
  # for easy, just a subset of s_grid
  # can also use nearest-neighbour mapping

nobs <- 50
sobs <- sample(s_grid, nobs)

# match the obs locations on the s_grid

Ht <- matrix(0, nobs, N)
for(i in 1:nobs) {
  idx <- which(sobs[i] == s_grid)
  Ht[i, idx] <- 1
}


#------------------#
# Simulate the data
#------------------#

# obs Zt = H_t * Y_t + e_t
  # e_t indpt Y_t
        # Gaussian rd vector whose entry are 
        # iid with mean zero and var sigma_e^2 = 1
  # H_t: same for each t
  # Y_t: latent process on the grid at time t


# ref:
## temp grid  (0~200)
nT <- 201
t_grid <- 0:(nT - 1) 


z_df <- data.frame()
for(j in 0:(nT-1)) {
  Yt <- filter(Y[[1]], t == j)$Y
  
  zt <- Ht %*% Yt + rnorm(nobs)
  
  z_df <- rbind(z_df, 
        data.frame(s = sobs, t = j, z = zt))
}


#--------------------------#
# Plot the simulated data z
#--------------------------#

ggplot(z_df) + 
  geom_point(aes(s, t, colour = z)) +
  col_scale(name = "simulated z") + 
  scale_y_reverse() + 
  theme_bw()

# noisy and sizeable gaps
# filling in these gaps by 1st estimating all parameters in IDE from data
# then predicting at unobserved locations. 















































