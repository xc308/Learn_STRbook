#***************************************#
# Lab2.3 Empirical Orthogonal Functions
#***************************************#

# empirical othorgonal functions can reveal 
# spatial structure in the data
# and can be used for subsequent dimensionality reduction


# EOFs can be obtained through either a spectral decomposition
# of the covariance matrix
# or a singular value decomposition (SVD)

# The data matrix has to be space-wide 
# (where space varies along the columns and time
# varies along the rows)

data("SSTlandmask", package = "STRbook")
data("SSTlonlat", package = "STRbook")
data("SSTdata", package = "STRbook")


install.packages("fields")
install.packages("CCA")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("gstat")
install.packages("sp")
install.packages("spacetime")


library(fields)
library(CCA)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gstat)
library(sp)
library(spacetime)
library(STRbook)




#============================#
# Preliminary data processing
#============================#

# since SSTdata contains readings over land
# we delete them using SSTlandmask
# also, consider whole years, so take the 33years
# that's 396 months of the data, constaining SST values
# spanning 1970-2002


delete_rows <- which(SSTlandmask == 1)

str(SSTdata) # 'data.frame':	2520 obs. of  399 variables:
SSTdata <- SSTdata[-delete_rows, 1:396]


# EOFs can reveal spatial sturcture in the data
# To calculate the EOFs, either spectral decomposition of the covariance matrix
# or SVD of the detrended space-time data matrix
# which has to be space-wide

str(SSTdata)
# 1st put this time-wide data sets into a space-wide one
dim(SSTdata) # 2261  396 Loaction by time 

Z <- t(SSTdata)
dim(Z) # [1]  396 2261



#=====================#
# Detrend and re-scale
#=====================#

spat_mean <- apply(SSTdata, 1, mean) # at each location, summing over time
nT <- ncol(SSTdata) # [1] 396 time points

Zspat_detrend <- Z - outer(rep(1, nT), spat_mean)

# standardize
Zt <-1/sqrt(nT - 1) * Zspat_detrend



#========#
# SVD Zt
#========#

E <- svd(Zt)
# E is a list contains U, V, D, 
# where V contains EOFs in space-wide format
# change the names of the col of V, and append lon. lat to it
# D diags of singular values

V <- E$v # 2261 396

colnames(V) <- paste0("EOF", 1:ncol(SSTdata))
head(V)

head(SSTlonlat)
EOFs <- cbind(SSTlonlat[-delete_rows, ], V)

head(EOFs[, 1:6])


# U contains the principle component time series in wide-tabel format
# i.e. each column corresponds to a time series associated with an EOFs

dim(E$u) # 396 396
U <- E$u
head(U)
str(U) # num [1:396, 1:396]


mutate() # add a new variable and preseves existing ones

TS<- data.frame(E$u) %>%
  mutate(t = 1:nrow(E$u)) %>%
  gather(key = EOF, value = PC, -t)








