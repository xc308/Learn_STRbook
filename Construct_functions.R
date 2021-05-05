## Initialize
alphahat <- M <- Q <- Q_eps <- Q_eta <- k <- betahat <-
  Qpost <- Qpostchol <- sigma2_eps <- sigma2_eta <- NULL
G_const <- new("Basis",
               manifold = plane(),
               fn = list(function(s) rep(1, nrow(s))),
               pars = list(),
               df = data.frame(),
               n = 1)



#' @title Create a single, constant basis function
#' @return Object of class \code{Basis}

constant_basis <- function() {
  new("Basis",
      manifold = plane(),
      fn = list(function(s) rep(1, nrow(s))),
      pars = list(),
      df = data.frame(),
      n = 1)
}


r <- nbasis(Y_basis) # 90
nrow(r) # NULL


#' @title construct_grid function
#' 


#-------------------
ngrid <- 41
ndim <- nrow(bbox)
ngrid <- rep(ngrid, ndim)
s <- lapply(1:nrow(bbox), function(i)
  seq(bbox[i,1], bbox[i,2], length.out = ngrid[i]))

str(s)
#----------------------------

construct_grid <- function(bbox, ngrid, coordnames = NULL) {
  ndim <- nrow(bbox)
  if(length(ngrid) == 1)
    ngrid <- rep(ngrid, ndim)
  if(!(length(ngrid) == ndim) | !is.numeric(ngrid))
    stop("ngrid needs to be a numeric (which will be rounded if not an integer)
         with length one or equal to the number of columns in bbox")
  ngrid <- round(ngrid)
  
  s <- lapply(1:nrow(bbox), function(i) seq(bbox[i,1], bbox[i,2], length.out = ngrid[i]))
  
  s_grid <- do.call("expand.grid", s)
  if(is.null(coordnames)) {
    names(s_grid) <- paste0("s",1:ndim)
  } else {
    names(s_grid) <- coordnames
  }
  list(s_grid_df = s_grid,
       s_grid_mat = s_grid %>% as.matrix(),
       ds = (bbox[,2] - bbox[,1])/ngrid,
       area = prod((bbox[,2] - bbox[,1])/ngrid))
}



repcol <- function(x,n){
  l <- lapply(1:n, function(i) x)
  y <- do.call("c", l)
  matrix(y, ncol = n, byrow = FALSE)
}




#---------------------------
K_basis <- list(G_const, G_const, G_const, G_const)
(eval_basis(K_basis[[1]], s$s_grid_mat) %*% k[[1]]) %>%
  as.numeric() %>%
  repcol(length(r)) # turn into col matrix

e <- eval_basis(K_basis[[1]], s$s_grid_mat)  %*% k[[1]] %>% as.numeric() %>% repcol(length(r))
str(e)

ndim <- dimensions(K_basis[[1]])  # 2

str(K_basis[[1]])


theta_s <- list()
for(i in 1:(2 + ndim)) {
  theta_s[[i]] <- (eval_basis(K_basis[[i]], s$s_grid_mat) %*% k[[i]]) %>%
    as.numeric() %>%
    repcol(length(r)) # turn into col matrix 
}

str(theta_s)
# List of 4
#$ : num [1:1681, 1] 150 150 150 150 150 150 150 150 150 150 ...
#$ : num [1:1681, 1] 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002 ...
#$ : num [1:1681, 1] -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 ...
#$ : num [1:1681, 1] 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 ...

theta_s_1 <- lapply(theta_s, function(x) x[,1])
str(theta_s_1)
#$ : num [1:1681] 150 150 150 150 150 150 150 150 150 150 ...
#$ : num [1:1681] 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002 ...
#$ : num [1:1681] -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 ...
#$ : num [1:1681] 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 ...


cbind_theta34 <- do.call("cbind", theta_s_1[3:(2 + ndim)])
str(cbind_theta34) # num [1:1681, 1:2]
head(cbind_theta34, 4)
#      [,1] [,2]
#[1,] -0.1  0.1
#[2,] -0.1  0.1
#[3,] -0.1  0.1
#[4,] -0.1  0.1


head(s$s_grid_mat, 4)
#        s1 s2
#[1,] 0.000  0
#[2,] 0.025  0
#[3,] 0.050  0
#[4,] 0.075  0

# horizontal vertical offsets
sum_s_theta34 <- do.call("cbind", theta_s_1[3:(2 + ndim)]) + s$s_grid_mat # horizontal vertical offsets

d_of_r <- distR(sum_s_theta34)
str(d_of_r) num [1:1681, 1:1681] 

D <- distR(sum_s_theta34)
D2 <- D^2
str(D2)

str(theta_s[[2]])
D_scaled <- -D^2/theta_s[[2]][1, 1]
exp_Dscale <- exp(-D^2/theta_s[[2]][1, 1])
str(exp_Dscale) # 1681:1681


kernel <- theta_s[[1]][1, 1] * exp(-D^2/theta_s[[2]][1, 1])
str(kernel)


View(distR)


#---------------------------
str(s$s_grid_mat) # num [1:1681, 1:2]
head(s$s_grid_mat) 
#         s1 s2
#[1,]  0.000  0
#[2,]  0.025  0 

head(s$s_grid_df)
#      s1 s2
# 1 0.000  0
# 2 0.025  0

construct_M <- function(Y_basis, s) {
  PHI <- eval_basis(Y_basis, s$s_grid_mat)
  GRAM <- crossprod(PHI)*s$area
  GRAM_inv <- solve(GRAM)
  ndim <- dimensions(Y_basis)
  
  function(K_basis, ki) {
    K <- construct_kernel(K_basis, ki)
    Kmat <- K(s$s_grid_mat, s$s_grid_mat)
    M <- GRAM_inv %*% crossprod(t(Kmat) %*% PHI, PHI)*s$area^2
  }
}


construct_kernel <- function(Basis, ki) {
  if(!is.list(Basis)) stop("Basis needs to be of class list")
  if(!all(sapply(Basis, function(x) is(x,"Basis"))))
    stop("All Basis functions need to be of class Basis")
  ndim <- dimensions(Basis[[1]])
  
  function(s, r) {
    if(1) { ## Actual basis
      theta_s <- list()
      for(i in 1:(2 + ndim)) {
        theta_s[[i]] <- (eval_basis(Basis[[i]], s) %*% ki[[i]]) %>%
          as.numeric() %>%
          repcol(nrow(r))
      }
      theta_s_1 <- lapply(theta_s, function(x) x[,1])
      D <- FRK::distR(s + do.call("cbind", theta_s_1[3:(2 + ndim)]), r)
      theta_s[[1]] * exp(-D^2/theta_s[[2]])
    } else {
      D <- FRK::distR(t(t(s) + c(ki[[3]], ki[[4]])), r)
      ki[[1]] * exp(-D^2/ki[[2]])
    }
  }
}






#' @title Create a single, constant basis function
#' @description Constructs an object of class \code{Basis} as defined in \code{FRK} that is constant over the entire spatial domain.
#' @return Object of class \code{Basis}
#' @seealso \code{\link{IDE}} for how to use basis functions to construct the IDE kernel
#' @export
#' @examples
#' basis1 <- constant_basis()
constant_basis <- function() {
  new("Basis",
      manifold = plane(),
      fn = list(function(s) rep(1, nrow(s))),
      pars = list(),
      df = data.frame(),
      n = 1)
}



get <- function(obj) {
  switch(obj, "sigma2_eps" = sigma2_eps,
         "sigma2_eta" = sigma2_eta,
         "Q_eps" = Q_eps,
         "Q_eta" = Q_eta,
         "alphahat" = alphahat,
         "betahat" = betahat,
         "coordnames" = coordnames,
         "data" = data,
         "PHI_obs" = PHI_obs,
         "plausible_ranges" = plausible_ranges,
         "process_basis" = process_basis,
         "kernel_basis" = kernel_basis,
         "Qpost" = Qpost,
         "time_points" = time_points,
         "X_obs" = X_obs,
         "Q" = Q,
         "M" = M,
         "k" = k,
         "f" = f,
         "nk"= nk,
         "r" = r,
         "s" = s,
         "m" = m,
         "T" = T,
         "Z" = Z)
}




FRK::distR()

la <- list(1, "a", 2:4)
la[1]
la[2]
do.call("cbind", la[3]) + 0.1

A <- matrix(rnorm(50),5,10)
D <- distR(A,A[-3,])


