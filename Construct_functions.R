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
nrow(r)


#' @title construct_grid function
#' 

ngrid <- 41

construct_grid <- function(bbox, ngrid, coordnames = NULL) {
  ndim <- nrow(bbox)
  if(length(ngrid) == 1)
    ngrid <- rep(ngrid, ndim)
  if(!(length(ngrid) == ndim) | !is.numeric(ngrid))
    stop("ngrid needs to be a numeric (which will be rounded if not an integer)
         with length one or equal to the number of columns in bbox")
  ngrid <- round(ngrid)
  
  s <- lapply(1:nrow(bbox), function(i)
    seq(bbox[i,1], bbox[i,2], length.out = ngrid[i]))
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



construct_kernel <- function(Basis, ki) {
  if(!is.list(Basis)) stop("Basis needs to be of class list") # if a list
  if(!all(sapply(Basis, function(x) is(x,"Basis")))) # if each component is a Basis
    stop("All Basis functions need to be of class Basis")
  ndim <- dimensions(Basis[[1]]) # n_alpha
  
  function(s, r) {
    if(1) { ## Actual basis
      theta_s <- list()
      for(i in 1:(2 + ndim)) {
        theta_s[[i]] <- (eval_basis(Basis[[i]], s) %*% ki[[i]]) %>%
          as.numeric() %>%
          repcol(nrow(r))
      }
      theta_s_1 <- lapply(theta_s, function(x) x[,1])
      D <- FRK::distR(s + do.call("cbind", theta_s_1[3:(2 + ndim)]), r) # why distR btw . and r
      theta_s[[1]] * exp(-D^2/theta_s[[2]])
    } else {
      D <- FRK::distR(t(t(s) + c(ki[[3]], ki[[4]])), r)
      ki[[1]] * exp(-D^2/ki[[2]])
    }
  }
}


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








FRK::distR()

la <- list(1, "a", 2:4)
la[1]
la[2]
do.call("cbind", la[3]) + 0.1

A <- matrix(rnorm(50),5,10)
D <- distR(A,A[-3,])


