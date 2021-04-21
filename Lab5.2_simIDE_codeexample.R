simIDE <- function(T = 9, nobs = 100, k_spat_invariant = 1, IDEmodel = NULL) {
  ## Suppress bindings warning
  timeind <- val <- s1 <- s2 <- z <- NULL
  
  if(is.null(IDEmodel)) {
    
    set.seed(1)
    zlocs <- data.frame(s1 = runif(100),
                        s2 = runif(100))
    
    ## Spatial decomposition
    Y_basis <- auto_basis(manifold = plane(),
                          data = SpatialPoints(zlocs),
                          regular = 1,
                          nres = 2)  # Large Basis
    r <- nbasis(Y_basis)
    
    ## Kernel decomposition
    G_const <- constant_basis()
    
    ## Regression coeffocients
    beta <- c(0.2,0.2,0.2)
    
    ## Other parameters
    sigma2_eta <- 0.01^2
    sigma2_eps <- 0.01^2
    
    ## Spatial domain
    bbox <- matrix(c(0,0,1,1),2,2)
    s <- construct_grid(bbox, 41)
    
    alpha <- matrix(0,r,T)
    
    ## Kernel
    if(k_spat_invariant) {
      K_basis <- list(G_const, G_const, G_const, G_const)
      k <- list(150, 0.002, -0.1, 0.1)
      alpha[65,1] <- 1  
    } else {
      G <- auto_basis(plane(), data = SpatialPoints(s$s_grid_df),nres = 1)
      nbk <- nbasis(G) 
      K_basis <- list(G_const, G_const, G, G)
      k <- list(200, 0.002, 0.1*rnorm(nbk), 0.1*rnorm(nbk))
      alpha[sample(1:r,10),1] <- 1
    }
    
    time_map <- data.frame(timeind = paste0("Y",0:(T-1)),
                           time = as.Date(0:(T-1), origin = "2017-12-01"),
                           stringsAsFactors = FALSE)
  } else {
    Y_basis <- IDEmodel$get("process_basis")
    r <- nbasis(Y_basis)
    beta <- coef(IDEmodel)
    sigma2_eta <- c(IDEmodel$get("sigma2_eta"))
    sigma2_eps <- c(IDEmodel$get("sigma2_eps"))
    s <- IDEmodel$get("s")
    T <- IDEmodel$get("T")
    nobs <- nrow(IDEmodel$get("data"))
    K_basis <- IDEmodel$get("kernel_basis")
    k <- IDEmodel$get("k")
    alpha <- matrix(0,r,T)
    alpha[,1] <- sqrt(sigma2_eta) * rnorm(r)
    time_map <- data.frame(timeind = paste0("Y",0:(T-1)),
                           time = IDEmodel$get("time_points"),
                           stringsAsFactors = FALSE)
  }
  
  ## Construct matrices
  Sigma_eta <- sigma2_eta * Diagonal(r)
  Sigma_eps <- sigma2_eps * Diagonal(nobs * T)
  Q_eta <- Sigma_eta %>% solve()
  Q_eps <- Sigma_eps %>% solve()
  
  Mfun <- construct_M(Y_basis, s)
  M <- Mfun(K_basis, k)
  
  PHI <- eval_basis(Y_basis, s$s_grid_mat)
  
  s$s_grid_df$Y0 <- (PHI %*% alpha[,1]) %>% as.numeric()
  
  for(i in 1:(T-1)) {
    alpha[,i+1] <- (M %*% alpha[,i]) %>% as.numeric() + sqrt(sigma2_eta)*rnorm(r)
    s$s_grid_df[paste0("Y",i)] <- (PHI %*% alpha[,i+1]) %>% as.numeric()
  }
  
  
  ## process_value in long format
    ## head(s$s_grid_df): s1 s2 Y0 Y1 ... Y8
  s_long <- gather(s$s_grid_df, timeind, val, -s1, -s2) %>%
    left_join(time_map, by = "timeind") %>%
    select(-timeind)
  
  if(is.null(IDEmodel)) 
    X_proc <-  cbind(1, s_long[,c("s1","s2")]) %>% as.matrix()
  
  
  ## simulate data
  if(is.null(IDEmodel)) {
    fixed_effects <- (X_proc %*% beta) %>% as.numeric()
    
    s_long$val <- s_long$val + fixed_effects
    
    zlocs <- data.frame(s1 = runif(nobs), 
                        s2 = runif(nobs)) # nobs = 100
    
    PHI_obs_1 <- eval_basis(Y_basis, zlocs[,1:2] %>% as.matrix())
    PHI_obs <- do.call("bdiag", lapply(1:T, function(x) PHI_obs_1))
    
    X_obs <-  cbind(1, do.call("rbind", lapply(1:T, function(x) zlocs))) %>% as.matrix()
    
    Z <- X_obs %*% beta + PHI_obs %*% c(alpha) +
      sqrt(sigma2_eps) * rnorm(nrow(PHI_obs))
    
    z_df <- data.frame(expand.grid.df(zlocs, data.frame(time = time_map$time)))
    z_df$z <- Z %>% as.numeric()
  } else {
    fixed_effects <- 0
    s_long$val <- s_long$val + fixed_effects
    PHI_obs <- IDEmodel$get("PHI_obs")
    X_obs <- IDEmodel$get("X_obs")
    Z <- X_obs %*% beta + PHI_obs %*% c(alpha) +
      sqrt(sigma2_eps) * rnorm(nrow(PHI_obs))
    z_df <- as.data.frame(IDEmodel$get("data"))
    z_df[[all.vars(IDEmodel$get("f"))[1]]] <- Z %>% as.numeric()
  }
  
  
  g_obs <- ggplot(z_df) + geom_point(aes(s1, s2, colour = z)) +
    facet_wrap(~time) +
    scale_colour_distiller(palette = "Spectral")
  if(is.null(IDEmodel)) g_obs <- g_obs + coord_fixed(xlim=c(0,1), ylim = c(0,1))
  
  g_truth <- ggplot(s_long) + geom_tile(aes(s1,s2,fill=val)) +
    facet_wrap(~time) +
    scale_fill_distiller(palette="Spectral",
                         limits = c(min(c(z_df$z,s_long$val)),
                                    max(z_df$z,s_long$val)))
  if(is.null(IDEmodel)) g_truth <- g_truth + coord_fixed(xlim=c(0,1), ylim = c(0,1))
  
  ## Data as STIDF
  if(is.null(IDEmodel)) {
    cnames <- c("s1","s2")
    z_STIDF <- STIDF(sp = SpatialPoints(z_df[,cnames]),
                     time = z_df$time,
                     data = select(z_df, -time, -s1, -s2))
  } else {
    z_STIDF <- IDEmodel$get("data")
    z_STIDF$z <- as.numeric(Z)
  }
  
  
  ## IDEmode used to generate data
  if(is.null(IDEmodel)) {
    IDEmodel <- IDE(f = z ~ s1 + s2 + 1,
                    data = z_STIDF,
                    dt = as.difftime(1, units = "days"),
                    grid_size = 41,
                    kernel_basis = K_basis)
    IDEmodel$set(sigma2_eps = sigma2_eps,
                 sigma2_eta = sigma2_eta,
                 k = k)
  }
  
  list(s_df = s_long,
       z_df = z_df,
       z_STIDF = z_STIDF,
       g_truth = g_truth,
       g_obs = g_obs,
       IDEmodel = IDEmodel)
}

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
