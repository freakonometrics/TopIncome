#' Density of the Pareto 1 distribution
#'
#' @param x a (positive) vector
#' @param mu a number (the lower bound)
#' @param alpha a number (the tail index)
#' @return the density of the Pareto 1 distribution at points \code{x}
#' @examples
#' dpareto1(2, 1, 1.5)
dpareto1 <- function (x, mu, alpha)  { alpha*mu^alpha/x^(alpha+1) }

#' Cumulative distribution function of the Pareto 1 distribution
#'
#' @param x a (positive) vector
#' @param mu a number (the lower bound)
#' @param alpha a number (the tail index)
#' @return the c.d.f. of the Pareto 1 distribution at points \code{x}
#' @examples
#' ppareto1(2, 1, 1.5)
ppareto1 <- function (x, mu, alpha)  { 1 - ( x/mu )^(-alpha) }

#' Quantile function of the Pareto 1 distribution
#'
#' @param p a vector of probabilities (with values in [0,1])
#' @param mu a number (the lower bound)
#' @param alpha a number (the tail index)
#' @return the quantile function of the Pareto 1 distribution at points \code{p}
#' @examples
#' qpareto1(.5, 1, 1.5)
qpareto1 <- function (p, mu, alpha)  { mu*(1-p)^(-1/alpha) }

#' Random generation of the Pareto 1 distribution
#'
#' @param n an integer
#' @param mu a number (the lower bound)
#' @param alpha a number (the tail index)
#' @return generates \code{n} values of the Pareto 1 distribution
#' @examples
#' set.seed(123)
#' rpareto1(6, 1, 1.5)
rpareto1 <- function (n, mu, alpha)  { mu*(1-runif(n))^(-1/alpha) }

#' Maximum Likelihood estimation of the Pareto 1 distribution, with weights
#'
#' @param data a vector of observations
#' @param weights a vector of weights (default = 1)
#' @param threhold the threshold parameter of the Pareto 1 distribution (\code{mu})
#' @return a list with the index \code{alpha} and \code{k}, the number of observations above \code{threshold}
#' @examples
#' set.seed(123)
#' x <- rpareto1(100, 1, 1.5)
#' w <- rgamma(100,10,10)
#' estim <- MLE.pareto1(data=x, weights=w, threshold=1)
#' estim
MLE.pareto1 <- function(data, weights=rep(1,length(x)), threshold=min(x))
{
  foo=cbind(data,weights)
  foo=foo[foo[,1]>threshold,]
  xx=foo[,1]
  ww=foo[,2]/sum(foo[,2])
  m <- as.numeric(threshold)
  a <- 1/(sum(ww*log(xx))-log(m))
  k <- NROW(xx)
  return(list(alpha=a,k=k))
}

#' Density of the Generalized Pareto distribution (GPD)
#'
#' @param x a (positive) vector
#' @param xi a number (the tail index)
#' @param mu a number (the lower bound)
#' @param beta a number (the scaling paramater, default = 1)
#' @return the c.d.f. of the Generalized Pareto distribution at points \code{x}
#' @examples
#' pgpd(2, 1/1.5, 1, 1)
pgpd <- function (x, xi, mu = 0, beta = 1)  { (1 - (1+(xi*(x-mu))/beta)^(-1/xi)) }

#' Density of the Generalized Pareto distribution (GPD)
#'
#' @param x a (positive) vector
#' @param xi a number (the tail index)
#' @param mu a number (the lower bound)
#' @param beta a number (the scaling paramater, default = 1)
#' @return the density of the Pareto 1 distribution at points \code{x}
#' @examples
#' dgpd(2, 1/1.5, 1, 1)
dgpd <- function (x, xi, mu = 0, beta = 1)  { (beta^(-1))*(1+(xi*(x-mu))/beta)^((-1/xi)-1) }

#' Random generation of the Generalized Pareto distribution (GPD)
#'
#' @param n an integer
#' @param xi a number (the tail index)
#' @param mu a number (the lower bound)
#' @param beta a number (the scaling paramater, default = 1)
#' @return generates \code{n} values of the Pareto 1 distribution at points \code{x}
#' @examples
#' rgpd(10, 1/1.5, 1, 1)
rgpd <- function (n, xi, mu = 0, beta = 1)  { mu + (beta/xi)*((1-runif(n))^(-xi)-1) }

#' Maximum Likelihood estimation of the Generalized Pareto distribution, with weights
#'
#' @param data a vector of observations
#' @param weights a vector of weights (default = 1)
#' @param threhold the threshold parameter of the Generalized Pareto distribution (\code{mu})
#' @param nextrmes the number of largest values considered (integer)
#' @param method method used for inference (\code{"ml"} for maximum likelihood)
#' @param information (not used)
#' @return a list with \code{n} the (total) number of observations, \code{threshold} the threshold, \code{p.less.thresh}, \code{n.exceed}, \code{k}, \code{method}, \code{converged}, \code{nllh.final} and \code{par.ests} a named vector with \code{"xi"} the tail index and \code{"beta"} the scaling coefficient
#' @examples
#' set.seed(123)
#' x <- rpareto1(100, 1, 1.5)
#' w <- rgamma(100,10,10)
#' estim <- MLE.gpd(data=x, weights=w, threshold=1)
#' estim$par.ests
MLE.gpd <- function (data, weights=rep(1,length(x)), threshold = NA, nextremes = NA, method="ml", information = c("observed", "expected"), ...) 
{
  n <- length(data)
  if (is.na(nextremes) && is.na(threshold)) 
    stop("Enter either a threshold or the number of upper extremes")
  if (!is.na(nextremes) && !is.na(threshold)) 
    stop("Enter EITHER a threshold or the number of upper extremes")
  if (!is.na(nextremes)) 
    
    data <- as.numeric(data)
  
  foo=cbind(data,weights)
  foo=foo[foo[,1]>threshold,]
  x=foo[,1]
  w=foo[,2]/sum(foo[,2])

  exceedances <- x
  excess <- exceedances - threshold
  Nu <- length(excess)
  xbar <- sum(w*excess)
  method <- "ml"
  s2 <- sum(w*(excess-xbar)^2)
  xi0 <- -0.5 * (((xbar * xbar)/s2) - 1)
  beta0 <- 0.5 * xbar * (((xbar * xbar)/s2) + 1)
  theta <- c(xi0, beta0)
  negloglik <- function(theta, tmp) {
    xi <- theta[1]
    beta <- theta[2]
    cond1 <- beta <= 0
    cond2 <- (xi <= 0) && (max(tmp) > (-beta/xi))
    if (cond1 || cond2) 
      f <- 1e+06
    else {
      y <- logb(1 + (xi * tmp)/beta)
      y <- w*y/xi
      f <- logb(beta) + (1 + xi) * sum(y)
    }
    f
  }
  fit <- optim(theta, negloglik, hessian = TRUE, ..., tmp = excess)
  if (fit$convergence) 
    warning("optimization may not have succeeded")
  par.ests <- fit$par
  converged <- fit$convergence
  nllh.final <- fit$value
  p.less.thresh <- 1 - Nu/n
  out <- list(n = length(data), threshold = threshold, 
              p.less.thresh = p.less.thresh, n.exceed = Nu, k= Nu, method = method, 
              par.ests = par.ests, converged = converged, nllh.final = nllh.final)
  names(out$par.ests) <- c("xi", "beta")
  return(out)
}

.EPDinput <- function(y, gamma, kappa, tau, kappaTau = TRUE) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #
    if (!is.numeric(gamma)) {
        stop("gamma should be numeric.")
    }
    
    if (!is.numeric(kappa)) {
        stop("kappa should be numeric.")
    }
    
    if (!is.numeric(tau)) {
        stop("tau should be numeric.")
    }
    

    if (any(tau >= 0)) {
        stop("tau should be strictly negative.")
    }
    
    if (any(gamma <= 0)) {
        stop("gamma should be strictly positive.")
    }
    if (kappaTau) {
        if (any(kappa <= pmax(-1, 1/tau))) {
            stop("kappa should be larger than max(-1,1/tau).")
        }
    }

    ly <- length(y)
    lg <- length(gamma)
    lk <- length(kappa)
    lt <- length(tau)
    
    l <- c(ly, lg, lk, lt)

    ind <- which(l > 1)
    
    if (length(ind) > 1) {

        if (!length(unique(l[ind])) == 1) {
            stop("All input arguments should have length 1 or equal length.")
        }
    }
}

#' Density of the Extended Pareto distribution
#'
#' @param x a (positive) vector
#' @param gamma a (strictly positive) number (the tail index)
#' @param kappa a number - must be larger than max{-1,1/tau}
#' @param tau a (negative) number (default is -1)
#' @log logical indicating if logarithm of density should be returned
#' @return the density of the Extended Pareto distribution at points \code{x}
#' @source \url{https://github.com/TReynkens/ReIns/blob/master/R/EPD.R} Tom Reynkens, ReIns package version 1.0.7
#' @examples
#' depd(2,.5,1,-1)
depd <- function(x, gamma, kappa, tau = -1, log = FALSE) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #

    .EPDinput(x, gamma, kappa, tau, kappaTau = TRUE)

    d <- 1 / (gamma*x^(1/gamma+1)) * (1+kappa*(1-x^tau))^(-1/gamma-1) *
    (1+kappa*(1-(1+tau)*x^tau))
    d[x <= 1] <- 0
    if (log) d <- log(d)
    return(d)
}

#' Cumulative Distribution Function of the Extended Pareto distribution
#'
#' @param x a (positive) vector
#' @param gamma a (strictly positive) number (the tail index)
#' @param kappa a number - must be larger than max{-1,1/tau}
#' @param tau a (negative) number (default is -1)
#' @log logical indicating if logarithm of density should be returned
#' @return the c.d.f. of the Extended Pareto distribution at points \code{x}
#' @source \url{https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R} Tom Reynkens, ReIns package version 1.0.7
#' @examples
#' pepd(2,.5,1,-1)
pepd <- function(x, gamma, kappa, tau = -1, lower.tail = TRUE, log.p = FALSE) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #
    
    .EPDinput(x, gamma, kappa, tau, kappaTau = FALSE)
    
    p <- 1 - (x * (1+kappa*(1-x^tau)))^(-1/gamma)
  
    p[x <= 1] <- 0

    if (any(kappa <= pmax(-1, 1/tau))) {
        if (length(kappa) > 1 | length(tau) > 1) {
            p[kappa <= pmax(-1, 1/tau)] <- NA
        } else {
            p <- NA
        }
    }
    
    if (!lower.tail) p <- 1-p
    
    if (log.p) p <- log(p)
    
    return(p)
}

#' Quantile Function of the Extended Pareto distribution
#'
#' @param p a vector of probabilities (in the interval [0,1])
#' @param gamma a (strictly positive) number (the tail index)
#' @param kappa a number - must be larger than max{-1,1/tau}
#' @param tau a (negative) number (default is -1)
#' @log logical indicating if logarithm of density should be returned
#' @return the c.d.f. of the Extended Pareto distribution at points \code{x}
#' @source \url{https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R} Tom Reynkens, ReIns package version 1.0.7
#' @examples
#' qepd(.5,.5,1,-1)
qepd <-  function(p, gamma, kappa, tau = -1, lower.tail = TRUE, log.p = FALSE) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #
    
    .EPDinput(p, gamma, kappa, tau, kappaTau = TRUE)
    
    
    if (log.p) p <- exp(p)
    
    if (!lower.tail) p <- 1-p
    
    if (any(p < 0 | p > 1)) {
        stop("p should be between 0 and 1.")
    }
    
    l <- length(p)
    Q <- numeric(l)
  
    endpoint <- 10
    
    if (any(p < 1)) {
        
        mx <- max(p[p < 1])
        
        while (pepd(endpoint, gamma, kappa, tau) <= mx) {
            endpoint <- endpoint*10
        }
    }

    for (i in 1:l) {
        
        if (p[i] < .Machine$double.eps) {
            # p=0 case
            Q[i] <- 1
            
        } else if (abs(p[i]-1) > .Machine$double.eps) {
            # 0<p<1 case
            
            # Function to minimise
            f <- function(x) {
                ((1-p[i])^(-gamma) - x*(1+kappa*(1-x^tau)))^2
            }
            # If minimising fails return NA
            Q[i] <- tryCatch(optimise(f, lower=1, upper=endpoint)$minimum, error=function(e) NA)
            
        } else {
            # p=1 case
            Q[i] <- Inf
        }
        
    }
    
    return(Q)
}


#' Random Generation of the Extended Pareto distribution
#'
#' @param n integer, number of generations
#' @param gamma a (strictly positive) number (the tail index)
#' @param kappa a number - must be larger than max{-1,1/tau}
#' @param tau a (negative) number (default is -1)
#' @param log logical indicating if logarithm of density should be returned
#' @return a vector of \code{n} values generated from an Extended Pareto distribution
#' @source \url{https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R} Tom Reynkens, ReIns package version 1.0.7
#' @examples
#' set.seed(123)
#' repd(6,.5,1,-1)
repd <-  function(n, gamma, kappa, tau = -1) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #
    
    return(qepd(runif(n), gamma=gamma, kappa=kappa, tau=tau))
}

#' Hill estimator of the tail index, with weights
#'
#' @param data the vector of observations
#' @param weights the vector of weights
#' @return Hill estimator of \code{gamma} (inverse of \code{alpha})
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- rpareto1(100, 1, 1.5)
#' w <- rgamma(100,10,10)
#' Hill(x,w)
# }
Hill = function(data,weights=rep(1,length(data))){
    w <- weights/sum(weights)
    n <- length(data)
    X <- as.numeric(sort(data))
    Hill <- sum(w[2:n]*(log(X[2:n])-log(X[1])))
    return(list(gamma = Hill))
}

#' Fit the Extended Pareto distribution to a vector of observations, with weights
#'
#' @param data vector of observations
#' @param weights vector of (positive) weights
#' @param rho parameter of Fraga Alves et al. (2003) estimate
#' @param start vector of length 2 containing the starting values for the optimisation
#' @param direct logical indicating if the parameters are obtained by directly maximising the log-likelihood function
#' @param warnings logical indicating if possible warnings from the optimisation function are shown
#' @return a list with \code{k} the vector of the values of the tail parameter, \code{gamma} the vector of the corresponding estimates for the tail parameter of the EPD, \code{kappa} the vector of the corresponding MLE estimates for the kappa parameter of the EPD and \code{tau} the vector of the corresponding estimates for the second order tail index parameter of the EPD using Hill estimates and values for \code{rho}
#' @source adapted from \url{https://github.com/TReynkens/ReIns/blob/master/R/EPD.R} Tom Reynkens, ReIns package version 1.0.7
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- repd(100,.5,1,-1)
#' w <- rgamma(100,10,10)
#' EPD(data=x, weight=w)
#' k
#' [1] 99
#'
#' gamma
#' [1] 0.3843927
#'
#' kappa
#' [1] 0.1819718
#'
#' tau
#' [1] -3.34675
#' }
EPD <- function(data, weights=rep(1,length(data)), rho = -1, start = NULL, direct = TRUE, warnings = FALSE, ...) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #
    
    df=data.frame(data,weights)
    df=df[order(df$data),]
    w <- df$weights/sum(df$weights)
    X <- as.numeric(df$data)
    n <- length(X)
    K <- (n-1)
    
    if (n == 1) {
        stop("We need at least two data points.")
    }
    
    if (direct) {
        EPD <- .EPDdirectMLE(data=data, weight=w, rho=rho, start=start, warnings=warnings)
    } else {
        # Select parameter using approach of Beirlant, Joosens and Segers (2009).
        EPD <- .EPDredMLE(data=data, weight=w, rho=rho)
    }
    
    if (length(rho) == 1) {
        EPD$gamma <- as.vector(EPD$gamma)
        EPD$kappa <- as.vector(EPD$kappa)
        EPD$tau <- as.vector(EPD$tau)
        
    }
    
    if (length(rho) == 1) {
        return(list(k=K, gamma=EPD$gamma[K], kappa=EPD$kappa[K], tau=EPD$tau[K]))
    } else {
        return(list(k=K, gamma=EPD$gamma[K,], kappa=EPD$kappa[K,], tau=EPD$tau[K,]))
    }
}


.EPDredMLE <- function(data, weights=rep(1,length(data)), rho = -1) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    #   original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #
    #   Fit EPD using approach of Beirlant, Joosens and Segers (2009).
    
    df=data.frame(data,weights)
    df=df[order(df$data),]
    w <- df$weights/sum(df$weights)
    X <- as.numeric(df$data)
    
    n <- length(X)
    K <- (n-1)
    
    if (n == 1) {
        stop("We need at least two data points.")
    }
    
    nrho <- length(rho)
    rho.orig <- rho
    
    H <- Hill(data, w)$gamma
    
    if (all(rho > 0) & nrho == 1) {
        rho <- .rhoEst(data, alpha=1, tau=rho, w=w)$rho
        beta <- -rho
        
    } else if (all(rho < 0)) {
        beta <- -rho
        
    } else {
        stop("rho should be a single positive number or a vector (of length >=1) of negative numbers.")
    }
    
    
    gamma <- matrix(0, n-1, nrho)
    kappa <- matrix(0, n-1, nrho)
    tau <- matrix(0, n-1, nrho)
    
    beta <- -rho
    
    for (j in 1:nrho) {
        
        if (nrho == 1 & all(rho.orig > 0)) {
            tau[K, 1] <- -beta[K]/H
            
            
            rhovec <- rho
        } else {
            
            tau[K, j] <- -beta[j]/H
            
            rhovec <- rep(rho[j], n-1)
        }
        
        
        E <- numeric(n-1)
        
        for (k in K) {
            i <- 1:k
            E[k] <- sum( w[n-k+i] * (X[n-k+i]/X[n-k])^tau[k,j] )
        }
        
        kappa[K,j] <- H * (1-2*rhovec[K]) * (1-rhovec[K])^3 / rhovec[K]^4 * (E[K] - 1 / (1-rhovec[K]))
        
        gamma[K,j] <- H - kappa[K,j] * rhovec[K] / (1 - rhovec[K])
        
    }
    
    return(list(gamma=gamma, kappa=kappa, tau=tau))
}

.EPDdirectMLE <- function(data, weights=rep(1,length(data)), rho = -1, start = NULL,  warnings = FALSE) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #
    
    
    df=data.frame(data,weights)
    df=df[order(df$data),]
    w <- df$weights/sum(df$weights)
    X <- as.numeric(df$data)
    n <- length(X)
    
    if (n == 1) {
        stop("We need at least two data points.")
    }
    
    nrho <- length(rho)
    rho.orig <- rho
    
    H <- Hill(data, w)$gamma
    
    if (all(rho > 0) & nrho == 1) {
        rho <- .rhoEst(data, alpha=1, tau=rho, w=w)$rho
        beta <- -rho
        
    } else if (all(rho < 0)) {
        beta <- -rho
        
    } else {
        stop("rho should be a single positive number or a vector (of length >=1) of negative numbers.")
    }
    
    gamma <- matrix(0, n-1, nrho)
    kappa <- matrix(0, n-1, nrho)
    tau <- matrix(0, n-1, nrho)
    
    for (j in 1:nrho) {
     
        for (k in (n-1):(n-1)) {
            
            epddf <- df[df$data > X[n-k],]
            epddata <- epddf$data/X[n-k]
            epdw    <- epddf$w
            
            if (nrho == 1 & all(rho.orig > 0)) {
                tau[k,1] <- -beta[k]/H
                
            } else {
                
                tau[k,j] <- -beta[j]/H
            }
            
            if (is.null(start)) {
                start2 <- numeric(2)
                start2[1] <- H
                start2[2] <- 0
            } else if (is.matrix(start)) {
                
                if (nrow(start >= n-1)) {
                    start2 <- numeric(2)
                    start2[1] <- start[k,1]
                    start2[2] <- start[k,2]
                } else {
                    stop("start does not contain enough rows.")
                }
                
            } else {
                start2 <- start
            }
            
            if (tau[k,j] < 0) {
                tmp <- EPDfit(epddata, start=start2, tau=tau[k,j], w=epdw)
                gamma[k,j] <- tmp[1]
                kappa[k,j] <- tmp[2]
            } else {
                gamma[k,j] <- kappa[k,j] <- NA
            }
            
        }
        
    }
    
    return(list(gamma=gamma, kappa=kappa, tau=tau))
}

#' Fit the Extended Pareto distribution to a vector of observations, with weights, using maximum likelihood estimation
#'
#' @param data vector of observations
#' @param tau the value for tau in the EPD distribution
#' @param weight vector of (positive) weights
#' @param rho parameter of Fraga Alves et al. (2003) estimate
#' @param start vector of length 2 containing the starting values for the optimisation (default are \code{c(.1,1})
#' @param warnings logical indicating if possible warnings from the optimisation function are shown
#' @return a vector with \code{gamma} the vector of the corresponding estimates for the tail parameter of the EPD, \code{kappa} the vector of the corresponding MLE estimates for the kappa parameter of the EPD
#' @source adapted from \url{https://github.com/TReynkens/ReIns/blob/master/R/EPD.R} Tom Reynkens, ReIns package version 1.0.7
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- repd(100,.5,1,-1)
#' w <- rgamma(100,10,10)
#' EPDfit(data=x, tau=-3.3, weight=w)
#' [1] 0.32989522 0.08325996
#' }
EPDfit <- function(data, tau, start = c(0.1, 1), warnings = FALSE, weights=rep(1,length(data))) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #
    
    w <- weights/sum(weights)
    if (is.numeric(start) & length(start) == 2) {
        gamma_start <- start[1]
        kappa_start <- start[2]
    } else {
        stop("start should be a 2-dimensional numeric vector.")
    }
    
    
    if (ifelse(length(data) > 1, var(data) == 0, 0)) {
        sg <- c(NA, NA)
    } else {

        fit <- optim(par=c(gamma_start, kappa_start), fn=.EPDneglogL, Y=data, tau=tau, w=w)

        sg <- fit$par
        
        if (fit$convergence > 0 & warnings) {
            warning("Optimisation did not complete succesfully.")
            if (!is.null(fit$message)) {
                print(fit$message)
            }
        }
    }
    return(sg)
}


.EPDneglogL <- function(theta, Y, tau, weights) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #
    
    w <- weights/sum(weights)
    gamma <- theta[1]
    kappa <- theta[2]
    
    if (kappa <= max(-1, 1/tau) | gamma <= 0) {
        logL <- -10^6
    } else {
        logL <- sum( w*log(depd(Y, gamma=gamma, kappa=kappa, tau=tau)) )
    }

    return(-logL)
}

.rhoEst <- function(data, alpha = 1, theta1 = 2, theta2 = 3, tau = 1, weights=rep(1,length(data))) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #
    
    
    if (alpha <= 0) {
        stop("alpha should be strictly positive.")
    }
    
    if (tau <= 0) {
        stop("tau should be strictly positive.")
    }
  
    df=data.frame(data,weights)
    df=df[order(df$data),]
    w <- df$weights/sum(df$weights)
    X <- as.numeric(df$data)
    
    n <- length(X)
    rho <- numeric(n)
    Tn <- numeric(n)
    K <- (n-1)
    
    M_alpha <- numeric(n)
    M_alpha_theta1 <- numeric(n)
    M_alpha_theta2 <- numeric(n)
    
    l <- log(X[n-K+1])
    for (k in K) {
        M_alpha[k] <- sum( (l[1:k]-log(X[n-k]))^alpha ) / k
        M_alpha_theta1[k] <- sum( (l[1:k]-log(X[n-k]))^(alpha*theta1) ) / k
        M_alpha_theta2[k] <- sum( (l[1:k]-log(X[n-k]))^(alpha*theta2) ) / k
    }
    
    Tn[K] <- ( (M_alpha[K]/gamma(alpha+1))^tau - (M_alpha_theta1[K]/gamma(alpha*theta1+1))^(tau/theta1)  ) /
    ( (M_alpha_theta1[K]/gamma(alpha*theta1+1))^(tau/theta1) - (M_alpha_theta2[K]/gamma(alpha*theta2+1))^(tau/theta2)  )
    
    rho[K] <- 1 - ( 2 * Tn[K] / ( 3 - Tn[K]) ) ^ (1/alpha)
    
    return(list(k=K, rho=rho[K], Tn=Tn[K]))
}


ProbEPD <- function(data, q, gamma, kappa, tau, ...) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #
    
    
    if ( length(gamma) != length(kappa) | length(gamma) != length(tau)) {
        stop("gamma, kappa and tau should have equal length.")
    }
    
    X <- as.numeric(sort(data))
    
    
    n <- length(X)
    prob <- numeric(n)
    K <- (n-1)
    
    K2 <- K[which(gamma[K] > 0)]
    
    prob[K2] <- (K2+1)/(n+1) * (1 - pepd(q/X[n-K2], gamma=gamma[K2], kappa=kappa[K2], tau=tau[K2]))
    prob[prob < 0 | prob > 1] <- NA
  
    
    return(list(k=K, P=prob[K], q=q))
    
}

#' Large Return Period associated to the Extended Pareto distribution
#'
#' @param data a vector of observations
#' @param q the used large quantile - to estimate 1/P[X>q]
#' @param gamma vector of \code{n-1} estimates for the EVD obtained from [EPD]
#' @param kappa vector of \code{n-1} estimates for the EVD obtained from [EPD]
#' @param tau vector of \code{n-1} estimates for the EVD obtained from [EPD]
#' @return a list with \code{k} the vector of the values of the tail parameter k, \code{R} the vector of the corresponding return period and \code{q} the used large quantile
#' @source \url{https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R} Tom Reynkens, ReIns package version 1.0.7
ReturnEPD <- function(data, q, gamma, kappa, tau, ...) {
    #
    #   Modified EPD function from the ReIns R package to have:
    #        - weighted ML estimation
    #        - results from only one cutoff
    #        - direct ML estimation by default
    # original R code :  ReIns package version 1.0.7
    #    https://github.com/TReynkens/ReIns/blob/master/R/EPD.R
    #    https://github.com/TReynkens/ReIns/blob/master/R/Distributions.R
    #
    
    if ( length(gamma) != length(kappa) | length(gamma) != length(tau)) {
        stop("gamma, kappa and tau should have equal length.")
    }
    
    X <- as.numeric(sort(data))
    n <- length(X)
    r <- numeric(n)
    K <- (n-1)
    
    K2 <- K[which(gamma[K] > 0)]
    
    r[K2] <- (n+1)/(K2+1) / (1 - pepd(q/X[n-K2], gamma=gamma[K2], kappa=kappa[K2], tau=tau[K2]))
    
    
    r[which(gamma[K] <= 0)] <- NA
    
    r[r <= 0] <- NA
    
    return(list(k=K, R=r[K], q=q))
    
}

#' Estimate Top Share
#'
#' @param data a vector of observations
#' @param weight a vector of weights
#' @param p (default \code{0.1})
#' @param q (default \code{0.1})
#' @param method (default \code{"edf"})
#' @param edp.direct logical (default \code{TRUE})
#' @return blah blah
TopShare <- function(data, weights=rep(1,length(x)), p=.1, q=.1, method="edf", epd.direct=TRUE) {
    require(Hmisc)
    #
    # p : top 100p% share
    # q : top 100q% of the distribution being Pareto
    #
    # method="edf" : sample top share, that is, based on the EDF
    # method="pareto1" : EDF+Pareto1
    #
    x = data
    if (p>1) stop('Error: p should be smaller than 1 \n\n')
    if (p<0) stop('Error: p should be greater than 0 \n\n')
    up=Hmisc::wtd.quantile(x, weights=weights, probs=1-p, normwt=TRUE) # weighted (1-p)-quantile
    up=as.numeric(up)
    
    ## Top Share based on the Empirical Distribution Function (EDF)
    
    if(method=="edf") {
        tis=sum(weights*x*(x>up))/sum(weights*x)
        return(list(index=tis,share=p,method="edf"))
    }
    
    ## Top Share based on Pareto I and GPD Models
    
    if (q>1) stop('Error: q should be smaller than 1 \n\n')
    if (q<0) stop('Error: q should be greater than 0 \n\n')
    
    
    u=Hmisc::wtd.quantile(x, weights=weights, probs=1-q, normwt=TRUE)
    u=as.numeric(u)  # threshold = weighted (1-q)-quantile
    data=cbind(x,weights)
    dataq=data[x<u,]
    xq = Hmisc::wtd.mean(dataq[,1], weights=dataq[,2])
    
    if(method=="pareto1" || method=="pareto2" || method=="gpd") {
        # Estimate the Pareto distribution with weighted data
        if(method=="pareto1") {
            coef=MLE.pareto1(x,w=weights,threshold=u)
            sigma=u
            alpha=coef$alpha
        }
        if(method=="pareto2" || method=="gpd") {
            coef=MLE.gpd(x,w=weights,threshold=u)
            sigma=coef$par.ests["beta"]/coef$par.ests["xi"]
            alpha=1/coef$par.ests["xi"]
        }
        # Top Income Shares with weighted data
        if (alpha<1) { tis=NaN
        } else if(p<=q) {
            num = (alpha/(alpha-1))*sigma*(p/q)^(-1/alpha) + u - sigma
            den = (1-q)*xq + q*sigma/(alpha-1) + q*u
            tis = p*num/den
        } else if(p>q)  {
            up=Hmisc::wtd.quantile(x, weights=weights, probs=1-p, normwt=TRUE)
            up=as.numeric(up)
            datap=data[x<up,]
            xp = Hmisc::wtd.mean(datap[,1], weights=datap[,2])
            den = (1-q)*xq + q*sigma/(alpha-1) + q*u
            tis = 1 - (1-p)*xp/den
        }
        return(list(index=tis,alpha=alpha,coef=coef,share.index=p,share.pareto=q,threshold=u))
    }
    
    ## Top Share based on EPD Model
    if(method=="epd") {
        
        dataqq=data[x>=u,]
        coef=EPD(dataqq[,1], w=dataqq[,2], direct=epd.direct)
        delta=coef$kappa
        tau=coef$tau
        alpha=1/coef$gamma
        
        up=Hmisc::wtd.quantile(x, weights=weights, probs=1-p, normwt=TRUE)
        up=as.numeric(up)
        datap=data[x<up,]
        xp = Hmisc::wtd.mean(datap[,1], weights=datap[,2])
        
        pextpareto=function(x, u=1, delta=0, tau=-1, alpha){#Compute CDF
            d=1-((x/u)*(1+delta-delta*(x/u)^tau))^(-alpha)
            d[x<=u] <- 0
            return(d) }
        ff_bis <- function(x) (1-pextpareto(1/x, u=u, delta=delta, tau=tau, alpha=alpha))/x^2
        
        if (alpha<=1) {tis=NaN    # infinite mean (tail=alpha=1/xi < 1)
        } else if(delta<max(-1,1/tau)) {tis=NaN   # kappa should be largen than max(-1,1/tau)
        } else if(p<=q) {
            uprim=u*qepd(1-p/q, gamma=coef$gamma, kappa=coef$kappa, tau=coef$tau)
            Eup = try( integrate(ff_bis, lower=0, upper=1/uprim)$value , TRUE)
            Eu = try( integrate(ff_bis, lower=0, upper=1/u)$value , TRUE)
            if (inherits(Eup, "try-error") && inherits(Eup, "try-error")) tis=NaN
            else tis=(p*uprim+q*Eup)/((1-q)*xq+q*(u+Eu))
        } else if(p>q)  {
            Eu = try( integrate(ff_bis, lower=0, upper=1/u)$value , TRUE)
            Ex = ((1-q)*xq+q*(u+Eu))
            if (inherits(Eu, "try-error")) tis=NaN
            else tis = 1-(1-p)*xp/Ex
        }
        return(list(index=tis,alpha=1/coef$gamma,coef=coef,share.index=p,share.pareto=q,threshold=u))
    }
    
}

#' Convert income/wealth Data
#'
#' @param income a vector of data (income or wealth)
#' @param weight a vector of weight (same length as \code{income}
#' @return a dataframe with 4 columns, \code{y} the vector income (or wealth), \code{weights} the vector of weights, \code{Fw} the cumulated proportion of people (with weights) and \code{Fx} the cumulated proportion of people (without weights)
tidy_income <- function(income, weights){
df=data.frame(w=weights, y=income)
df$w=df$weights/sum(df$weights)
df=df[order(df$y),]
Fw=cumsum(df$weights)/(sum(df$weights)+df$weights[1])
n=length(df$y)
Fx=(1:n)/(n+1)
data = data.frame(y=df$y, weights=df$w, Fw=Fw, Fx=Fx)
return(data)
}

#' Pareto diagrams - Pareto 1, GPD and EPD
#'
#' @param data dataframe obtained from \code{tidy_income} function
#' @param p numeric, the probability level (default 1%)
#' @param q numeric, the probability level to model a Pareto distribution (default 10% - top 10%)
#' @param viz logical \code{TRUE} to plot the estimates
#' @return a table with estimations of top share and a graph
Pareto_diagram = function(data, p=.01, q=.1, viz=TRUE){

res1=TopShare(data$y, weights=data$weights, p=p, q=q, method="pareto1")
res2=TopShare(data$y, weights=data$weights, p=p, q=q, method="gpd")
res3=TopShare(data$y, weights=data$weights, p=p, q=q, method="epd", epd.direct=epd.direct)

pot=data[data$y>0,]   # Keep positive data

if(viz) par(mfrow=c(1,1), mar=c(4, 4, 4, 1))  # bottom, left, top, right
if(viz) plot(log(pot$y), log(1-pot$Fw), main=mytitle, xlab="log(x)", ylab="log(1-F(x))", cex=.6, col="gray", xlim=PDxlim)

u=seq(log(res1$threshold), 30, length.out=500)
yhat.par1=ppareto1(exp(u),mu=res1$threshold,alpha=res1$coef$alpha)
yhat.par2=pgpd(exp(u),xi=res2$coef$par.ests["xi"],mu=res2$coef$threshold,beta=res2$coef$par.ests["beta"])
yhat.epd=pepd(exp(u)/res3$threshold,gamma=res3$coef$gamma,kappa=res3$coef$kappa,tau=res3$coef$tau)
if(viz){
lines(u,log(1-yhat.par1)+log(q), col="blue", lty=2, lwd=1.5)
lines(u,log(1-yhat.epd)+log(q), col="red", lty=1, lwd=1.5)
lines(u,log(1-yhat.par2)+log(q),col="green", lty=3, lwd=1.5)
legend("topright", legend=c("Pareto 1", "GPD", "EPD"), col=c("blue","green", "red"), lty=c(2,3,1))
}

# plot percentile as vertical dashed lines

res90=TopShare(data$y, weights=data$weights, p=p, q=.10, method="pareto1")
if(viz) abline(v=log(res90$threshold), col="lightgrey", lty=2)  # percentile 90
#legend(log(res90$threshold)-top.x, top.y, legend=c("top10%"), cex=.82, bty="n")
if(viz) legend(log(res90$threshold)-top.x, top.y, legend=expression(italic('q')[90]), cex=.9, bty="n")

res95=TopShare(data$y, weights=data$weights, p=p, q=.05, method="pareto1")
if(viz) abline(v=log(res95$threshold), col="lightgrey", lty=2)  # percentile 95
if(viz) legend(log(res95$threshold)-top.x, top.y, legend=expression(italic('q')[95]), cex=.9, bty="n")

if(viz) res99=TopShare(data$y, weights=data$weights, p=p, q=.01, method="pareto1")
if(viz) abline(v=log(res99$threshold), col="lightgrey", lty=2)  # percentile 99
legend(log(res99$threshold)-top.x, top.y, legend=expression(italic('q')[99]), cex=.9, bty="n")
}

#' Table of top shares (using three thresholds)
#'
#' @param data dataframe obtained from \code{tidy_income} function
#' @p probability level (default 1%)
#' @param q1 numeric, the probability level to model a Pareto distribution (default 10% - top 10%)
#' @param q2 numeric, the probability level to model a Pareto distribution (default 5% - top 5%)
#' @param q3 numeric, the probability level to model a Pareto distribution (default 1% - top 1%)
Table_Top_Share = function(data, p=.01, q1=.1 , q2=.05 , q3=.01){

res90=TopShare(data$y, weights=data$weights, p=p, q=q1, method="pareto1")
res95=TopShare(data$y, weights=data$weights, p=p, q=q2, method="pareto1")
res99=TopShare(data$y, weights=data$weights, p=p, q=q3, method="pareto1")
pareto1.index=cbind(res90$index, res95$index, res99$index)
pareto1.alpha=cbind(res90$alpha, res95$alpha, res99$alpha)

res90=TopShare(data$y, weights=data$weights, p=p, q=q1, method="pareto2")
res95=TopShare(data$y, weights=data$weights, p=p, q=q2, method="pareto2")
res99=TopShare(data$y, weights=data$weights, p=p, q=q3, method="pareto2")
gpd.index=cbind(res90$index, res95$index, res99$index)
gpd.alpha=cbind(res90$alpha, res95$alpha, res99$alpha)

res90=TopShare(data$y, weights=data$weights, p=p, q=q1, method="epd")
res95=TopShare(data$y, weights=data$weights, p=p, q=q2, method="epd")
res99=TopShare(data$y, weights=data$weights, p=p, q=q3, method="epd")
epd.index=cbind(res90$index, res95$index, res99$index)
epd.alpha=cbind(res90$alpha, res95$alpha, res99$alpha)

cutoff=c(1-q1,1-q2,1-q3)
cat("----- index ----------\n")
M=rbind(cutoff,pareto1.index,gpd.index,epd.index)
colnames(M)=c("index1","index2","index3")
cat(M,"\n")
cat("----- alpha ----------\n")
M=rbind(cutoff,pareto1.alpha,gpd.alpha,epd.alpha)
colnames(M)=c("alpha1","alpha2","alpha3")
cat(M,"\n")
cat("----- top share ------\n")
    T=TopShare(data$y, weights=data$w, p=p)
cat(T)
    return(T)}

#' Top Income plot
#'
#' @param data dataframe obtained from \code{tidy_income} function
#' @p probability level (default 1%)
#' @param thr numeric vector of probability levels to model a Pareto distribution (from 85% up to 99.9%)
#' @param TSlim numeric 2-vector, range of y for the plot (default \cote{NULL})
#' @param tail logical to plot the tail index (default \cote{TRUE})
Top_Income = function(data, p=.01, thr=seq(.85,.999,by=.001), TSlim=NULL, tail = TRUE){

thr=round(thr,10)
tail=matrix(0,NROW(thr),7)
tis.index=matrix(0,NROW(thr),7)
tis.alpha=matrix(0,NROW(thr),7)
for(i in 1:NROW(thr)) {
    
    res1=TopShare(data$y, weights=data$weights, p=p, q=1-thr[i], method="pareto1")
    res2=TopShare(data$y, weights=data$weights, p=p, q=1-thr[i], method="gpd")
    res3=TopShare(data$y, weights=data$weights, p=p, q=1-thr[i], method="epd", epd.direct=epd.direct)
    res4=TopShare(data$y, weights=data$weights, p=p, method="edf")
    
    tis.index[i,1]=res1$threshold     # threshold y0
    tis.index[i,2]=res1$coef$k          # k largest observations
    tis.index[i,3]=thr[i]             # quantile threshold
    tis.index[i,4]=res1$index
    tis.index[i,5]=res2$index
    tis.index[i,6]=res3$index
    tis.index[i,7]=res4$index
    
    tis.alpha[i,1]=res2$threshold           # threshold y0
    tis.alpha[i,2]=res2$coef$k          # k largest observations
    tis.alpha[i,3]=thr[i]             # quantile threshold
    tis.alpha[i,4]=res1$alpha
    tis.alpha[i,5]=res2$alpha
    tis.alpha[i,6]=res3$alpha
    tis.alpha[i,7]=0
    
}

if(tail){
plot(tis.alpha[,2],tis.alpha[,4], ylim=c(0,ysup), type="b", cex=.75, pch=3, main="MLE estimates of the tail index", xlab="k largest values", ylab="tail index (alpha)", col="blue")
lines(tis.alpha[,2],tis.alpha[,4], col="blue", type="l", cex=.75)
lines(tis.alpha[,2],tis.alpha[,5], col="green", type="p", cex=.75, pch=2)
lines(tis.alpha[,2],tis.alpha[,5], col="green", type="l", cex=.75)
lines(tis.alpha[,2],tis.alpha[,6], col="red", type="b", cex=.75, pch=1)
lines(tis.alpha[,2],tis.alpha[,6], col="red", type="l", cex=.75)
abline(v=tis.alpha[(tis.alpha[,3]==.90),2], col="lightgray", lty=2) # 10% top obs
abline(v=tis.alpha[(tis.alpha[,3]==.95),2], col="lightgray", lty=2) #  5% top obs
abline(v=tis.alpha[(tis.alpha[,3]==.99),2], col="lightgray", lty=2) #  1% top obs
legend("topright", legend=c("Pareto 1 (Hill estimator)","GPD", "EPD"), col=c("blue", "green", "red"), pch=c(3,2,1), lty=1)

legend(tis.alpha[(tis.alpha[,3]==.90),2]-top.xx,top.yy, legend=expression(italic('q')[90]), cex=.9, bty="n")
legend(tis.alpha[(tis.alpha[,3]==.95),2]-top.xx,top.yy, legend=expression(italic('q')[95]), cex=.9, bty="n")
legend(tis.alpha[(tis.alpha[,3]==.99),2]-top.xx,top.yy, legend=expression(italic('q')[99]), cex=.9, bty="n")
}

if(is.null(TSlim)) TSlim = range(tis.index)

plot(tis.index[,2],tis.index[,4], ylim=TSlim, type="b", cex=.75, pch=3, main="Top 1% share", xlab="k largest values", ylab="share", col="blue")
lines(tis.index[,2],tis.index[,4], col="blue", type="l", cex=.75)
lines(tis.index[,2],tis.index[,5], col="green", type="p", cex=.75, pch=2)
lines(tis.index[,2],tis.index[,5], col="green", type="l", cex=.75)
lines(tis.index[,2],tis.index[,6], col="red", type="b", cex=.75, pch=1)
lines(tis.index[,2],tis.index[,6], col="red", type="l", cex=.75)
lines(tis.index[,2],tis.index[,7], col="gray", type="l", cex=.75)
abline(v=tis.index[(tis.index[,3]==.90),2], col="lightgray", lty=2) # 10% top obs
abline(v=tis.index[(tis.index[,3]==.95),2], col="lightgray", lty=2) #  5% top obs
abline(v=tis.index[(tis.index[,3]==.99),2], col="lightgray", lty=2) #  1% top obs

legend("topright", legend=c("Pareto 1","GPD", "EPD"), col=c("blue", "green", "red"), pch=c(3,2,1),lty=1)

legend(tis.index[(tis.index[,3]==.90),2]-top.xx,top.yy, legend=expression(italic('q')[90]), cex=.9, bty="n")
legend(tis.index[(tis.index[,3]==.95),2]-top.xx,top.yy, legend=expression(italic('q')[95]), cex=.9, bty="n")
legend(tis.index[(tis.index[,3]==.99),2]-top.xx,top.yy, legend=expression(italic('q')[99]), cex=.9, bty="n")
}
