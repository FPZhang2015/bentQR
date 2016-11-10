
#' @title   Grid search algorithm for bent line quantile regression
#' @description This function estimates the bent line regression model
#' by grid search algorithm. 
#'
#'
#'
#' @rdname bentQRgrid
#' @param y A vector of response
#' @param x A numeric variable with change point
#' @param z A vector of covariates
#' @param tau quantile level

#'
#' @return A list with the elements
#' \item{bet.est}{The estimated regression coefficients with intercept.}
#' \item{bp.est}{The estimated change point.}
#' \item{bet.se}{The estimated standard error of the regression coefficients.}
#' \item{bp.se}{The estimated standard error of the change point.}
#'
#' @author Feipeng Zhang
#' @keywords bentQRgrid
#'
#' @importFrom  quantreg rq
#' @importFrom stats sd quantile
#' @export
#' @references Li, C., Wei, Y., Chappell, R., and He, X. (2011). 
#' Bent line quantile regression with application to
#' an allometric study of land mammals' speed and mass.Biometrics, 67:242--249.
#' 
#'
#' @examples
#' \dontrun{
#' 
#' ## simulated data
#' n <- 200
#' t0 <- 1.5    # true change point
#' bet0 <- c(1, 3, -2, 1)     # true regression coefficients
#' tau <- 0.7
#' modtype <- 1
#' etype <- 1
#' dat <- bentQRsimdat(n, bet0, t0, tau, modtype, etype)
#' head(dat) 
#' y <- dat$y
#' x <- dat$x
#' z <- dat$z
#' 
#' bentQRgrid(y, x, z, tau)
#'
#' ## MRS data
#' data(data_mrs)
#' x <- log(data_mrs$mass)
#' y <- log(data_mrs$speed)
#' z <- data_mrs$hopper
#' tau <- 0.5
#' bentQRgrid(y, x, z, tau)

#' }



bentQRgrid <- function(y, x, z, tau){

  #### estimate
  ## check function

  checkfun <- function(u){
    u*(tau-ifelse(u<0, 1, 0))
  }

  qnregprof <- function(y, x, z, tau)
  {

    tt <- seq(quantile(x,0.01), quantile(x,0.99), length=100)   # grid search interval

    p <- ifelse(is.null(z), 0, ncol(as.matrix(z)))
    sst<- NULL
    bet <- matrix(0, length(tt), p+3)
    for (kk in 1:length(tt)){

      if(p == 0){ # without z
        datnew<- data.frame(
          y=y, x1=x, x2=pmax(x-tt[kk],0)
        )
      }else{
        datnew<- data.frame(
          y=y, x1=x, x2=pmax(x-tt[kk],0), z
        )
        names(datnew) = c("y", "x1", "x2", paste0("z", 1:p, seq = ""))
      }

      fit <- rq(y~., tau=tau, data= datnew, method="br")
      bet[kk, ] <- fit$coef
      sst[kk] <- sum(checkfun(fit$residuals))
    }
    bp.grid <- min(tt[sst== min(sst)])
    bet.grid <- bet[tt==bp.grid, ]

    return(list(bet= bet.grid, bp=bp.grid))
  }


  # Epanechnikov kernel
  Ku = function(u)
  {
    3/4*(abs(u)<=1)*(1-u^2)
  }



  #### standard error estimate
  qnregEse <- function(y, x, z, bet.est, t.est)
  {
    n <- length(y)
    beta2 <- bet.est[3]

    gest <- as.vector(cbind(1, x, pmax(x-t.est, 0), z) %*% bet.est)
    res <- y - gest
    hn <- 1.06* n^(-1/5) * sd(res)  ## the bandwidth by Silverman's rule of thumb
    ker <- Ku(res/hn)/hn

    gdev <- cbind(1, x, pmax(x-t.est, 0), z, (-1)*beta2*ifelse(x>t.est, 1, 0))
    Jn <- tau*(1-tau)/n * t(gdev)%*% gdev
    Ln <- n^(-1) * t(gdev) %*% diag(ker) %*% gdev

    Ln.inv <- solve(Ln+1e-8)
    Sig_thet <- n^(-1) * Ln.inv %*% Jn %*% Ln.inv
    ese <- as.vector(sqrt(diag(Sig_thet)))
    names(ese) <- NULL
    return(ese)
  }


  fit.grid <- qnregprof(y, x, z, tau)
  bet.est <- fit.grid$bet
  bp.est <- fit.grid$bp
  thet.se <- qnregEse(y, x, z, fit.grid$bet, fit.grid$bp)
  bet.se <- thet.se[-length(thet.se)]
  bp.se <- thet.se[length(thet.se)]

  return(list(bet.est = bet.est, bp.est = bp.est,
              bet.se = bet.se, bp.se = bp.se))

}
