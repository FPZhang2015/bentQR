
#' @title Segmented estimate for bent line quantile regression
#'
#' @description  This function estimates the bent line quantile regression model
#' by segmented estimation.
#'
#' @rdname bentQRseg
#' @param y A vector of response
#' @param x A numeric variable with change point
#' @param z A vector of covariates
#' @param tau quantile level
#' @param bet.ini  A initial vector of regression coefficients
#' @param bp.ini  A initial value of change point
#' @param modtype type of error, modtype =1 for iid, modtype = 2 for nid
#' @param tol  tolerance value, 1e-4 for default
#' @param max.iter the maximum iteration steps, 100 for default
#'
#' @return A list with the elements
#' \item{bet.est}{The estimated regression coefficients with intercept.}
#' \item{bp.est}{The estimated change point.}
#' \item{bet.se}{The estimated standard error of the regression coefficients.}
#' \item{bp.se}{The estimated standard error of the change point.}
#' \item{iter}{The iteration steps.}
#'
#' @author Feipeng Zhang
#' @keywords bentQRseg
#'
#' @importFrom  quantreg rq
#' @export
#'
#' @seealso \code{bentQRgrid} 
#' @examples
#'
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
#' fit.grid <- bentQRgrid(y, x, z, tau)
#' bet.ini <- fit.grid$bet.est
#' bp.ini <- fit.grid$bp.est
#' bentQRseg(y, x, z, tau, bet.ini, bp.ini, modtype = 1)
#'
#' ## MRS data
#' data(data_mrs)
#' x <- log(data_mrs$mass)
#' y <- log(data_mrs$speed)
#' z <- data_mrs$hopper
#' tau <- 0.5
#'
#' fit.grid <- bentQRgrid(y, x, z, tau)
#' bet.ini <- fit.grid$bet.est
#' bp.ini <- fit.grid$bp.est
#' bentQRseg(y, x, z, tau, bet.ini, bp.ini, modtype = 1)

#' }


bentQRseg <- function(y, x, z, tau,
                      bet.ini, bp.ini, modtype, tol=1e-4, max.iter=100)
{

  p <- ifelse(is.null(z), 0, ncol(as.matrix(z)))

  bet0 = c(bet.ini[c(1:3)], 0.01, bet.ini[-c(1:3)])
  t0 = bp.ini

  iter = 1
  while(iter <= max.iter) {

    if(p == 0){ # without z
      datnew<- data.frame(
        y=y, x1=x, x2=pmax(x-t0,0), x3= (-1)*ifelse(x > t0, 1, 0)
      )
    }else{
      datnew<- data.frame(
        y=y, x1=x, x2=pmax(x-t0,0), x3= (-1)*ifelse(x > t0, 1, 0), z
      )
      names(datnew) = c("y", "x1", "x2", "x3", paste0("z", 1:p, seq = ""))
    }


    if(modtype==1){ ## modtype =1, iid
      fit1 <- rq(y~., tau = tau, data = datnew, method ="br")
      bet1 <- fit1$coef
      bet1.se <- summary(fit1, se ="iid")$coef[, 2]
    } else if (modtype==2){ ## modtype =2, nid
      fit1 <- rq(y~., tau = tau, data = datnew, method ="br")
      bet1 <- fit1$coef
      bet1.se <- summary(fit1, se ="nid")$coef[, 2]
    }

    names(bet1) <- names(bet1.se) <- NULL
    t1 <-  t0 + 0.7*bet1[4]/(bet1[3]+1e-10)
    t1.se <- bet1.se[4]/(abs(bet1[3])+1e-10)

    if(abs(bet1[4]) < tol && max(abs(bet1 - bet0))<tol) break

    # check if psi1 is admissible
    a<- (max(x) <= t1)
    b<- (min(x) >= t1)
    isErr<- sum(a+b)!=0 || is.na(sum(a+b))
    if(isErr){stop("estimated tau out of its range")}

    bet0 = bet1
    t0 = t1
    iter = iter+1  # increment index
  } # end while

  bet.est <- bet1[-4]
  bp.est <- t1
  bet.se <- bet1.se[-4]
  bp.se <- t1.se

  return(list(bet.est = bet.est, bp.est = bp.est,
              bet.se = bet.se, bp.se = t1.se,
              iter = iter))
}
