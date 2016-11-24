#' Simulated data from the bent line quantile regression
#'
#' The function generates data from the bent line quantile regression model. 
#' 
#'  The bent line quantile regression model: 
#'  Q_tau (Y|x,z) = [1, x, (x-t_0)_+, z] * bet0. 
#'
#' @rdname bentQRsimdat
#' @param n sample size.
#' @param bet0 the vecotr of true regression coefficients.
#' @param t0 the true location of the change point.
#' @param tau the quantile level. 
#' @param modtype type of model, where
#'   modtype = 1, IID case,
#'   modetype = 2, heteroscedasticity case
#' @param etype type of error, where 
#'   etype = 1 for N(0,1),
#'   etype = 2 for t_3, 
#'   etype = 3 for 0.9N(0,1) + 0.1t_3, 
#'   etype = 4 for 0.9N(0,1) + 0.1Cauchy(0,1). 

#' @return A matrix with the elements
#' \item{y}{The response variable.}
#' \item{x}{The scalar covariate with a change point}
#' \item{z}{A vector of covariates.}
#'
#' @author Feipeng Zhang 
#' @keywords bentQRsimdat

#' @importFrom stats runif rbinom rnorm qnorm  rt qt  rcauchy qcauchy 
#'
#' @export
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
#'}
#'


bentQRsimdat <- function(n, bet0, t0, tau, modtype, etype){
  
  x <- runif(n, -2, 5)
  z <- rbinom(n, 1, 0.5)
  xz <- cbind(1, x, pmax(x-t0, 0), z)
  
  
  if(etype ==1){
    e <- rnorm(n, 0, 1) - qnorm(tau, 0, 1)    # noraml distribution
  }else if(etype ==2){
    e <- rt(n, 3) - qt(tau, 3)            # t_3 distribution
  }else if(etype ==3){
    e <- 0.9*rnorm(n, 0, 1) + 0.1*rt(n, 3) -
      (0.9*qnorm(tau, 0, 1) + 0.1*qt(tau, 3))    ## mixture 0.9*N(0,1) + 0.1 t_3
  }else if(etype ==4){
    e <- 0.9*rnorm(n, 0, 1) + 0.1*rcauchy(n, 0, 1) - 
      (0.9*qnorm(tau, 0, 1) + 0.1*qcauchy(tau, 0, 1))   ## Cauchy distribution
  }
  
  if (modtype == 1){## iid  case
    err <- e
  }  else if(modtype==2){ ## heteroscedasticity case
    err <- e*(1 + 0.3*x)  
  }
  
  
  y <- as.vector(xz %*% bet0) + err
  return(data.frame(y=y, x=x, z=z))
}

