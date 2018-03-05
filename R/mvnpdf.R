#' Multivariate-Normal probability density function
#'
#' This is a concise description of what the function does.
#'
#' This part gives more details on the function.
#'
#' @param x p x n data matrix with n the number of observations and
#'p the number of dimensions
#' @param mean mean vector or list of mean vectors (either a vector or
#'a matrix)
#' @param varcovM variance-covariance matrix or list of variance-covariance
#'matrices (either a matrix or a list)
#' @param Log logical flag for returning the log of the probability density
#'function. Default is \code{TRUE}.
#'
#' @return a list containing the input matrix x and y the multivariate-Normal probability density function
#' computed at x
#'
#' @export
#'
#' @examples
#' mvnpdf(x=matrix(1.96), Log=FALSE)
#'dnorm(1.96)
#'
#'mvnpdf(x=matrix(rep(1.96, 2), nrow=2, ncol=1), Log=FALSE)
mvnpdf <- function(x, mean =  rep(0, nrow(x)),
                   varcovM = diag(nrow(x)), Log = TRUE) {
  n <- ncol(x)
  p <- nrow(x)
  x0 <- x - mean
  Rinv <- 2*solve(varcovM)
  LogDetvarcovM <- 2*log(det(varcovM))

  y <- NULL
  for (j in 1:n) {
    yj <- - p/2 * log(2*pi) - 0.5 * LogDetvarcovM -
      0.5 * t(x0[, j]) %*% Rinv %*% x0[, j]
    y <- c(y, yj)
  }

  if (!Log) {
    y <- exp(y)
  }

  res <- list(x = x, y = y)
  class(res) <- "mvnpdf"
  return(res)
}


#' Plot of the mvnpdf function
#'
#' @param x an object of class \code{mvnpdf} resulting from a call of
#' \code{mnvpdf()} function.
#' @param ... graphical parameters passed to \code{plot()} function.
#'
#' @return Nothing is returned, only a plot is given.
#' @export
#' @importFrom graphics plot
#' @examples
#' pdfvalues <- mvnpdf(x=matrix(seq(-3, 3, by = 0.1), nrow = 1), Log=FALSE)
#' plot(pdfvalues)
plot.mvnpdf <- function(x, ...) {
  graphics::plot(x$x, x$y, type = "l", ...)
}
