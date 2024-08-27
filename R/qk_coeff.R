#' Coefficient Matrix \eqn{\mathcal{C}_k}
#'
#' This function computes the coefficient matrix \eqn{\mathcal{C}_k}, which
#' is a matrix of constants that allows us to obtain \eqn{E[p_{\lambda}(W)W^r]},
#' where \eqn{r+|\lambda|=k} and \eqn{W \sim W_m^{\beta}(n, \Sigma)}.
#'
#' @param k The order of the \eqn{\mathcal{C}_k} matrix
#' @param alpha The type of Wishart distribution (\eqn{\alpha=2/\beta}):
#'   \itemize{
#'     \item 1/2: Quaternion Wishart
#'     \item 1: Complex Wishart
#'     \item 2: Real Wishart (default)
#'   }
#'
#' @return \eqn{\mathcal{C}_k}, a matrix that allows us to obtain \eqn{E[p_{\lambda}(W)W^r]}, 
#' where \eqn{r+|\lambda|=k} and \eqn{W \sim W_m^{\beta}(n, \Sigma)}.
#' The matrix is represented as a 3-dimensional array where each slice along the third 
#' dimension represents a coefficient matrix of the polynomial in descending powers of \eqn{n}.
#'
#' @export
#'
#' @examples
#' # Example 1:
#' qk_coeff(2) # For real Wishart distribution with k = 2
#'
#' # Example 2:
#' qk_coeff(3, 1) # For complex Wishart distribution with k = 3
#'
#' # Example 3:
#' qk_coeff(2, 1/2) # For quaternion Wishart distribution with k = 2
#'
qk_coeff <- function(k, alpha = 2) {
  if (k == 1) return(1)

  # Initialize c and h
  c <- array(1, dim = c(1, 1, 1))     # C_1 = n
  h <- array(1, dim = c(1, 1, 1))     # H_1 = n
  D <- dkmap(1, alpha)
  l2 <- dim(c)[1]

  for (i in 1:(k - 1)) {
    l1 <- dim(D)[1]
    c1 <- c
    Da <- D[, 1:l2, drop=FALSE]
    Db <- D[, (l2 + 1):ncol(D), drop=FALSE]
    c <- array(0, dim = c(l1, l1, i + 1))
    c[1:l2, 1:l2, 1:i] <- c1
    c[(l2 + 1):dim(c)[1], (l2 + 1):dim(c)[2], 1:i] <- h
    for (j in 1:i) {
      c[, , j + 1] <- c[, , j + 1] + cbind(Da %*% c1[, , j], Db %*% h[, , j])
    }
    if (i < (k - 1)) {
      l2 <- dim(c)[1]
      D <- dkmap(i + 1, alpha)
      D12 <- D[1:l2, (l2 + 1):ncol(D)]
      D21 <- D[(l2 + 1):nrow(D), 1:l2]
      h <- array(0, dim = c(dim(D21)[1], dim(D12)[2], i + 1))
      for (j in 1:(i + 1)) {
        h[, , j] <- D21 %*% c[, , j] %*% D12
      }
      h <- h / ((i + 1) * alpha)
    }
  }
  return(c)
}
