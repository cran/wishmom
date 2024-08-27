#' Inverse of a Coefficient Matrix \eqn{\tilde{\mathcal{C}}_k}
#'
#' This function computes the inverse of the coefficient matrix \eqn{\tilde{\mathcal{C}}_k}
#'
#' @param k The order of the \eqn{\tilde{\mathcal{C}}_k} matrix
#' @param alpha The type of beta-Wishart distribution (\eqn{\alpha=2/\beta}):
#'   \itemize{
#'     \item 1/2: Quaternion Wishart
#'     \item 1: Complex Wishart
#'     \item 2: Real Wishart (default)
#'   }
#'
#' @return Inverse of a coefficient matrix \eqn{\tilde{\mathcal{C}}_k} that allows us to
#'         obtain \eqn{E[p_{\lambda}(W^{-1})W^{-r}]}, where \eqn{r+|\lambda|=k}
#'         and \eqn{W ~ W_m^{\beta}(n,\Sigma)}.  The matrix is represented as a 
#'         3-dimensional array where each slice along the third dimension represents 
#'         a coefficient matrix of the polynomial in descending powers of \eqn{\tilde{n}}.
#'
#' @export
#'
#' @examples
#' # Example 1:
#' qkn_coeff(2) # For real Wishart distribution with k = 2

#' # Example 2:
#' qkn_coeff(3, 1) # For complex Wishart distribution with k = 3
#'
#' # Example 3:
#' qkn_coeff(2, 1/2) # For quaternion Wishart distribution with k = 2
#'

qkn_coeff <- function(k, alpha = 2) {
  if (k == 1) return(1)

  c <- array(1, dim = c(1, 1, 1))  # C_1^{-1} = n1
  h <- array(1, dim = c(1, 1, 1))  # H_1^{-1} = n1
  D <- -dkmap(1, alpha)
  l2 <- dim(c)[1]

  for (i in 1:(k-1)) {
    l1 <- dim(D)[1]
    c1 <- c
    Da <- D[, 1:l2, drop=FALSE]
    Db <- D[, (l2+1):ncol(D), drop=FALSE]
    c <- array(0, dim = c(l1, l1, i + 1))
    c[1:l2, 1:l2, 1:i] <- c1
    c[(l2+1):l1, (l2+1):l1, 1:i] <- h
    for (j in 1:i) {
      c[, , j + 1] <- c[, , j + 1] + cbind(Da %*% c1[, , j], Db %*% h[, , j])
    }
    if (i < (k-1)) {
      l2 <- dim(c)[1]
      D <- -dkmap(i + 1, alpha)
      D12 <- D[1:l2, (l2+1):ncol(D)]
      D21 <- D[(l2+1):nrow(D), 1:l2]
      h <- array(0, dim = c(dim(D21)[1], dim(D12)[2], i + 1))
      for (j in 1:(i + 1)) {
        h[, , j] <- D21 %*% c[, , j] %*% D12
      }
      h <- h / (i + 1) / alpha
    }
  }

  return(c)
}
