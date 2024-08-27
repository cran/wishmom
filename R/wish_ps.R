#' Coefficient Matrix \eqn{\mathcal{H}_k}
#'
#' This function computes the coefficient matrix \eqn{\mathcal{H}_k} that allows us to compute
#' the expected value of a power-sum symmetric function of \eqn{W}, where
#' \eqn{W \sim W_m^{\beta}(n,\Sigma)}.
#'
#' @param k The order of the \eqn{\mathcal{H}_k} matrix
#' @param alpha The type of Wishart distribution (\eqn{\alpha = 2/\beta}):
#'   \itemize{
#'     \item 1/2: Quaternion Wishart
#'     \item 1: Complex Wishart
#'     \item 2: Real Wishart (default)
#'   }
#'
#' @return A coefficient matrix \eqn{\mathcal{H}_k} that allows us to compute
#'         the expected value of a power-sum symmetric function of \eqn{W},
#'         where \eqn{W \sim W_m^{\beta}(n,\Sigma)}.  The matrix is represented as a 
#'         3-dimensional array where each slice along the third dimension represents 
#'         a coefficient matrix of the polynomial in descending powers of \eqn{n}.
#'
#' @export
#'
#' @examples
#' # Example 1:
#' wish_ps(3) # For real Wishart distribution with k = 3
#'
#' # Example 2:
#' wish_ps(4, 1) # For complex Wishart distribution with k = 4
#'
#' # Example 3:
#' wish_ps(2, 1/2) # For quaternion Wishart distribution with k = 2
#'
wish_ps <- function(k, alpha = 2) {
  if (k == 1) return(1)

  c <- array(1, dim = c(1, 1, 1))  # C_1 = n
  h <- array(1, dim = c(1, 1, 1))  # H_1 = n
  D <- dkmap(1, alpha)
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
      c[, , (j + 1)] <- c[, , (j + 1)] + cbind(Da %*% c1[, , j], Db %*% h[, , j])
    }
    l2 <- dim(c)[1]
    D <- dkmap(i + 1, alpha)
    D12 <- D[1:l2, (l2+1):ncol(D)]
    D21 <- D[(l2+1):nrow(D), 1:l2]
    h <- array(0, dim = c(dim(D21)[1], dim(D12)[2], i + 1))
    for (j in 1:(i + 1)) {
      h[, , j] <- D21 %*% c[, , j] %*% D12
    }
    h <- h / (i + 1) / alpha
  }

  return(h)
}
