#' Coefficient Matrix \eqn{\mathcal{H}_k}
#'
#' This function computes the coefficient matrix \eqn{\mathcal{H}_k} that allows us to compute
#' the expected value of a power-sum symmetric function of \eqn{W}, where
#' \eqn{W \sim W_m^{\beta}(n,\Sigma)}.
#'
#' @param k The order of the \eqn{\mathcal{H}_k} matrix
#' @param n The degrees of freedom of \eqn{W}
#' @param alpha The type of Wishart distribution (\eqn{\alpha = 2/\beta}):
#'   \itemize{
#'     \item 1/2: Quaternion Wishart
#'     \item 1: Complex Wishart
#'     \item 2: Real Wishart (default)
#'   }
#'
#' @return A coefficient matrix \eqn{\mathcal{H}_k} that allows us to compute
#'         the expected value of a power-sum symmetric function of \eqn{W},
#'         where \eqn{W \sim W_m^{\beta}(n,\Sigma)}.
#'
#' @export
#'
#' @examples
#' # Example 1:
#' wish_psn(3, 10) # For real Wishart distribution with k = 3 and n = 10
#'
#' # Example 2:
#' wish_psn(4, 10, 1) # For complex Wishart distribution with k = 4 and n = 10
#'
#' # Example 3:
#' wish_psn(2, 10, 1/2) # For quaternion Wishart distribution with k = 2 and n = 10

wish_psn <- function(k, n, alpha = 2) {
  if (k == 1) return(n)
  c <- array(n, dim = c(1, 1))  # C_1 = n
  h <- array(n, dim = c(1, 1))  # H_1 = n
  D <- dkmap(1, alpha)
  dk2 <- 1
  for (i in 1:(k - 1)) {
    l1 = nrow(D)
    for (i1 in 1:l1) {
      D[i1, i1] = D[i1, i1] + n
    }
    c = cbind(D[, 1:dk2] %*% c, D[, (dk2 + 1):ncol(D)] %*% h)
    D = dkmap(i + 1,alpha)
    dk2 <- nrow(c)
    D12 = D[1:dk2, (dk2 + 1):ncol(D)]
    D21 = D[(dk2 + 1):nrow(D), 1:dk2]
    h = (D21 %*% c %*% D12)/(i + 1)/alpha
  }
  return(h)
}
