#' Coefficient Matrix \eqn{\tilde{\mathcal{H}}_k}
#'
#' This function computes the coefficient matrix \eqn{\tilde{\mathcal{H}}_k} for \eqn{W \sim W_m^{\beta}(n, \Sigma)}.
#'
#' @param k The order of the \eqn{\tilde{\mathcal{H}}_k} matrix (a positive integer)
#' @param alpha The type of Wishart distribution (\eqn{\alpha = 2/\beta}):
#' \itemize{
#'   \item 1/2: Quaternion Wishart
#'   \item 1: Complex Wishart
#'   \item 2: Real Wishart (default)
#' }
#' @return A list with two elements:
#' \itemize{
#'   \item \code{c}: A 3-dimensional array containing the coefficient matrices of the numerator of \eqn{\tilde{\mathcal{H}}_k} in descending powers of \eqn{n1}, where \eqn{n1 = n - m + 1 - \alpha}
#'   \item \code{den}: A vector containing the coefficients of the denominator of \eqn{\tilde{\mathcal{H}}_k}, in descending powers of \eqn{n1}
#' }
#'
#' @examples
#' # Example 1:
#' iwish_psr(2) # For real Wishart distribution with k = 2
#'
#' # Example 2:
#' iwish_psr(4, 1) # For complex Wishart distribution with k = 4
#'
#' # Example 3:
#' iwish_psr(2, 1/2) # For quaterion Wishart distribution with k = 2
#'
#' @export

iwish_psr <- function(k, alpha = 2) {
  if (k == 1) {
    c <- array(1, dim = c(1, 1, 1))
    den <- 1
    return(list(c = c, den = den))
  }

  c <- iwish_ps(k, alpha)
  r <- dim(c)[1]
  a <- denpoly(k, alpha)
  m1 <- length(a) - k
  B <- array(0, dim = c(r, r, m1 + 1))
  B[,,1] <- diag(r)

  for (i in 1:m1) {
    for (j in 1:min(k - 1, i)) {
      B[,,i + 1] <- B[,,i + 1] - c[,,j + 1] %*% B[,,i - j + 1]
    }
    B[,,i + 1] <- B[,,i + 1] + a[i + 1] * diag(r)
  }
  return(list(c = B, den = a))
}

