#' Coefficient Matrix \eqn{\mathcal{C}_k}
#'
#' This function computes the coefficient matrix \eqn{\mathcal{C}_k}, which
#' is a matrix of constants that allows us to obtain \eqn{E[p_{\lambda}(W)W^r]},
#' where \eqn{r+|\lambda|=k} and \eqn{W \sim W_m^{\beta}(n, \Sigma)}.
#'
#' @param k The order of the \eqn{\mathcal{C}_k} matrix
#' @param n The degrees of freedom of the beta-Wishart matrix \eqn{W}
#' @param alpha The type of Wishart distribution (\eqn{\alpha=2/\beta}):
#'   \itemize{
#'     \item 1/2: Quaternion Wishart
#'     \item 1: Complex Wishart
#'     \item 2: Real Wishart (default)
#'   }
#'
#' @return \eqn{\mathcal{C}_k}, a matrix that allows us to obtain
#' \eqn{E[p_{\lambda}(W)W^r]}, where \eqn{r+|\lambda|=k} and \eqn{W \sim W_m^{\beta}(n, \Sigma)}.
#'
#' @export
#'
#' @examples
#' # Example 1:
#' qk_coeffn(2, 2) # For real Wishart distribution with k = 2 and n = 2
#'
#' # Example 2:
#' qk_coeffn(3, 2, 1) # For complex Wishart distribution with k = 3 and n = 2
#'
#' # Example 3:
#' qk_coeffn(2, 2, 1/2) # For quaternion Wishart distribution with k = 2 and n = 2

qk_coeffn <- function(k, n, alpha = 2) {
  if (k == 1) return(n)

  # Initialize c and h
  c <- array(n, dim = c(1, 1))     # C_1 = n
  h <- array(n, dim = c(1, 1))     # H_1 = n
  D <- dkmap(1, alpha)
  dk2 <- 1

  for (i in 1:(k - 1)) {
    l1 <- nrow(D)
    for (i1 in 1:l1) {
      D[i1, i1] = D[i1, i1] + n
    }
    c = cbind(D[, 1:dk2] %*% c, D[, (dk2 + 1):l1] %*% h)
    if (i < k-1){
      D = dkmap(i+1,alpha)
      dk2 = nrow(c)
      D12 = D[1:dk2, (dk2 + 1):ncol(D)]
      D21 = D[(dk2 + 1):nrow(D), 1:dk2]
      h = (D21 %*% c %*% D12)/(i + 1)/alpha
    }
  }
  return(c)
}
