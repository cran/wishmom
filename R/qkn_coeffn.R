#' Inverse of a Coefficient Matrix \eqn{\tilde{\mathcal{C}}_k}
#'
#' This function computes the inverse of the coefficient matrix \eqn{\tilde{\mathcal{C}}_k}
#'
#' @param k The order of the \eqn{\tilde{\mathcal{C}}_k} matrix
#' @param n1 The parameter \eqn{n - m + 1 - \alpha}, where:
#' \itemize{
#'   \item \eqn{n} is the degrees of freedom of \eqn{W}
#'   \item \eqn{m} is the number of rows of \eqn{W}
#' }
#' @param alpha The type of beta-Wishart distribution (\eqn{\alpha=2/\beta}):
#'   \itemize{
#'     \item 1/2: Quaternion Wishart
#'     \item 1: Complex Wishart
#'     \item 2: Real Wishart (default)
#'   }
#'
#' @return Inverse of a coefficient matrix \eqn{\tilde{\mathcal{C}}_k} that allows us to
#'         obtain \eqn{E[p_{\lambda}(W^{-1})W^{-r}]}, where \eqn{r+|\lambda|=k}
#'         and \eqn{W ~ W_m^{\beta}(n,\Sigma)}
#'
#' @export
#'
#' @examples
#' # Example 1:
#' qkn_coeffn(2, 20) # For real Wishart distribution with k = 2 and n1 = 20

#' # Example 2:
#' qkn_coeffn(3, 20, 1) # For complex Wishart distribution with k = 3 and n1 = 20
#'
#' # Example 3:
#' qkn_coeffn(2, 20, 1/2) # For quaternion Wishart distribution with k = 2 and n1 = 20
#'

qkn_coeffn <- function(k, n1, alpha = 2) {
  if (k == 1) return(n1)

  c <- array(n1, dim = c(1, 1))  # C_1^{-1} = n1
  h <- array(n1, dim = c(1, 1))  # H_1^{-1} = n1
  D <- -dkmap(1, alpha)
  dk2 <- 1

  for (i in 1:(k-1)) {
    l1 <- nrow(D)
    for (i1 in 1:l1) {
      D[i1, i1] = D[i1, i1] + n1
    }
    c = cbind(D[, 1:dk2] %*% c, D[, (dk2 + 1):l1] %*% h)
    if (i < k-1){
      D = -dkmap(i+1,alpha)
      dk2 = nrow(c)
      D12 = D[1:dk2, (dk2 + 1):ncol(D)]
      D21 = D[(dk2 + 1):nrow(D), 1:dk2]
      h = (D21 %*% c %*% D12)/(i + 1)/alpha
    }
  }
  return(c)
}
