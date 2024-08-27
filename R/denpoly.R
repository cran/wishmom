#' Coefficients of the Denominator Polynomial for \eqn{\tilde{H}_k} and \eqn{\tilde{C}_k}
#'
#' This function computes the coefficients of the denominator polynomial for the elements of
#' \eqn{\tilde{H}_k} and \eqn{\tilde{C}_k}.
#' The function returns a vector containing the coefficients in descending powers of
#' \eqn{\tilde{n}}, with the last element being the coefficient of \eqn{\tilde{n}}.
#'
#' @param k The order of the polynomial (a positive integer)
#' @param alpha The type of Wishart distribution \eqn{(\alpha=2/\beta)}:
#'   \itemize{
#'     \item 1/2: Quaternion Wishart
#'     \item 1: Complex Wishart
#'     \item 2: Real Wishart (default)
#'   }
#'
#' @return A vector containing the coefficients of the denominator
#'         polynomial in descending powers of \eqn{\tilde{n}} for the elements of
#'         \eqn{\tilde{H}_k} and \eqn{\mathcal{C}_k}.
#'
#' @export
#'
#' @examples
#' # Example 1: Compute the denominator polynomial for k = 3, alpha = 2
#' # Output corresponds to the polynomial n1^5-3n1^4-8n1^3+12n1^2+16n1,
#' # where n1 is \eqn{\tilde{n}}
#' denpoly(3)
#'
#' # Example 2: Compute the denominator polynomial for k = 2, alpha = 1
#' # Output corresponds to the polynomial n1^3-n1, where n1 is \eqn{\tilde{n}}
#' denpoly(2, alpha = 1)
#'

denpoly <- function(k, alpha = 2) {
  if (k == 1) return(1)
  pp <- ip_desc(k)
  dk <- nrow(pp)

  y1 <- c((-(k-1):0) * alpha, 1:(k-1))

  for (i in 2:(dk - 1)) {
    pp1 <- pp[i, ]
    pp1 <- pp1[pp1 > 0]

    q <- NULL
    for (j in 1:length(pp1)) {
      q <- c(q, j - 1 - alpha * (0:(pp1[j] - 1)))
    }

    a1 <- unique(q)
    for (j in 1:length(a1)) {
      nn <- sum(q == a1[j]) - sum(y1 == a1[j])
      if (nn > 0) y1 <- c(y1, rep(a1[j], nn))
    }
  }

  y1 <- y1[-k]  # Get rid of the term n1 in |\tilde{H}_k^{-1}|
  y <- rep(0, length(y1) + 1)
  y[1] <- 1
  for (i in 1:length(y1)) {
    y[2:(i + 1)] <- y[2:(i + 1)] + y1[i] * y[1:i]
  }

  return(y)
}
