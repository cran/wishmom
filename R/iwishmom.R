#' Expectation of a Matrix-valued Function of an Inverse beta-Wishart Distribution
#'
#' When \code{iw = 0}, the function calculates \eqn{E[\prod_{j=1}^r \mbox{tr}(W^{-j})^{f_j}]},
#' where \eqn{W \sim W_m^{\beta}(n, S)}.  When \code{iw != 0},
#' the function calculates \eqn{E[\prod_{j=1}^r \mbox{tr}(W^{-j})^{f_j}W^{-iw}]}.
#'
#' @param n The degrees of freedom of the beta-Wishart matrix \eqn{W}
#' @param S The covariance matrix of the beta-Wishart matrix \eqn{W}
#' @param f A vector of nonnegative integers \eqn{f_j} that represents
#'          the power of \eqn{\mbox{tr}(W^{-j})}, where \eqn{j=1, \ldots, r}
#' @param iw The power of the inverse beta-Wishart matrix \eqn{W^{-1}} (0 by default)
#' @param alpha The type of Wishart distribution \eqn{(\alpha=2/\beta)}:
#'   \itemize{
#'     \item 1/2: Quaternion Wishart
#'     \item 1: Complex Wishart
#'     \item 2: Real Wishart (default)
#'   }
#'
#' @return When \code{iw = 0}, it returns \eqn{E[\prod_{j=1}^r \mbox{tr}(W^{-j})^{f_j}]}.
#' When \code{iw != 0}, it returns \eqn{E[\prod_{j=1}^r \mbox{tr}(W^{-j})^{f_j}W^{-iw}]}.
#'
#' @export
#'
#' @examples
#' # Example 1: For E[tr(W^{-1})^2] with W ~ W_m^1(n,S),
#' # where n and S are defined below:
#' n <- 20
#' S <- matrix(c(25, 49,
#'               49, 109), nrow=2, ncol=2)
#' iwishmom(n, S, 2) # iw = 0, for real Wishart distribution
#'
#' # Example 2: For E[tr(W^{-1})^2*tr(W^{-3})W^{-2}] with W ~ W_m^1(n,S),
#' # where n and S are defined below:
#' n <- 20
#' S <- matrix(c(25, 49,
#'               49, 109), nrow=2, ncol=2)
#' iwishmom(n, S, c(2, 0, 1), 2, 2) # iw = 2, for real Wishart distribution
#'
#' # Example 3: For E[tr(W^{-1})^2*tr(W^{-3})] with W ~ W_m^2(n,S),
#' # where n and S are defined below:
#' # Hermitian S for the complex case
#' n <- 20
#' S <- matrix(c(25, 49 + 2i,
#'               49 - 2i, 109), nrow=2, ncol=2)
#' iwishmom(n, S, c(2, 0, 1), 0, 1) # iw = 0, for complex Wishart distribution
#'
#' # Example 4: For E[tr(W^{-1})*tr(W^{-2})^2*tr(W^{-3})^2*W^{-1}] with W ~ W_m^2(n,S),
#' # where n and S are defined below:
#' n <- 30
#' S <- matrix(c(25, 49 + 2i,
#'               49 - 2i, 109), nrow=2, ncol=2)
#' iwishmom(n, S, c(1, 2, 2), 1, 1) # iw = 1, for complex Wishart distribution

iwishmom <- function(n, S, f, iw = 0, alpha = 2) {
  m <- nrow(S)
  n1 <- n - m + 1 - alpha
  lf <- length(f)
  k1 <- sum(f * seq_len(lf))
  k <- k1 + iw

  #  Check existence of moments
  if (n1 <= 2*k) {
    if (iw == 0) {
      Q <- NaN
    } else {
      Q <- matrix(NaN, nrow = m, ncol = m)
    }
    return(Q)
  }

  Si <- solve(S)
  diag(Si) <- Re(diag(Si))

  #  Fast return for k=1

  if (k == 1) {
    if (iw == 0) {
      Q <- sum(diag(Si))/n1
    } else {
      Q <- Si/n1
    }
    return(Q)
  }

  tS <- k
  S1 <- Si
  tS[1] <- sum(diag(S1))
  for (i in 2:k) {
    S1 <- S1 %*% Si
    tS[i] <- sum(diag(S1))
  }
  tS <- Re(tS)

  # Initializing F matrix
  F <- matrix(1, k, k)
  for (j in 2:k) {
    F[, j] <- F[, j - 1]
    F[j, j] <- F[j, j] + 1
    if (j < k) {
      for (i in (j + 1):k) {
        F[i, j] <- F[i, j] + F[i - j, j]
      }
    }
  }
  dd <- cumsum(c(1, diag(F)))

  ind <- function(f, f1 = 0) {
    if (k1 == 0) return(1)
    ind1 <- F[k1, k1]
    ss <- k1
    if (length(f) >= 2) {
      for (jj in length(f):2) {
        if (1 <= f[jj]) {
          for (kk in 1:f[jj]) {
            ind1 <- ind1 - F[ss, jj - 1]
            ss <- ss - jj
          }
        }
      }
    }
    if (f1 > 0) ind1 <- ind1 + dd[k1]
    return(ind1)
  }

  if (iw == 0) {  # The output is a scalar for this case
    if (k == 0) return(m)
    c <- iwish_ps(k, alpha)
    temp <- matrix(0, nrow = dim(c)[1], ncol = dim(c)[1])
    for (i in 1:k) {
      temp <- temp + c[,,i]*n1^(k-i+1)
    }
    c <- solve(temp)
    c <- c[ind(f), ]
    pp <- ip_desc(k)
    Q <- c[1] * tS[k]
    for (i in 2:nrow(pp)) {
      pp1 <- pp[i, ]
      pp1 <- pp1[pp1 > 0]
      Q <- Q + c[i] * prod(tS[pp1])
    }
  } else {
    c <- qkn_coeff(k, alpha)
    temp <- matrix(0, nrow = dim(c)[1], ncol = dim(c)[1])
    for (i in 1:k) {
      temp <- temp + c[,,i]*n1^(k-i+1)
    }
    c <- solve(temp)
    c <- rev(c[ind(f,iw), ])
    S1 <- Si
    Q <- matrix(0, m, m)
    count <- 1
    for (i in 1:(k - 1)) {
      y <- 0
      pp <- ip_desc(k-i)
      for (j in nrow(pp):1) {
        pp1 <- pp[j, ]
        pp1 <- pp1[pp1 > 0]
        y <- y + c[count] * prod(tS[pp1])
        count <- count + 1
      }
      Q <- Q + y * S1
      S1 <- S1 %*% Si
      diag(S1) <- Re(diag(S1))
    }
    Q <- Q + c[count] * S1
  }
  return(Q)
}
