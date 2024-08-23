#' Mapping Matrix that maps \eqn{Q_{k+1}} to \eqn{Q_k} for a beta-Wishart Distribution, but without \eqn{n} on the diagonal
#'
#' This function computes the matrix that maps \eqn{Q_{k+1}} to \eqn{Q_k} when \eqn{W \sim W_m^{\beta}(n, \Sigma)}.
#'
#' @param k The order of the mapping matrix \eqn{D_k} (a positive integer)
#' @param alpha The type of beta-Wishart distribution (\eqn{\alpha=2/\beta}):
#'   \itemize{
#'     \item 1/2: Quaternion Wishart
#'     \item 1: Complex Wishart
#'     \item 2: Real Wishart (default)
#'   }
#'
#' @return A matrix that maps \eqn{Q_{k+1}} to \eqn{Q_k}, but without \eqn{n} on the diagonal.
#'
#' @export
#'
#' @examples
#' # Example 1: Compute the mapping matrix for k = 2, real Wishart
#' dkmap(2)

#' # Example 2: Compute the mapping matrix for k = 1, complex Wishart
#' dkmap(1, 1)

#' # Example 3: Compute the mapping matrix for k = 2, quaternion Wishart
#' dkmap(2, 1/2)

dkmap <- function(k, alpha = 2) {

  # F(i, j) is the number of integer partitions of i with j parts or less,
  # it is also the number of integer partitions of i whose first part
  # is smaller than or equal to j.
  if (k==1) {
    D <- matrix(1, nrow=2, ncol=2)
    D[1, 1] <- alpha - 1
    D[2, 1] <- alpha
    D[2, 2] <- 0
    return(D)
  }

  F <- matrix(1, nrow = k, ncol = k)
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
  lp2 <- F[k, k]
  dk1 <- dd[k + 1]
  dk2 <- dd[k]
  D <- matrix(0, nrow = dk1, ncol = dk1)

  ind <- function(p, f1) {
    # This function takes an integer partition p and an integer f1
    # and returns the location of W^f1*prod_{i=1}^l(p)tr(W^p_i)
    # within the list of Q_{k-1}, where k-1=f1+|p|.
    k1 <- sum(p)
    ind1 <- F[k1, k1]
    if (!missing(f1) && f1 > 0) {
      ind1 <- ind1 + dd[k1]
    }
    if (length(p) >= 1) {
      for (i2 in 1:length(p)) {
        kk <- p[i2]
        if (kk > 1) {
          ind1 <- ind1 - F[k1, kk - 1]
          k1 <- k1 - kk
        } else {
          break
        }
      }
    }
    return(ind1)
  }

  # Recursion for E[W^{k+1}]
  D[1, 1] <- (alpha - 1) * k
  D[1, dd[1:k] + 1] <- 1
  #  Recursion for E[W^k*tr(W)]
  D[2,1] = alpha
  D[2,2] = (alpha-1)*(k-1)
  count = 2;

  for (i in 2:k) {
    # These are terms for W^{k+1-i}
    count <- count + 1
    D[count, 1] <- alpha * i
    D[count, count] <- (alpha - 1) * (k - i)
    x <- c(i, matrix(1, nrow = 1, ncol = i-1))
    m <- 1
    h <- 1
    for (ii in 2:F[i, i]) {
      if (x[h] == 2) {
        m <- m + 1
        x[h] <- 1
        h <- h - 1
      } else {
        r <- x[h] - 1
        t <- m - h + 1
        x[h] <- r
        while (t >= r) {
          h <- h + 1
          x[h] <- r
          t <- t - r
        }
        if (t == 0) {
          m <- h
        } else {
          m <- h + 1
          if (t > 1) {
            h <- h + 1
            x[h] <- t
          }
        }
      }

      count <- count + 1
      # Construct D_{k,11}
      if (count <= dk2) {
        D[count, count] <- (alpha - 1) * (k - i)
      }
      # Construct D_{k,21} and D_{k,12}
      pp <- x[1:m]
      for (j in 1:x[1]) {
        v <- which(pp == j)
        lj <- length(v)
        if (lj > 0) {
          pp1 <- pp
          pp1 <- pp1[-v[1]]
          i1 <- ind(pp1, k - i + j)
          D[count, i1] <- alpha * j * lj
          D[i1, count] <- 1
        }
      }
    }
  }
  return(D)
}

