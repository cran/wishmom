#' Generate Integer Partitions in Reverse Lexicographical Order
#'
#' This function generates all integer partitions of a given integer \code{k} in reverse lexicographical order.
#' The function is adapted from "Algorithm ZS1" described in Zoghbi and Stojmenovic (1998),
#' "Fast Algorithms for Generating Integer Partitions", International Journal of Computer Mathematics,
#' Volume 70, Issue 2, pages 319-332.
#'
#' @param k An integer to be partitioned
#'
#' @return A matrix where each row represents an integer partition of \code{k}.
#' The partitions are listed in reverse lexicographical order.
#'
#' @references
#' Zoghbi, A., & Stojmenovic, I. (1998). Fast Algorithms for Generating Integer Partitions.
#' International Journal of Computer Mathematics, 70(2), 319-332. DOI: 10.1080/00207169808804755
#'
#' @examples
#' # Example 1:
#' ip_desc(3)
#'
#' # Example 2:
#' ip_desc(5)
#'
#' @export

ip_desc <- function(k) {
  # Adapted from "Algorithm ZS1" in
  # ANTOINE ZOGHBI and IVAN STOJMENOVIC (1998), "FAST ALGORITHMS FOR
  # GENERATING INTEGER PARTITIONS", International Journal of Computer
  # Mathematics, Volume 70, Issue 2, pages 319-332
  # DOI: 10.1080/00207169808804755

  x <- rep(1, k)
  x[1] <- k
  m <- 1
  h <- 1
  S <- matrix(c(k, rep(0, k - m)), nrow = 1, byrow = TRUE)

  while (x[1] != 1) {
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
    S <- rbind(S, c(x[1:m], rep(0, k - m)))
  }

  return(S)
}
