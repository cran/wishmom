#' Symbolic Expectation of a Matrix-valued Function of a beta-Wishart Distribution
#'
#' When \code{iw = 0}, the function returns an analytical expression of
#' \eqn{E[\prod_{j=1}^r \mbox{tr}(W^j)^{f_j}]}, where \eqn{W \sim W_m^{\beta}(n, S)}.
#' When \code{iw != 0}, the function returns an analytical expression of
#' \eqn{E[\prod_{j=1}^r \mbox{tr}(W^j)^{f_j}W^{iw}]}.
#' For a given \code{f}, \code{iw}, and \code{alpha}, this function provides the aforementioned
#' expectations in terms of the variables \eqn{n} and \eqn{\Sigma}.
#'
#' @param f A vector of nonnegative integers \eqn{f_j} that represents
#'          the power of \eqn{\mbox{tr}(W^{-j})}, where \eqn{j=1, \ldots, r}
#' @param iw The power of the beta-Wishart matrix \eqn{W^{-1}} (0 by default)
#' @param alpha The type of Wishart distribution \eqn{(\alpha=2/\beta)}:
#'   \itemize{
#'     \item 1/2: Quaternion Wishart
#'     \item 1: Complex Wishart
#'     \item 2: Real Wishart (default)
#'   }
#' @param latex A Boolean indicating whether the output will be a LaTeX string or a dataframe
#'   (FALSE by default)
#'
#' @return When \code{iw = 0}, it returns an analytical expression of
#'         \eqn{E[\prod_{j=1}^r \mbox{tr}(W^j)^{f_j}]}.
#'         When \code{iw != 0}, it returns an analytical expression of
#'         \eqn{E[\prod_{j=1}^r \mbox{tr}(W^j)^{f_j}W^{iw}]}.
#'         If \code{latex = FALSE}, the output is a data frame that stores the coefficients
#'         for calculating the result. If \code{latex = TRUE}, the output is a LaTeX
#'         formatted string of the result in terms of \eqn{n} and \eqn{\Sigma}.
#'
#' @export
#'
#' @examples
#' # Example 1: For E[tr(W)^4] with W ~ W_m^1(n,Sigma), represented as a dataframe:
#' wishmom_sym(4) # iw = 0, for real Wishart distribution
#'
#' # Example 2: For E[tr(W)*tr(W^2)W] with W ~ W_m^1(n,S), represented as a dataframe:
#' wishmom_sym(c(1, 1), 1) # iw = 1, for real Wishart distribution
#'
#' # Example 3: For E[tr(W)^4] with W ~ W_m^2(n,S), represented as a LaTeX string:
#' # Using writeLines() to format
#' writeLines(wishmom_sym(4, 0, 1, latex=TRUE)) # iw = 0, for complex Wishart distribution
#'
#' # Example 4: For E[tr(W)*tr(W^2)W] with W ~ W_m^2(n,S), represented as a LaTeX string:
#' # Using writeLines() to format
#' writeLines(wishmom_sym(c(1, 1), 1, 1, latex=TRUE)) # iw = 1, for real Wishart distribution

wishmom_sym <- function(f, iw = 0, alpha = 2, latex=FALSE) {
  lf <- length(f)
  k1 <- sum(f * seq_len(lf))
  k <- k1 + iw

  #  Fast return for k=1

  if (k == 1) {
    if (latex) {
      if (iw == 0) {
        return("n p_{(1)}(\\Sigma)");
      } else {
        return("n\\Sigma")
      }
    }
    if (iw ==0) {
      results <- data.frame(
        kappa = "1",
        h_kappa = "n"
      )
    } else {
        results <- data.frame(
          i = 1,
          rho = c(NA),
          c = "n"
        )
    }
    return(results)
  }

  #  F(i,j) is the number of integer partitions of i with j parts or less,
  #  it is also the number of integer partitions of i whose first part
  #  is smaller than or equal to j.

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

#
#  This function produces a string expression that represents
#  the power of i.  When i=1, the function returns a null string.
#  When i>=10, the function returns a pair of braces around
#  the power.
#
  powers <- function(i)
    if (i >= 10) {
      return(paste0("^{", i, "}"))
    } else if (i==1) {
      return("")
    } else {
      return(paste0("^", i))
    }

#
#  This function produces the string expression of a
#  polynomial of n, with the coefficients of declining
#  powers of n stored in v, with the last element being the
#  constant term.
#
  n_poly <- function(v) {
    expr <- ""
    lv = length(v)
    for (i in 1:lv) {
      coeff <- v[i]
      if (coeff != 0) {
        if (expr != "" && coeff>0) expr <- paste0(expr, "+")
        if (i == lv) {
          expr <- paste0(expr, coeff)
        } else {
          if (coeff == 1) {
            expr <- paste0(expr, "n", powers(lv-i))
          } else if (coeff == -1) {
            expr <- paste0(expr, "-n", powers(lv-i))
          } else {
            expr <- paste0(expr, coeff, "n", powers(lv-i))
          }
        }
      }
    }
    return(expr)
  }

  if (iw == 0) {
    c <- wish_ps(k, alpha)
    c <- c[ind(f), ,]
    partitions <- ip_desc(k)
    kappa_lst <- c()
    coeff_lst <- c()
    for (i in 1:nrow(c)) {
      partition <- partitions[i,]
      kappa_lst <- c(kappa_lst, paste0(partition[partition !=0 ], collapse=","))
      coeff_lst <- c(coeff_lst, n_poly(c[i, ]))
    }

    results <- data.frame(
      kappa = kappa_lst,
      h_kappa = coeff_lst
    )
    if (!latex) return(results)

    latex_expr <- paste0("(", results[1, "h_kappa"], ")p_{(", results[1, 1], ")}(\\Sigma)")
    if (nrow(results) == 1) return(latex_expr)
    for (i in 2:nrow(results)) {
      latex_expr <- paste0(latex_expr, "\n+(", results[i, "h_kappa"], ")p_{(", results[i, 1], ")}(\\Sigma)")
    }
  } else {
    c <- qk_coeff(k, alpha)
    c <- c[ind(f,iw), , ]
    c <- c[nrow(c):1, ]
    count <- 1

    i_lst <- c()
    rho_lst <- c()
    coeff_lst <- c()

    for (i in 1:(k - 1)) {
      pp2 <- ip_desc(k-i)
      np2 = nrow(pp2)
      if (np2 == 1) {
        i_lst <- c(i_lst, i)
        rho_lst <- c(rho_lst, "1")
        coeff_lst <- c(coeff_lst, toString(n_poly(c[count,])))
      } else {
        for (j in 1:np2) {
          coeff <- n_poly(c[count+np2-j, ])
          i_lst <- c(i_lst, i)
          partition <- pp2[j, ]
          rho_lst <- c(rho_lst, paste0(partition[partition !=0 ], collapse=","))
          coeff_lst <- c(coeff_lst, toString(coeff))
        }
      }
      count <- count+np2
    }
    coeff <- n_poly(c[count, ])
    i_lst <- c(i_lst, k)
    rho_lst <- c(rho_lst, NA)
    coeff_lst <- c(coeff_lst, toString(coeff))

    results <- data.frame(
      i = i_lst,
      rho = rho_lst,
      c = coeff_lst
    )

    if (!latex) return(results)

    formatted_list <- list()
    i <- 1
    group_string <- ""
    for (index in 1:nrow(results)) {
      row <- results[index, ]
      if (row["i"] > i) {
        if (i == (k-1)) {
          group_string <- paste0(group_string, "\\Sigma", powers(k-1), "\n")
        } else {
          group_string <- paste0("[", group_string, "]\\Sigma", powers(i), "\n")
        }
        formatted_list <- c(formatted_list, group_string)
        i <- i + 1
        group_string <- ""
      }
      if (i == k){
        group_string <- paste0("(", row["c"], ")\\Sigma", powers(k))
        formatted_list <- c(formatted_list, group_string)
      } else {
        if (nchar(group_string) > 0) group_string <- paste0(group_string, "+")
        group_string <- paste0(group_string, "(", row["c"], ")p_{(", row["rho"], ")}(\\Sigma)")
      }
    }
    latex_expr <- paste(formatted_list, collapse = "+")
  }
  return(latex_expr)
}
