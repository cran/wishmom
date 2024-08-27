#' Symbolic Expectation of a Matrix-valued Function of an Inverse beta-Wishart Distribution
#'
#' When \code{iw = 0}, the function returns an analytical expression of
#' \eqn{E[\prod_{j=1}^r \mbox{tr}(W^{-j})^{f_j}]}, where \eqn{W \sim W_m^{\beta}(n, S)}.
#' When \code{iw != 0}, the function returns an analytical expression of
#' \eqn{E[\prod_{j=1}^r \mbox{tr}(W^{-j})^{f_j}W^{-iw}]}.
#' For a given \code{f}, \code{iw}, and \code{alpha}, this function provides the aforementioned
#' expectations in terms of the variables \eqn{\tilde{n}} and \eqn{\Sigma}.
#'
#' @param f A vector of nonnegative integers \eqn{f_j} that represents
#'          the power of \eqn{\mbox{tr}(W^{-j})}, where \eqn{j=1, \ldots, r}
#' @param iw The power of the inverse beta-Wishart matrix \eqn{W^{-1}} (0 by default)
#' @param alpha The type of Wishart distribution \eqn{(\alpha=2/\beta)}:
#'   \itemize{
#'     \item 1/2: Quaternion Wishart
#'     \item 1: Complex Wishart
#'     \item 2: Real Wishart (default)
#'   }
#' @param latex A Boolean indicating whether the output will be a LaTeX string or dataframe
#'   (FALSE by default)
#'
#' @return When \code{iw = 0}, it returns an analytical expression of
#'         \eqn{E[\prod_{j=1}^r \mbox{tr}(W^{-j})^{f_j}]}.
#'         When \code{iw != 0}, it returns an analytical expression of
#'         \eqn{E[\prod_{j=1}^r \mbox{tr}(W^{-j})^{f_j}W^{-iw}]}.
#'         If \code{latex = FALSE}, the output is a data frame that stores the
#'         coefficients for calculating the result. If \code{latex = TRUE}, the
#'         output is a LaTeX formatted string of the result in terms of
#'         \eqn{\tilde{n}} and \eqn{\Sigma}.
#'
#' @export
#'
#' @examples
#' # Example 1: For E[tr(W^{-1})^4] with W ~ W_m^1(n,Sigma), represented as a dataframe:
#' iwishmom_sym(4) # iw = 0, for real Wishart distribution
#'
#' # Example 2: For E[tr(W^{-1})*tr(W^{-2})W^{-1}] with W ~ W_m^1(n,S), represented as a dataframe:
#' iwishmom_sym(c(1, 1), 1) # iw = 1, for real Wishart distribution
#'
#' # Example 3: For E[tr(W^{-1})^4] with W ~ W_m^2(n,S), represented as a LaTeX string:
#' # Using writeLines() to format
#' writeLines(iwishmom_sym(4, 0, 1, latex=TRUE)) # iw = 0, for complex Wishart distribution
#'
#' # Example 4: For E[tr(W^{-1})*tr(W^{-2})W^{-1}] with W ~ W_m^2(n,S), represented as a LaTeX string:
#' # Using writeLines() to format
#' writeLines(iwishmom_sym(c(1, 1), 1, 1, latex=TRUE)) # iw = 1, for real Wishart distribution

iwishmom_sym <- function(f, iw = 0, alpha = 2, latex = FALSE) {
  lf <- length(f)
  k1 <- sum(f * seq_len(lf))
  k <- k1 + iw

  #  Fast return for k=1
  if (k == 1) {
    if (iw == 1) {
      if (latex) {
        return("\\Sigma^{-1}/\\tilde{n}")
      }
      results <- data.frame(
        i = 1,
        rho = 1,
        c_numerator = 1
      )
      return(list(c = results, denominator = "n1"))
    } else {
      if (latex) {
        return("p_{(1)}(\\Sigma^{-1})/\\tilde{n}")
      }
      results <- data.frame(
        kappa = "1",
        h_kappa_numerator = "1"
      )
      return(list(c = results, denominator = "n1"))
    }
  }

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
#  polynomial of n1, with the coefficients of declining
#  powers of n1 stored in v, with the last element being the
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
            expr <- paste0(expr, "n1", powers(lv-i))
          } else if (coeff == -1) {
            expr <- paste0(expr, "-n1", powers(lv-i))
          } else {
            expr <- paste0(expr, coeff, "n1", powers(lv-i))
          }
        }
      }
    }
    return(expr)
  }

  if (iw == 0) {  # The output is a scalar for this case
    c <- iwish_psr(k, alpha)
    den <- c[["den"]]
    den <- n_poly(c(den,0))
    c <- c[["c"]]
    nH <- dim(c)[3]
    c <- c[ind(f), , ]
    partitions <- ip_desc(k)
    kappa_lst <- c()
    num_lst <- c()
    for (i in 1:nrow(c)) {
      partition <- partitions[i,]
      kappa_lst <- c(kappa_lst, paste0(partition[partition!=0], collapse=","))
      num_lst <- c(num_lst, n_poly(c[i, ]))
    }

    results <- data.frame(
      kappa = kappa_lst,
      h_kappa_numerator = num_lst
    )
    if (!latex) return(list(dataframe = results, denominator = den))

    latex_expr <- paste0("(", results[1, "h_kappa_numerator"], ")p_{(", results[1, 1], ")}(\\Sigma^{-1})")
    if (nrow(results) == 1) return(latex_expr)
    for (i in 2:nrow(results)) {
      latex_expr <- paste0(latex_expr, "\n+(", results[i, "h_kappa_numerator"], ")p_{(", results[i, 1], ")}(\\Sigma^{-1})")
    }
  } else {
    c <- qkn_coeffr(k, alpha)
    den <- c[["den"]]
    den <- n_poly(c(den,0))
    c <- c[["c"]]
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
      c_numerator = coeff_lst
    )

    if (!latex) return(list(dataframe = results, denominator = den))

    formatted_list <- list()
    i <- 1
    group_string <- ""
    for (index in 1:nrow(results)) {
      row <- results[index, ]
      if (row["i"] > i) {
        if (i == (k-1)) {
          group_string <- paste0(group_string, "\\Sigma^{-", k-1, "}\n")
        } else {
          group_string <- paste0("[", group_string, "]\\Sigma^{-", i, "}\n")
        }
        formatted_list <- c(formatted_list, group_string)
        i <- i + 1
        group_string <- ""
      }
      if (i == k){
        group_string <- paste0("(", row["c_numerator"], ")\\Sigma^{-", k, "}")
        formatted_list <- c(formatted_list, group_string)
      } else {
        if (nchar(group_string) > 0) group_string <- paste0(group_string, "+")
        group_string <- paste0(group_string, "(", row["c_numerator"], ")p_{(", row["rho"], ")}(\\Sigma^{-1})")
      }
    }
    latex_expr <- paste(formatted_list, collapse = "+")
  }
  latex_expr <- gsub("n1", "\\\\tilde{n}", latex_expr)
  den <- gsub("n1", "\\\\tilde{n}", den)
  return(paste0("[", latex_expr, "]\n/(", den, ")"))
}
