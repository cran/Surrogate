#' Assess surrogacy in the information-theoretic causal-inference framework (Individual Causal Association, ICA)
#' for a continuous surrogate and true endpoint measured repeatedly over time in a single-trial setting
#'
#'
#' The function quantifies surrogacy under linear model
#' with subject-specific random intercepts with separable residual covariance structure. See Details below.
#'
#' @param V0,V1 Treatment-specific residual 2 x 2 covariance matrices for T,S.
#' @param D0,D1 Treatment-specific random-intercept 2 x 2 covariance matrices for T,S.
#' @param times Clinical evaluation time points.
#' @param correlation Temporal correlation structures such as exponential, compound_symmetry, and independence.
#' @param temporal_rho Temporal-correlation parameter. Ignored when \code{correlation = "independence"}
#' @param grid_V Sensitivity grid for unidentified residual correlations.
#' @param grid_D Sensitivity grid for unidentified random-intercept correlations.
#' @param tol Numerical tolerance used when checking positive definiteness and validating ICA values.
#' @param return_candidates Logical. If \code{TRUE}, the retained admissible covariance
#'  completion candidates for the residual and random-intercept components are returned.
#'
#' @return
#' * Num.Pos.Def.V: An object of class numeric that contains the number of
#'   positive definite `V` matrices.
#' * Num.Pos.Def.D: An object of class numeric that contains the number of
#'   positive definite `D` matrices.
#' * Num.Pos.Def.Pairs: An object of class numeric that contains the number of
#'   valid `V-D` pairs used in the final calculations.
#' * R2_Lambda: A scalar or vector that contains the individual causal association \eqn{R_{\Lambda}^2}.
#' * V_candidates: retained residual covariance candidates, if requested.
#' * D_candidates: retained random-intercept covariance candidates, if requested.
#'
#' @details
#' # Random intercepts with separable residual covariance structure:
#' The evolution of \eqn{y_j} over time is defined using a linear model with subject-specific random intercept
#'  for each potential outcome.
#' \deqn{T_{0j} = x_{T0j} \gamma_{T0} + b_{T0}+ \varepsilon_{T0j}}
#' \deqn{T_{1j} = x_{T1j} \gamma_{T1} + b_{T1}+ \varepsilon_{T1j}}
#' \deqn{S_{0j} = x_{S0j} \gamma_{S0} + b_{S0}+ \varepsilon_{S0j}}
#' \deqn{S_{1j} = x_{S1j} \gamma_{S1} + b_{S1}+ \varepsilon_{S1j}}
#'
#'It can be written as \eqn{y=x \gamma + b+ \varepsilon} where \eqn{b \sim N(0, D)}
#' with \eqn{D} an unstructured \eqn{4 \times 4} covariance matrix and \eqn{\varepsilon \sim N(0, \Sigma)}.
#' The marginal model follows that \eqn{y \sim N(0, \Sigma_y)} where \eqn{\Sigma_y= zDz+ \Sigma}
#' with \eqn{\Sigma = R \otimes V}.
#' \eqn{V} is assumed to be unstructured and
#' invariant across time, and \eqn{R} reflects a general correlation matrix.
#'
#' The marginal covariance matrix of the individual causal treatment effects:
#' \deqn{\Sigma_\Delta= \Sigma_u \otimes \R + \Sigma_D\otimes 1_p 1_p^\top}
#'
#'with \eqn{\Sigma_u=QVQ^\top} and \eqn{\Sigma_D=QDQ^\top} a \eqn{2\times 2} matrix.
#' Equivalently, \eqn{\Sigma_{\Delta}} has blocks
#' \deqn{\Sigma_{\Delta T\Delta T}=d_{11}(1_{p}1^\top_{p})+u_{11}R,\quad
#'  \Sigma_{\Delta T\Delta S}=d_{12}(1_{p}1^\top_{p})+u_{12}R,\quad
#'  \Sigma_{\Delta S\Delta S}=d_{22}(1_{p}1^\top_{p})+u_{22}R }
#'
#' Under the random-intercept model specified, the ICA metric is given by:
#' \deqn{ R_{\Lambda}^{2}=1-(1-\rho_{u}^{2})^{p-1}\,(1-\rho_{ud}^{2}) }
#' where
#'
#'\deqn{\rho_{u}^{2}
#'=\frac{u_{12}^{2}}{u_{11}\,u_{22}},
#'  \qquad
#'  \rho_{ud}^{2}
#'  =\frac{(u_{12}+\alpha\,d_{12})^{2}}
#'  {(u_{11}+\alpha\,d_{11})(u_{22}+\alpha\,d_{22})},
#'  \qquad\mbox{and}\qquad
#'  \alpha=1_{p}^{\top}R^{-1}1_{p}}
#'
#'Here \eqn{u_{ij}} are the entries of \eqn{\Sigma_{u}}
#'and \eqn{d_{ij}} those of \eqn{\Sigma_{D}}.
#'
#'
#' @references
#' Deliorman, G., Pardo, M. C., Van der Elst W., Molenberghs G., & Alonso, A. (submitted in 2026).
#' Assessing Surrogate Endpoints in Longitudinal Studies: An Information-Theoretic
#' Causal Inference Approach.
#'
#'
#' @examples
#' #Example outputs
#'
#' V0 <- matrix(c(0.02579325, 0.02514841,
#'               0.02514841, 0.02872256),
#'               nrow = 2, byrow = TRUE,
#'               dimnames = list(c("T", "S"), c("T", "S")))
#'
#'
#' V1 <- matrix(c(0.02907611, 0.02962312,
#'                0.02962312, 0.03400179),
#'                nrow = 2, byrow = TRUE,
#'                dimnames = list(c("T", "S"), c("T", "S")))
#'
#' D0 <- matrix(c(0.02587092, 0.02458997,
#'                0.02458997, 0.02405581),
#'                nrow = 2, byrow = TRUE,
#'                dimnames = list(c("T", "S"), c("T", "S")))
#'
#' D1 <- matrix(c(0.03649246, 0.03225987,
#'               0.03225987, 0.02897932),
#'               nrow = 2, byrow = TRUE,
#'               dimnames = list(c("T", "S"), c("T", "S")))
#'
#' # Example run
#' ica_ri<- ICA_contcont_long_ri(V0 = V0,  V1 = V1,
#'                      D0 = D0,  D1 = D1,
#'                      times=c(1,2,4,6,8),
#'                      correlation  = "exponential",
#'                      temporal_rho = 0.7117544 ,
#'                      grid_V = seq(-1, 1, by = 0.5),
#'                      grid_D = seq(-1, 1, by = 0.5),
#'                      return_candidates = FALSE,
#'                      tol = 1e-12)
#'
#' summary_longitudinal_ica(ica_ri)
#'
#' @export
ICA_contcont_long_ri <- function(V0,
                                 V1,
                                 D0,
                                 D1,
                                 times,
                                 correlation = c("exponential", "compound_symmetry", "independence"),
                                 temporal_rho = NULL,
                                 grid_V = seq(-1, 1, by = 0.1),
                                 grid_D = seq(-1, 1, by = 0.1),
                                 return_candidates = FALSE,
                                 tol = 1e-12) {
  times <- as.numeric(times)
  correlation <- match.arg(correlation)

  if (length(times) < 1L || any(!is.finite(times)) || anyDuplicated(times)) {
    stop("`times` must contain distinct finite numeric values.")
  }


  if (correlation %in% c("exponential", "compound_symmetry")) {
    if (is.null(temporal_rho) ||
        length(temporal_rho) != 1L ||
        !is.finite(temporal_rho)) {
      stop("`temporal_rho` must be supplied as a single finite value when correlation is 'exponential' or 'compound_symmetry'.")
    }
  }

  if (correlation == "exponential") {
    if (temporal_rho <= 0 || temporal_rho >= 1) {
      stop("For exponential correlation, `temporal_rho` must lie in (0,1).")
    }
  }

  if (correlation == "compound_symmetry") {
    lower <- if (length(times) <= 1L) -1 else -1 / (length(times) - 1L)
    if (temporal_rho <= lower || temporal_rho >= 1) {
      stop("For compound symmetry, `temporal_rho` must lie in (", signif(lower, 5), ", 1).")
    }
  }

  if (correlation == "independence") {
    temporal_rho <- NA_real_
  }


  Vc <- .surr_completion_components(V0, V1, grid_V, tol = tol,
                                    return_candidates = return_candidates)
  Dc <- .surr_completion_components(D0, D1, grid_D, tol = tol,
                                    return_candidates = return_candidates)
  U <- Vc$components
  DD <- Dc$components

  p <- length(times)


  if (correlation == "exponential") {
    R <- make_longitudinal_R(times, rho = temporal_rho, type = "exponential")
  }

  if (correlation == "compound_symmetry") {
    R <- make_longitudinal_R(times, rho = temporal_rho, type = "cs")
  }

  if (correlation == "independence") {
    R <- make_longitudinal_R(times, rho = NULL, type = "independence")
  }

  one <- rep(1, p)
  alpha <- as.numeric(crossprod(one, solve(R, one)))

  max_pairs <- nrow(U) * nrow(DD)
  r2 <- numeric(max_pairs)
  out_index <- 0L

  if (max_pairs > 0L) {
    d11 <- DD$c11; d12 <- DD$c12; d22 <- DD$c22
    for (i in seq_len(nrow(U))) {
      u11 <- U$c11[i]; u12 <- U$c12[i]; u22 <- U$c22[i]
      rho_u2 <- U$rho2[i]
      denom <- (u11 + alpha * d11) * (u22 + alpha * d22)
      keep <- is.finite(denom) & denom > tol
      if (!any(keep)) next

      rho_ud2 <- (u12 + alpha * d12[keep])^2 / denom[keep]
      vals <- 1 - (1 - rho_u2)^(p - 1L) * (1 - rho_ud2)
      vals <- .surr_validate_R2(vals)
      idx <- out_index + seq_along(vals)
      r2[idx] <- vals
      out_index <- out_index + length(vals)
    }
  }
  if (out_index == 0L) r2 <- numeric(0) else r2 <- r2[seq_len(out_index)]

  out <- list(
    Num.Pos.Def.V = nrow(U),
    Num.Pos.Def.D = nrow(DD),
    Num.Pos.Def.Pairs = length(r2),
    temporal_rho = temporal_rho,
    alpha = alpha,
    correlation = correlation,
    R2_Lambda = r2,
    V_candidates = if (return_candidates) U else NULL,
    D_candidates = if (return_candidates) DD else NULL,
    Call = match.call()
  )
  class(out) <- c("ICA_contcont_long_ri", "surrogate_longitudinal_ica")
  out
}

