compute_entropy = function(probs) {
  return(-1 * sum(probs * log(probs)))
}


estimate_mutual_information_BinCont = function(delta_S, delta_T) {
  # Estimate three conditional densities for the three possible values of delta
  # T.
  lower_S = min(delta_S)
  upper_S = max(delta_S)
  range_S = upper_S - lower_S
  lower_S = lower_S - 0.2 * range_S
  upper_S = upper_S + 0.2 * range_S

  dens_min1 = density(delta_S[delta_T == -1],
                      from = lower_S,
                      to = upper_S,
                      n = 1024)
  dens_0 = density(delta_S[delta_T == 0],
                   from = lower_S,
                   to = upper_S,
                   n = 1024)
  dens_plus1 = density(delta_S[delta_T == 1],
                       from = lower_S,
                       to = upper_S,
                       n = 1024)

  # Convert the estimated densities to R functions.
  densfun_min1 = approxfun(dens_min1$x, dens_min1$y)
  densfun_0 = approxfun(dens_0$x, dens_0$y)
  densfun_plus1 = approxfun(dens_plus1$x, dens_plus1$y)

  # Compute marginal probabilities for distribution of Delta T.
  pi_min1 = mean(delta_T == -1)
  pi_0 = mean(delta_T == 0)
  pi_plus1 = 1 - pi_min1 - pi_0

  # Construct marginal density function of S.
  densfun_marg = function(x) {
    pi_min1 * densfun_min1(x) +
      pi_0 * densfun_0(x) +
      pi_plus1 * densfun_plus1(x)
  }

  # Compute integral for each conditional density using the trapezoidal rule.
  integral_min1 = pi_min1 * cubature::cubintegrate(
    f = function(x) {
      y = densfun_min1(x)
      ifelse(y == 0,
             0,
             y * log(y))
    },
    lower = lower_S,
    upper = upper_S
  )$integral
  integral_0 = pi_0 * cubature::cubintegrate(
    f = function(x) {
      y = densfun_0(x)
      ifelse(y == 0,
             0,
             y * log(y))
    },
    lower = lower_S,
    upper = upper_S
  )$integral
  integral_plus1 = pi_plus1 * cubature::cubintegrate(
    f = function(x) {
      y = densfun_plus1(x)
      ifelse(y == 0,
             0,
             y * log(y))
    },
    lower = lower_S,
    upper = upper_S
  )$integral

  # Compute differential entropy of Delta S
  diff_entropy = cubature::cubintegrate(
    f = function(x) {
      y = densfun_marg(x)
      ifelse(y == 0,
             0,
             y * log(y))
    },
    lower = lower_S,
    upper = upper_S
  )$integral

  mutual_information = integral_min1 + integral_0 + integral_plus1 - diff_entropy
  return(mutual_information)
}

#' Estimate ICA in Binary-Continuous Setting
#'
#' `estimate_ICA_BinCont()` estimates the individual causal association (ICA)
#' for a sample of individual causal treatment effects with a continuous
#' surrogate and a binary true endpoint. The ICA in this setting is defined as
#' follows, \deqn{R^2_H = \frac{I(\Delta S; \Delta T)}{H(\Delta T)}} where
#' \eqn{I(\Delta S; \Delta T)} is the mutual information and \eqn{H(\Delta T)}
#' the entropy.
#'
#' @param delta_S (numeric) Vector of individual causal treatment effects on the
#'   surrogate.
#' @param delta_T (integer) Vector of individual causal treatment effects on the true
#'   endpoint. Should take on one of the following values: `-1L`, `0L`, or `1L`.
#'
#' @return (numeric) Estimated ICA

estimate_ICA_BinCont = function(delta_S, delta_T) {
  # Compute marginal probabilities for distribution of Delta T.
  pi_min1 = mean(delta_T == -1)
  pi_0 = mean(delta_T == 0)
  pi_plus1 = 1 - pi_min1 - pi_0
  # Compute ICA
  ICA = estimate_mutual_information_BinCont(delta_S, delta_T) /
    compute_entropy(c(pi_min1, pi_0, pi_plus1))
  return(ICA)
}


#' Sample copula data from a given four-dimensional D-vine copula
#'
#' `sample_dvine()` is a helper function that samples copula data from a given
#' D-vine copula. See details for more information on the parameterization of
#' the D-vine copula.
#'
#' @details
#'
#' # D-vine Copula
#'
#' Let \eqn{\boldsymbol{U} = (U_1, U_2, U_3, U_4)'} be a random vector with
#' uniform margins. The corresponding distribution function is then a
#' 4-dimensional copula. A D-vine copula as a family of \eqn{k}-dimensional
#' copulas. Indeed, a D-vine copula is a \eqn{k}-dimensional copula that is
#' constructed from a particular product of bivariate copula densities. In this
#' function, only 4-dimensional copula densities are considered. Under the
#' simplifying assumption, the 4-dimensional D-vine copula density is the
#' product of the following bivariate copula densities:
#' * \eqn{c_{12}}, \eqn{c_{23}}, and \eqn{c_{34}}
#' * \eqn{c_{13;2}} and \eqn{c_{24;3}}
#' * \eqn{c_{14;23}}
#'
#'
#' @param copula_par Parameter vector for the sequence of bivariate copulas that
#'   define the D-vine copula. The elements of `copula_par` correspond to
#'   \eqn{(c_{12}, c_{23}, c_{34}, c_{13;2}, c_{24;3}, c_{14;23})}.
#' @param rotation_par Vector of rotation parameters for the sequence of
#'   bivariate copulas that define the D-vine copula. The elements of
#'   `rotation_par` correspond to \eqn{(c_{12}, c_{23}, c_{34}, c_{13;2},
#'   c_{24;3}, c_{14;23})}.
#' @param copula_family1 Copula family of \eqn{c_{12}} and \eqn{c_{34}}. For the
#'   possible options, see `loglik_copula_scale()`.
#' @param copula_family2 Copula family of the other bivariate copulas. For the
#'   possible options, see `loglik_copula_scale()`.
#' @param n Number of samples to be taken from the D-vine copula.
#'
#'
#' @return A \eqn{n \times 4} matrix where each row corresponds to one sampled
#'   vector and the columns correspond to \eqn{U_1}, \eqn{U_2}, \eqn{U_3}, and
#'   \eqn{U_4}.
#'
#' @examples
sample_dvine = function(copula_par,
                        rotation_par,
                        copula_family1,
                        copula_family2 = copula_family1,
                        n) {
  # D-vine copula parameters
  c12 = copula_par[1]
  c34 = copula_par[2]
  c23 = copula_par[3]
  c13_2 = copula_par[4]
  c24_3 = copula_par[5]
  c14_23 = copula_par[6]

  # D-vine copula rotations
  r12 = rotation_par[1]
  r34 = rotation_par[2]
  r23 = rotation_par[3]
  r13_2 = rotation_par[4]
  r24_3 = rotation_par[5]
  r14_23 = rotation_par[6]

  pair_copulas = list(
    list(
      rvinecopulib::bicop_dist(
        family = copula_family1,
        rotation = r12,
        parameters = c12
      ),
      rvinecopulib::bicop_dist(
        family = copula_family2,
        rotation = r23,
        parameters = c23
      ),
      rvinecopulib::bicop_dist(
        family = copula_family1,
        rotation = r34,
        parameters = c34
      )
    ),
    list(
      rvinecopulib::bicop_dist(
        family = copula_family2,
        rotation = r13_2,
        parameters = c13_2
      ),
      rvinecopulib::bicop_dist(
        family = copula_family2,
        rotation = r24_3,
        parameters = c24_3
      )
    ),
    list(
      rvinecopulib::bicop_dist(
        family = copula_family2,
        rotation = r14_23,
        parameters = c14_23
      )
    )
  )
  copula_structure = rvinecopulib::dvine_structure(order = 1:4)
  vine_cop = rvinecopulib::vinecop_dist(pair_copulas = pair_copulas, structure = copula_structure)

  U = rvinecopulib::rvinecop(n = n, vinecop = vine_cop, cores = 1)


  return(U)
}

#' Sample individual casual treatment effects from given D-vine copula model in
#' binary continuous setting
#'
#' @param q_S0 Quantile function for the distribution of \eqn{S_0}.
#' @param q_S1 Quantile function for the distribution of \eqn{S_1}.
#' @inheritParams sample_dvine
#'
#' @return A data frame
sample_deltas_BinCont = function(copula_par,
                                 rotation_par,
                                 copula_family1,
                                 copula_family2 = copula_family1,
                                 n,
                                 q_S0,
                                 q_S1){
  # Sample data on the copula scale.
  U = sample_dvine(copula_par,
                   rotation_par,
                   copula_family1,
                   copula_family2,
                   n)
  # Convert copula data to the original scale. For the true endpoints, we
  # convert from to the copula scale to the latent standard normal scale to the
  # binary scale. For the surrogate endpoint, we convert from the copula scale
  # to the continuous scale, as defined by the corresponding marginal
  # distributions.
  T0 = ifelse(qnorm(U[, 1]) < 0, 0L, 1L)
  T1 = ifelse(qnorm(U[, 4]) < 0, 0L, 1L)
  S0 = q_S0(U[, 2])
  S1 = q_S1(U[, 3])

  return(
    data.frame(
      DeltaS = S1 - S0,
      DeltaT = T1 - T0
    )
  )
}

#' Compute Individual Causal Association for a given D-vine copula model in the
#' Binary-Continuous Setting
#'
#' The [compute_ICA_BinCont()] function computes the individual causal
#' association for a fully identified D-vine copula model in the setting with a
#' continuous surrogate endpoint and a binary true endpoint.
#'
#' @param minfo_prec Number of Monte Carlo samples for the computation of the
#'   mutual information.
#' @inheritParams sample_deltas_BinCont
#'
#' @return (numeric) A Named vector with the following elements:
#'  * ICA
#'  * Spearman's rho, \eqn{\rho_s (\Delta S, \Delta T)} (if asked)
#'  * Kendall's tau, \eqn{\tau (\Delta S, \Delta T)} (if asked)
#'  * Marginal association parameters in terms of Spearman's rho:
#'  \deqn{(\rho_s(S_0, S_1), \rho_s(S_0, S_T_0), \rho_s(S_0, T_1),
#'  \rho_s(S_1, T_0), \rho_s(S_0, S_1),
#'  \rho_s(T_0, T_1)}
#'
#' @examples
compute_ICA_BinCont = function(copula_par,
                               rotation_par,
                               copula_family1,
                               copula_family2 = copula_family1,
                               minfo_prec,
                               q_S0,
                               q_S1)
{
  delta_df = sample_deltas_BinCont(
    copula_par,
    rotation_par,
    copula_family1,
    copula_family2 = copula_family1,
    n = minfo_prec,
    q_S0,
    q_S1
  )
  # Compute mutual information between Delta S and Delta T.
  mutual_information = estimate_mutual_information_BinCont(delta_df$DeltaS, delta_df$DeltaT)
  # Compute marginal probabilities for distribution of Delta T.
  pi_min1 = mean(delta_df$DeltaT == -1)
  pi_0 = mean(delta_df$DeltaT == 0)
  pi_plus1 = 1 - pi_min1 - pi_0
  entropy_DeltaT = compute_entropy(c(pi_min1, pi_0, pi_plus1))
  # Compute ICA
  ICA = mutual_information / entropy_DeltaT
  sp_rho = cor(delta_df$DeltaS, delta_df$DeltaT, method = "spearman")
  kendall_tau = cor(delta_df$DeltaS, delta_df$DeltaT, method = "kendall")
  return(c(ICA = ICA, sp_rho = sp_rho))
}