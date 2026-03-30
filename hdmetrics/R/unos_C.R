#' Uno's C statistic for survival models
#'
#' This function returns the time-varying and global Uno's C statistics for an input survival model's predictions.
#'
#' @param evaluate_at scalar or vector, numeric times at which Uno's C will be calculated.
#' @param time string variable name, either time to event or end of follow-up.
#' @param event string variable name, indicator of whether \code{time} corresponds to an event (\code{event=1}) or censoring (\code{event=0}).
#' @param log_risk string variable name, log risk score predicted from model.
#' @param data dataframe, containing at least \code{time} and \code{event}.
#' @return A list with the following named slots:
#' \item{model_coeff}{dataframe with final model coefficients and standard error estimates (where applicable) for the analysis model.}
#' @export
#' @importFrom survival Surv
#' @importFrom survival survfit

unos_C = function(evaluate_at, time, event, log_risk, data) {
  # Fit Kaplan-Meier for the censoring variable C
  km_form = as.formula(paste0("Surv(time =", time, ", event = (1 - ", event, ")) ~ 1"))
  km_cens = survfit(formula = km_form,
                    data = data)

  # Write function to interpolate from Kaplan-Meier fit
  km_cens_ext = function(t) {
    s <- summary(object = km_cens,
                 times = unique(c(0, km_cens$time)),
                 extend = TRUE)
    approx(
      x = s$time,
      y = s$surv,
      xout = t - 1e-8,
      method = "constant",
      f = 0,
      rule = 2
    )$y
  }

  # Calculate time-varying Uno's C at all times in evaluate_at
  all_time_uno_C = sapply(
    X = evaluate_at,
    FUN = single_time_uno_C,
    time = time,
    event = event,
    log_risk = log_risk,
    data = data,
    G_hat = km_cens_ext)

  # Create event weights for global Uno's C values across times by fold
  weights = sapply(
    X = evaluate_at,
    FUN = function(t) {
      sum(data[, time] == t & data[, event] == 1)
      }
    )
  weights = weights / sum(weights)

  # Summarize Uno's C global values and combine with time-specific data
  all_time_uno_C = weighted.mean(
    x = all_time_uno_C,
    w = weights,
    na.rm = TRUE)

  # Return dataframe with (i) times and Uno's C and (ii) Global
  return(
    list(
      tv_unos_c = data.frame(x = evaluate_at,
                             C = all_time_uno_C),
      global_unos_c = all_time_uno_C
    )
  )
}
