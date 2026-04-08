
#' Kaplan-Meier adjusted AUROC
#'
#' This function returns the KM-adjusted area under the receiver operating characteristic curves for survival models.
#'
#' @param evaluate_at scalar, vector, numeric times at which Uno's C will be calculated.
#' @param time string variable name, either time to event or end of follow-up.
#' @param event string variable name, indicator of whether \code{time} corresponds to an event (\code{event=1}) or censoring (\code{event=0}).
#' @param log_risk string variable name, log risk score predicted from model.
#' @param data dataframe, containing at least \code{time} and \code{event}.
#' @return numeric AUC value
#' @export
#' 
#'
km_auroc <- function(evaluate_at, time, event, log_risk, data) {
  
  Stime <- data[, time]
  status <- data[, event]
  marker <- data[, log_risk]
  
  # marginal survival at "evaluate_at"
  km_fit <- survfit(Surv(Stime, status) ~ 1)
  
  all_time_km_auroc <- sapply(
    X = evaluate_at,
    FUN = single_time_km_auroc,
    km_fit = km_fit,
    marker = marker,
    Stime = Stime,
    status = status
  )
  
  
  return(data.frame(x = evaluate_at,
                    auc = all_time_km_auroc))
  
}
