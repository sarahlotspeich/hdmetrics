#' Kaplan-Meier adjusted AUROC
#'
#' This function returns the KM-adjusted area under the receiver operating characteristic curves for survival models.
#'
#' @param evaluate_at scalar or vector, numeric times at which Uno's C will be calculated.
#' @param time string variable name, either time to event or end of follow-up.
#' @param event string variable name, indicator of whether \code{time} corresponds to an event (\code{event=1}) or censoring (\code{event=0}).
#' @param log_risk string variable name, log risk score predicted from model.
#' @param data dataframe, containing at least \code{time} and \code{event}.
#' @return numeric AUC value
#' @export
#' 
#' @importFrom survivalROC survivalROC
#'
km_auroc = function(evaluate_at, time, event, log_risk, data) {
  
  auc_val = survivalROC(
    Stime = data[, time],
    status = data[, event],
    marker = data[, log_risk],
    predict.time = evaluate_at,
    method = "KM"
  )$AUC
  
  
}