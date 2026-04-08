
single_time_km_auroc <- function(km_fit, evaluate_at, marker, Stime, event) {
  
  km_step <- function(fit, t) {
    idx <- which(fit$time <= t)
    if (length(idx) == 0) return(1)
    fit$surv[max(idx)]
  }
  
  S_t <- km_step(km_fit, evaluate_at)
  F_t <- 1 - S_t
  
  # cut points
  cut_vals <- sort(unique(marker)) # sort the unique log risk values
  n_cuts <- length(cut_vals)
  TP <- numeric(n_cuts) # true positive allocation
  FP <- numeric(n_cuts) # false positive allocation
  
  # loop over cutpoints  from unique log risk values
  for (k in seq_len(n_cuts)) {
    in_group <- (marker > cut_vals[k])
    p_c <- mean(in_group)
    
    if (sum(in_group) == 0 || F_t == 0) {
      TP[k] <- 0
      FP[k] <- 0
      next
    }
    
    km_c <- survfit(Surv(Stime[in_group], event[in_group]) ~ 1)
    S_t_c <- km_step(km_c, evaluate_at) # survival given log-risk is greater than cut point
    
    TP[k] <- (p_c - p_c*S_t_c) / F_t
    FP[k] <- p_c
  }
  
  # sort cutpoints and add 0, 1
  FP <- c(0, FP, 1)
  TP <- c(0, TP, 1)
  
  ordFP <- order(FP) # order of FP
  FP <- FP[ordFP] # reorder FP
  TP <- TP[ordFP] # reorder in same way as FP to keep paired values together
  
  # get AUC via trapezoidal rule
  auc <- sum(diff(FP) * (head(TP, -1) + tail(TP, -1))/2)
  
  return(auc)
}

