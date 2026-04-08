single_time_km_auroc_prodlim <- function(km_fit, evaluate_at, marker, Stime, status) {
  
  km_step <- function(fit, t) {
    idx <- which(fit$time <= t)
    if (length(idx) == 0) return(1)
    fit$surv[max(idx)]
  }
  
  S_t <- km_step(km_fit, evaluate_at)
  F_t <- 1 - S_t
  
  cut_vals <- sort(unique(marker))
  n_cuts   <- length(cut_vals)
  
  # event times up to evaluate_at
  event_times <- sort(unique(Stime[status == 1 & Stime <= evaluate_at]))
  
  # precompute at-risk and event counts at each event time (full sample)
  # then subset per cutpoint using logical vectors
  
  # product limit for a given subset defined by logical vector `in_grp`
  km_pl <- function(in_grp) {
    if (sum(in_grp) == 0) return(0)
    t_sub <- Stime[in_grp]
    s_sub <- status[in_grp]
    pl    <- 1
    for (t in event_times) {
      d <- sum(s_sub == 1 & t_sub == t)
      n <- sum(t_sub >= t)
      if (n > 0) pl <- pl * (1 - d / n)
    }
    pl
  }
  
  TP <- numeric(n_cuts)
  FP <- numeric(n_cuts)
  
  for (k in seq_len(n_cuts)) {
    in_grp <- marker > cut_vals[k]
    p_c    <- mean(in_grp)
    
    if (sum(in_grp) == 0 || F_t == 0) {
      TP[k] <- 0
      FP[k] <- 0
      next
    }
    
    S_t_c <- km_pl(in_grp)
    TP[k] <- (p_c - p_c * S_t_c) / F_t
    FP[k] <- p_c
  }
  
  FP  <- c(0, FP, 1)
  TP  <- c(0, TP, 1)
  ord <- order(FP)
  FP  <- FP[ord]
  TP  <- TP[ord]
  
  sum(diff(FP) * (head(TP, -1) + tail(TP, -1)) / 2)
}
