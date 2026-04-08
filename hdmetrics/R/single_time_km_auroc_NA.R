single_time_km_auroc_NA <- function(km_fit, evaluate_at, marker, Stime, status) {
  
  km_step <- function(fit, t) {
    idx <- which(fit$time <= t)
    if (length(idx) == 0) return(1)
    fit$surv[max(idx)]
  }
  
  S_t <- km_step(km_fit, evaluate_at)
  F_t <- 1 - S_t
  
  cut_vals <- sort(unique(marker))
  n        <- length(marker)
  
  # precompute KM quantities once for all observations
  # sort event times
  ord       <- order(Stime)
  t_ord     <- Stime[ord]
  s_ord     <- status[ord]
  m_ord     <- marker[ord]
  
  # unique event times up to evaluate_at
  event_times <- unique(t_ord[s_ord == 1 & t_ord <= evaluate_at])
  event_times <- sort(event_times)
  
  # for each cutpoint, compute S(t* | marker > c) via KM product-limit
  # vectorized over cutpoints using outer()
  
  # for each event time and each cutpoint:
  # at_risk[i, k] = number at risk at time t_i among marker > cut_vals[k]
  # events[i, k]  = number of events at time t_i among marker > cut_vals[k]
  
  marker_mat  <- outer(marker, cut_vals, ">")        # n x n_cuts logical
  atrisk_mat  <- outer(Stime, event_times, ">=")     # n x n_events logical
  event_mat   <- outer(Stime, event_times, "==")     # n x n_events logical
  event_ind   <- matrix(status == 1, nrow = n, 
                        ncol = length(event_times))  # n x n_events logical
  
  # at risk among marker > c: colSums over individuals
  # dim: n_events x n_cuts
  denom <- t(atrisk_mat) %*% marker_mat              # n_events x n_cuts
  numer <- t(event_mat * event_ind) %*% marker_mat   # n_events x n_cuts
  
  # KM product limit: prod(1 - d/n) over event times, per cutpoint
  # avoid 0/0
  haz        <- ifelse(denom > 0, numer / denom, 0)  # n_events x n_cuts
  S_t_c      <- apply(1 - haz, 2, prod)              # n_cuts vector
  
  # TP and FP
  p_c <- colMeans(marker_mat)                        # P(marker > c)
  TP  <- ifelse(F_t > 0, (p_c - p_c * S_t_c) / F_t, 0)
  FP  <- p_c
  
  # endpoints, sort, trapezoid
  FP  <- c(0, FP, 1)
  TP  <- c(0, TP, 1)
  ord <- order(FP)
  FP  <- FP[ord]
  TP  <- TP[ord]
  
  sum(diff(FP) * (head(TP, -1) + tail(TP, -1)) / 2)
}