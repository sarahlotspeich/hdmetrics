single_time_uno_C = function(t_star, time, event, log_risk, data, G_hat) {
  # Extract min times, event indicators, linear predictors, and sample size
  W  = data[, time]
  d  = data[, event]
  lp = data[, log_risk]
  n  = length(W)

  # Create all pairwise combinations, dropping i == j
  pair_df = expand.grid(i = 1:n, j = 1:n)
  pair_df = subset(pair_df, i != j)

  # Extract values for each pair
  Ti = W[pair_df$i]
  Tj = W[pair_df$j]
  di = d[pair_df$i]
  lpi = lp[pair_df$i]
  lpj = lp[pair_df$j]

  # Identify valid comparable pairs (Equation (6) in Uno et al. (2011))
  valid_pairs = which(Ti < Tj & Ti <= t_star & di == 1)
  if (length(valid_pairs) == 0) return(NA)

  # Calculate inverse probability of censoring weights
  G_Ti = G_hat(Ti[valid_pairs])
  wi = 1 / (G_Ti ^ 2)

  # Concordance: a pair is concordant if the LP[i] > LP[j]
  concordant = as.numeric(lpi[valid_pairs] > lpj[valid_pairs])

  # Compute Uno's C (Equation (6) in Uno et al. (2011))
  sum(wi * concordant) / sum(wi)
}
