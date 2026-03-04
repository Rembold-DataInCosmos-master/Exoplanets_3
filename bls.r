library(tidyverse)

bls_search <- function(t, y, periods, n_bins = 200, d_min = 0.01, d_max = 0.1) {
  N <- length(y)
  num_periods <- length(periods)
  
  # Pre-calculate bin counts for duration L
  L_min <- max(1, floor(d_min * n_bins))
  L_max <- max(L_min + 1, floor(d_max * n_bins))
  L_range <- L_min:L_max
  
  # Pre-allocate vectors for speed
  best_powers <- numeric(num_periods)
  best_L      <- integer(num_periods)
  best_start  <- integer(num_periods)
  best_s      <- numeric(num_periods)
  best_r      <- numeric(num_periods)
 
  # Loop over periods
  for (p_idx in seq_along(periods)) {
    P <- periods[p_idx]
    
    # Binning
    phase <- (t %% P) / P
    bin_idx <- floor(phase * n_bins) + 1
    bin_idx[bin_idx > n_bins] <- n_bins
    
    # Counts (r) and Weighted Sums (s)
    r_bins <- tabulate(bin_idx, nbins = n_bins)
    s_bins <- as.numeric(tapply(y, bin_idx, sum))
    s_bins[is.na(s_bins)] <- 0
    
    # 2. Tiling & Cumulative Sums
    S <- cumsum(c(0, s_bins, s_bins))
    R <- cumsum(c(0, r_bins, r_bins))
    
    # Initialize trackers for this specific Period
    max_pwr_P <- -1
    
    # Loop over durations
    for (L in L_range) {
      i_idx <- 1:n_bins # We vectorize the position, to compute all at once
      j_idx <- i_idx + L
      
      s_sum <- S[j_idx + 1] - S[i_idx]
      r_sum <- R[j_idx + 1] - R[i_idx]
      
      # Calculate Power, ensuring only dips counting
      mask <- (r_sum > 0) & (r_sum < N) & (s_sum < 0)
      
      if(any(mask)) {
        pwr <- rep(0, n_bins)
        pwr[mask] <- (s_sum[mask]^2) / (r_sum[mask] * (N - r_sum[mask]))
        
        # Identify the best start position for this L
        loc_max_idx <- which.max(pwr)
        loc_max_pwr <- pwr[loc_max_idx]
        
        # If this is the best L we've seen for this Period, save the metadata
        if (loc_max_pwr > max_pwr_P) {
          max_pwr_P <- loc_max_pwr
          best_powers[p_idx] <- loc_max_pwr
          best_L[p_idx]      <- L
          best_start[p_idx]  <- loc_max_idx
          best_s[p_idx]      <- s_sum[loc_max_idx]
          best_r[p_idx]      <- r_sum[loc_max_idx]
        }
      }
    }
  }
  
  # 4. Return the Final Results Table
  return(tibble(
    period = periods,
    power = best_powers,
    duration_bins = best_L,
    start_bin = best_start,
    s = best_s,
    r = best_r
  ))
}
