import numpy as np
import pandas as pd

def bls_search(t, y, periods, n_bins=200, d_min=0.01, d_max=0.1):
    """Conducts a Box Least Squares search to find power across a
    provided grid.

    Args:
        t (np.array): Array of the times of observation
        y (np.array): Array of flattened and zero-based observations
        periods (np.array): Array of periods to check
        n_bins (int): The number of bins to divide the folded signal into
        d_min (float): The minimum size of a duration, in normalized phase units
        d_max (float): The maximum size of a duration, in normalized phase units

    Returns
        (pd.DataFrame): A dataframe containing best values for each period
    """

    N = len(y)
    results = []
   
    # For binning to work well, we constrain our starting and ending durations into
    # multiples of bins. This might not work out EXACTLY to the desired start and
    # stops, but it is close enough and provides big speed gains
    L_min = max(1, int(d_min * n_bins))
    L_max = max(L_min + 1, int(d_max * n_bins))
   
    # Looping over all the periods
    for P in periods:
        # Binning & Cumulative Summing for increased speed
        phase = (t % P) / P
        bin_idx = np.clip((phase * n_bins).astype(int), 0, n_bins - 1)
        s_bins = np.bincount(bin_idx, weights=y, minlength=n_bins)
        r_bins = np.bincount(bin_idx, minlength=n_bins)
        
        S = np.cumsum(np.insert(np.tile(s_bins, 2), 0, 0))
        R = np.cumsum(np.insert(np.tile(r_bins, 2), 0, 0))
       
        best_pwr_for_P = -1
        best_metadata = {}

        # Looping over all the possible durations
        for L in range(L_min, L_max + 1):
            # We vectorize to compute all the possible box positions simultaneously
            i = np.arange(n_bins)
            j = i + L
            s_sum = S[j] - S[i]
            r_sum = R[j] - R[i]
            
            # Calculate Power, tracking only dips
            mask = (r_sum > 0) & (r_sum < N) & (s_sum < 0)
            pwr = np.zeros(n_bins)
            pwr[mask] = (s_sum[mask]**2) / (r_sum[mask] * (N - r_sum[mask]))
            
            # Find the best start position 'i' for THIS specific L
            current_max_idx = np.argmax(pwr)
            current_max_pwr = pwr[current_max_idx]
            
            if current_max_pwr > best_pwr_for_P:
                best_pwr_for_P = current_max_pwr
                # Store the raw ingredients
                best_metadata = {
                    'period': P,
                    'power': best_pwr_for_P,
                    'duration_bins': L,
                    'start_bin': current_max_idx,
                    's': s_sum[current_max_idx],
                    'r': r_sum[current_max_idx]
                }
        
        results.append(best_metadata)
        
    return pd.DataFrame(results)
