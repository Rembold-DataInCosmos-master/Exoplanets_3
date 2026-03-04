# This file compiled a C++ executable that is then linked to R through Rcpp
# It SUBSTANTIALLY improves runtimes, but your system must have a method to
# compile C++ code (Rtools in Windows or Xcode in macOS) This puts the R 
# performance on the same level as the Python implementation, since NumPy 
# already is abstracting out most of the hard word to C code

library(Rcpp)

# Define the C++ Engine
cppFunction('
List bls_search(NumericVector t, NumericVector y, NumericVector periods, 
                      int n_bins, double d_min, double d_max) {
    int n_p = periods.size();
    int n_t = t.size();
    
    // Convert fractional durations to integer bin widths
    int L_min = (int)(d_min * n_bins);
    if (L_min < 1) L_min = 1;
    int L_max = (int)(d_max * n_bins);
    if (L_max <= L_min) L_max = L_min + 1;

    NumericVector out_pwr(n_p);
    IntegerVector out_L(n_p);
    IntegerVector out_start(n_p);
    NumericVector out_s(n_p);
    NumericVector out_r(n_p);

    for (int p = 0; p < n_p; ++p) {
        double P = periods[p];
        NumericVector s_bins(n_bins);
        NumericVector r_bins(n_bins);

        for (int i = 0; i < n_t; ++i) {
            int b = (int)((fmod(t[i], P) / P) * n_bins);
            if (b >= n_bins) b = n_bins - 1;
            if (b < 0) b = 0;
            s_bins[b] += y[i];
            r_bins[b] += 1.0;
        }

        double best_pwr_P = -1.0;

        for (int L = L_min; L <= L_max; ++L) {
            for (int i = 0; i < n_bins; ++i) {
                double s_sum = 0;
                double r_sum = 0;
                
                for (int k = 0; k < L; ++k) {
                    int idx = (i + k) % n_bins;
                    s_sum += s_bins[idx];
                    r_sum += r_bins[idx];
                }

                if (r_sum > 0 && r_sum < n_t && s_sum < 0) {
                    double pwr = (s_sum * s_sum) / (r_sum * (n_t - r_sum));
                    if (pwr > best_pwr_P) {
                        best_pwr_P = pwr;
                        out_pwr[p] = pwr;
                        out_L[p] = L;
                        out_start[p] = i;
                        out_s[p] = s_sum;
                        out_r[p] = r_sum;
                    }
                }
            }
        }
    }

    return List::create(
        Named("period") = periods,
        Named("power") = out_pwr,
        Named("duration_bins") = out_L,
        Named("start_bin") = out_start,
        Named("s") = out_s,
        Named("r") = out_r
    );
}')

print("Note! You will need to wrap the output of this bls_search function in `as_tibble` to convert it back into a dataframe.")
