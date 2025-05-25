functions{
  matrix calc_rph_dedup1(real deg, int n_rph_sum_rep1, vector rph_base_use, array[,] int index_dedup1_rep1, vector deg_size_rep1){
    matrix[n_rph_sum_rep1, 2] rph_dedup1;
    vector[n_rph_sum_rep1] rph_param = deg^deg_size_rep1;
    for(j in 1:2){
      rph_dedup1[, j] = rph_base_use[index_dedup1_rep1[j, ]] .* rph_param;
    }
    return rph_dedup1;
  }

  matrix calc_gamma_param(matrix rph_dedup1, matrix cv_square_rep1, int n_roop, int n_rph_sum_rep1, real tph){
    vector[n_rph_sum_rep1] tmp1_sum = rep_vector(0, n_rph_sum_rep1);
    vector[n_rph_sum_rep1] tmp2_sum = rep_vector(0, n_rph_sum_rep1);
    for(i in 1:n_roop){
      tmp1_sum += rph_dedup1[, i];
      tmp2_sum += rph_dedup1[, i] .* cv_square_rep1[, i];
    }

    matrix[n_rph_sum_rep1, 2] gamma_param;
    gamma_param[, 1] = tmp1_sum^2 * tph ./ tmp2_sum;
    gamma_param[, 2] = tmp1_sum ./ tmp2_sum;

    return gamma_param;
  }

  vector calc_lpd_general(int n_general_sum_rep1, vector oph_rep1_general, vector gamma_alpha, vector gamma_beta){
    vector[n_general_sum_rep1] lpd_general;
    for(i in 1:n_general_sum_rep1){
      lpd_general[i] = gamma_lpdf(oph_rep1_general[i] | gamma_alpha[i], gamma_beta[i]);
    }
    return lpd_general;
  }

  vector calc_lpd_dropout(int n_dropout_sum_rep1, vector at_dropout_rep1, vector gamma_alpha, vector gamma_beta){
    vector[n_dropout_sum_rep1] lpd_dropout;
    for(i in 1:n_dropout_sum_rep1){
      lpd_dropout[i] = gamma_lcdf(at_dropout_rep1[i] | gamma_alpha[i], gamma_beta[i]);
    }
    return lpd_dropout;
  }

  row_vector calc_loglike_gtc(int n_peak, int n_gtc, matrix lpd_dedup2, row_vector gtc_logprob_mk){
    row_vector[n_peak] rv = rep_row_vector(1, n_peak);
    row_vector[n_gtc] loglike_gtc = rv * lpd_dedup2 + gtc_logprob_mk;
    return loglike_gtc;
  }
}

data{
  int n_rep;
  int n_mk;
  real tph_mean;
  real tph_sd;

  int n_rph_base_sum;
  vector[n_rph_base_sum] rph_base_use;

  array[n_rep] int n_dedup1_sum_per_rep;

  array[n_rep] int n_rph_sum_per_rep;
  int n_rph_sum;
  array[2, n_rph_sum] int index_dedup1_mcmc;
  matrix[n_rph_sum, 2] cv_square_mcmc;
  vector[n_rph_sum] deg_size_mcmc;
  vector[n_rph_sum] oph_mcmc;

  array[n_rep] int n_general_sum_per_rep;
  array[n_rep] int n_dropout_sum_per_rep;
  int n_dropout_sum;
  vector[n_dropout_sum] at_dropout;
  array[n_rep] int n_dropin_sum_per_rep;
  int n_dropin_sum;
  vector[n_dropin_sum] loglike_di;

  int n_dedup2_sum;
  array[n_rep] int n_dedup2_sum_per_rep;
  array[n_mk] int n_dedup2_use;
  array[n_dedup2_sum] int index_dedup2_use;

  array[n_mk] int n_peak_use;
  array[n_mk] int n_gtc_use;
  row_vector[sum(n_gtc_use)] gtc_logprob_use;
}

parameters{
  vector<lower = 0>[n_rep] tph;
  real<lower = 0, upper = 1> deg;
}

model{
  // Prior distribution for each parameter
  tph ~ normal(tph_mean, tph_sd);
  deg ~ beta(1, 1);

  // Initial indices to access first replicate
  int index_rph_start = 1;
  int index_dropout_start = 1;
  int index_dropin_start = 1;
  int index_dedup2_start = 1;

  for(i in 1:n_rep){
    // Extract data for one replicate
    int n_rph_sum_rep1 = n_rph_sum_per_rep[i];
    array[2, n_rph_sum_rep1] int index_dedup1_rep1 = index_dedup1_mcmc[, index_rph_start:(index_rph_start + n_rph_sum_rep1 - 1)];
    matrix[n_rph_sum_rep1, 2] cv_square_rep1 = cv_square_mcmc[index_rph_start:(index_rph_start + n_rph_sum_rep1 - 1), ];
    vector[n_rph_sum_rep1] oph_rep1 = oph_mcmc[index_rph_start:(index_rph_start + n_rph_sum_rep1 - 1)];

    // Calculate relative peak heights
    matrix[n_rph_sum_rep1, 2] rph_dedup1 = calc_rph_dedup1(deg, n_rph_sum_rep1, rph_base_use, index_dedup1_rep1, deg_size_mcmc[index_rph_start:(index_rph_start + n_rph_sum_rep1 - 1)]);

    // Calculate gamma parameters
    matrix[n_rph_sum_rep1, 2] gamma_param = calc_gamma_param(rph_dedup1, cv_square_rep1, 2, n_rph_sum_rep1, tph[i]);

    // Calculate log-probability densities
    int n_dedup1_sum_rep1 = n_dedup1_sum_per_rep[i];
    vector[n_dedup1_sum_rep1] lpd_dedup1_use_rep1 = rep_vector(0, n_dedup1_sum_rep1);

    // Log probability density for general peaks
    int n_general_sum_rep1 = n_general_sum_per_rep[i];
    if(n_general_sum_rep1 != 0){
      lpd_dedup1_use_rep1[1:n_general_sum_rep1] = calc_lpd_general(n_general_sum_rep1,
                                                                   oph_rep1[1:n_general_sum_rep1],
                                                                   gamma_param[1:n_general_sum_rep1, 1],
                                                                   gamma_param[1:n_general_sum_rep1, 2]);
    }

    // Log probability density for drop-out peaks
    int n_dropout_sum_rep1 = n_dropout_sum_per_rep[i];
    if(n_dropout_sum_rep1 != 0){
      lpd_dedup1_use_rep1[(n_general_sum_rep1 + 1):(n_general_sum_rep1 + n_dropout_sum_rep1)] = calc_lpd_dropout(n_dropout_sum_rep1,
                                                                                                                 at_dropout[index_dropout_start:(index_dropout_start + n_dropout_sum_rep1 - 1)],
                                                                                                                 gamma_param[(n_general_sum_rep1 + 1):(n_general_sum_rep1 + n_dropout_sum_rep1), 1],
                                                                                                                 gamma_param[(n_general_sum_rep1 + 1):(n_general_sum_rep1 + n_dropout_sum_rep1), 2]);
    }

    // Log probability density for drop-in peaks
    int n_dropin_sum_rep1 = n_dropin_sum_per_rep[i];
    if(n_dropin_sum_rep1 != 0){
      vector[n_dropin_sum_rep1] loglike_di_rep1 = loglike_di[index_dropin_start:(index_dropin_start + n_dropin_sum_rep1 - 1)];
      lpd_dedup1_use_rep1[(n_general_sum_rep1 + n_dropout_sum_rep1 + 1):(n_general_sum_rep1 + n_dropout_sum_rep1 + n_dropin_sum_rep1)] = loglike_di_rep1;
    }

    // Assign each log probability density to each peak for the genotype combination
    int n_dedup2_sum_rep1 = n_dedup2_sum_per_rep[i];
    array[n_dedup2_sum_rep1] int index_dedup2_use_rep1 = index_dedup2_use[index_dedup2_start:(index_dedup2_start + n_dedup2_sum_rep1 - 1)];
    vector[n_dedup2_sum_rep1] lpd_dedup2_use_rep1 = lpd_dedup1_use_rep1[index_dedup2_use_rep1];

    // Calculate overall likelihood
    int index_lpd_start = 1;
    int index_gtc_start = 1;
    for(j in 1:n_mk){
      int n_peak = n_peak_use[j];
      int n_gtc = n_gtc_use[j];
      int n_dedup2 = n_dedup2_use[j];
      row_vector[n_gtc] loglike_gtc = calc_loglike_gtc(n_peak,
                                                       n_gtc,
                                                       to_matrix(lpd_dedup2_use_rep1[index_lpd_start:(index_lpd_start + n_dedup2 - 1)], n_peak, n_gtc),
                                                       gtc_logprob_use[index_gtc_start:(index_gtc_start + n_gtc - 1)]);

      target += log_sum_exp(loglike_gtc);

      index_lpd_start += n_dedup2;
      index_gtc_start += n_gtc;
    }

    // Update indices for next replicate
    index_rph_start += n_rph_sum_rep1;
    index_dropout_start += n_dropout_sum_rep1;
    index_dropin_start += n_dropin_sum_rep1;
    index_dedup2_start += n_dedup2_sum_rep1;
  }
}
