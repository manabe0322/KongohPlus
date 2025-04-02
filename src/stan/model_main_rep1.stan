functions{
  matrix calc_rph_dedup1(vector mp, vector deg, int hnc, int n_rph_sum_rep1, vector rph_base_use, array[,] int index_dedup1_rep1, vector size_rep1, real min_size){
    matrix[n_rph_sum_rep1, 2 * hnc] rph_dedup1;
    for(i in 1:hnc){
      real mp_c1 = mp[i];
      real deg_c1 = deg[i];
      vector[n_rph_sum_rep1] rph_param = mp_c1 * deg_c1^((size_rep1 - min_size) / min_size);
      for(j in 1:2){
        rph_dedup1[, 2 * i - 2 + j] = rph_base_use[index_dedup1_rep1[2 * i - 2 + j, ]] .* rph_param;
      }
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
  int hnc;
  int n_mk;
  real tph_mean;
  real tph_sd;

  int n_rph_base_sum;
  vector[n_rph_base_sum] rph_base_use;
  real min_size;

  int n_dedup1_sum;

  int n_rph_sum;
  array[2 * hnc, n_rph_sum] int index_dedup1_rep1;
  matrix[n_rph_sum, 2 * hnc] cv_square_rep1;
  vector[n_rph_sum] size_rep1;
  vector[n_rph_sum] oph_rep1;

  int n_general_sum;
  int n_dropout_sum;
  vector[n_dropout_sum] at_dropout_rep1;
  int n_dropin_sum;
  vector[n_dropin_sum] loglike_di_rep1;

  int n_dedup2_sum;
  array[n_mk] int n_dedup2_use;
  array[n_dedup2_sum] int index_dedup2_use_rep1;

  array[n_mk] int n_peak_use;
  array[n_mk] int n_gtc_use;
  row_vector[sum(n_gtc_use)] gtc_logprob_use;
}

parameters{
  real<lower = 0> tph;
  ordered[hnc] lambda;
  vector<lower = 0, upper = 1>[hnc] deg;
}

transformed parameters{
  simplex[hnc] mp;
  mp = softmax(lambda);
}

model{
  tph ~ normal(tph_mean, tph_sd);
  lambda ~ normal(0, 1);
  deg ~ beta(1, 1);

  // Calculate relative peak heights
  matrix[n_rph_sum, 2 * hnc] rph_dedup1 = calc_rph_dedup1(mp, deg, hnc, n_rph_sum, rph_base_use, index_dedup1_rep1, size_rep1, min_size);

  // Calculate gamma parameters
  matrix[n_rph_sum, 2] gamma_param = calc_gamma_param(rph_dedup1, cv_square_rep1, 2 * hnc, n_rph_sum, tph);

  // Calculate log-probability densities
  vector[n_dedup1_sum] lpd_dedup1_use_rep1 = rep_vector(0, n_dedup1_sum);

  // For general peaks
  if(n_general_sum != 0){
    vector[n_general_sum] oph_rep1_general = oph_rep1[1:n_general_sum];
    vector[n_general_sum] gamma_alpha = gamma_param[1:n_general_sum, 1];
    vector[n_general_sum] gamma_beta = gamma_param[1:n_general_sum, 2];
    lpd_dedup1_use_rep1[1:n_general_sum] = calc_lpd_general(n_general_sum, oph_rep1_general, gamma_alpha, gamma_beta);
  }

  // For drop-out peaks
  if(n_dropout_sum != 0){
    vector[n_dropout_sum] gamma_alpha = gamma_param[(n_general_sum + 1):(n_general_sum + n_dropout_sum), 1];
    vector[n_dropout_sum] gamma_beta = gamma_param[(n_general_sum + 1):(n_general_sum + n_dropout_sum), 2];
    lpd_dedup1_use_rep1[(n_general_sum + 1):(n_general_sum + n_dropout_sum)] = calc_lpd_dropout(n_dropout_sum, at_dropout_rep1, gamma_alpha, gamma_beta);
  }

  // For drop-in peaks
  if(n_dropin_sum != 0){
    lpd_dedup1_use_rep1[(n_general_sum + n_dropout_sum + 1):(n_general_sum + n_dropout_sum + n_dropin_sum)] = loglike_di_rep1;
  }

  vector[n_dedup2_sum] lpd_dedup2_use_rep1 = lpd_dedup1_use_rep1[index_dedup2_use_rep1];

  int index_lpd_start = 1;
  int index_gtc_start = 1;
  for(j in 1:n_mk){
    int n_peak = n_peak_use[j];
    int n_gtc = n_gtc_use[j];
    int n_dedup2 = n_dedup2_use[j];
    matrix[n_peak, n_gtc] lpd_dedup2 = to_matrix(lpd_dedup2_use_rep1[index_lpd_start:(index_lpd_start + n_dedup2 - 1)], n_peak, n_gtc);
    row_vector[n_gtc] gtc_logprob_mk = gtc_logprob_use[index_gtc_start:(index_gtc_start + n_gtc - 1)];
    row_vector[n_gtc] loglike_gtc = calc_loglike_gtc(n_peak, n_gtc, lpd_dedup2, gtc_logprob_mk);

    target += log_sum_exp(loglike_gtc);
    index_lpd_start += n_dedup2;
    index_gtc_start += n_gtc;
  }
}
