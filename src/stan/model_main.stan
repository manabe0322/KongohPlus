functions{
  matrix calc_rph_dedup1(vector mp, vector deg, int hnc, int n_dedup1_sum_rep1, vector rph_base_use, array[,] int index_dedup1_use_rep1, vector size_dedup1_use_rep1, real min_size, vector ae_dedup1_use_rep1){
    matrix[2 * hnc, n_dedup1_sum_rep1] rph_dedup1_use_rep1;
    for(i in 1:hnc){
      real mp_c1 = mp[i];
      real deg_c1 = deg[i];
      vector[n_dedup1_sum_rep1] rph_param = ae_dedup1_use_rep1 .* deg_c1^((size_dedup1_use_rep1 - min_size) / min_size) * mp_c1;
      rph_dedup1_use_rep1[2 * i - 1, ] = to_row_vector(rph_base_use[index_dedup1_use_rep1[2 * i - 1, ]] .* rph_param);
      rph_dedup1_use_rep1[2 * i, ] = to_row_vector(rph_base_use[index_dedup1_use_rep1[2 * i, ]] .* rph_param);
    }
    return rph_dedup1_use_rep1;
  }

  matrix calc_gamma_param(matrix rph_dedup1_use_rep1, matrix cv_square_dedup1_use, int n_roop, int n_dedup1_sum_rep1, real tph){
    matrix[n_roop, n_dedup1_sum_rep1] tmp1 = rph_dedup1_use_rep1 * tph;
    matrix[n_roop, n_dedup1_sum_rep1] tmp2 = rph_dedup1_use_rep1 .* cv_square_dedup1_use * tph;
    row_vector[n_dedup1_sum_rep1] tmp1_sum = to_row_vector(rep_vector(0, n_dedup1_sum_rep1));
    row_vector[n_dedup1_sum_rep1] tmp2_sum = to_row_vector(rep_vector(0, n_dedup1_sum_rep1));
    for(i in 1:n_roop){
      tmp1_sum += tmp1[i, ];
      tmp2_sum += tmp2[i, ];
    }

    matrix[2, n_dedup1_sum_rep1] gamma_param;
    gamma_param[1, ] = tmp1_sum^2 ./ tmp2_sum;
    gamma_param[2, ] = tmp1_sum ./ tmp2_sum;

    return gamma_param;
  }

  vector calc_lpd_general(int n_general_sum, vector oph_general_rep1, vector gamma_alpha, vector gamma_beta){
    vector[n_general_sum] lpd_general;
    for(i in 1:n_general_sum){
      lpd_general[i] = gamma_lpdf(oph_general_rep1[i] | gamma_alpha[i], gamma_beta[i]);
    }
    return lpd_general;
  }

  vector calc_lpd_dropout(int n_dropout_sum, vector at_dropout_rep1, vector gamma_alpha, vector gamma_beta){
    vector[n_dropout_sum] lpd_dropout;
    for(i in 1:n_dropout_sum){
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
  int n_rep;
  int n_mk;
  real ltph;
  real htph;

  int n_rph_base_sum;
  vector[n_rph_base_sum] rph_base_use;
  vector[n_mk] ae_alpha_use;
  real min_size;

  array[n_rep] int n_dedup1_sum_per_rep;
  int n_dedup1_sum;
  array[2 * hnc, n_dedup1_sum] int index_dedup1_use;
  matrix[2 * hnc, n_dedup1_sum] cv_square_dedup1_use;
  vector[n_dedup1_sum] size_dedup1_use;
  vector[n_dedup1_sum] oph_dedup1_use;
  array[n_dedup1_sum] int id_mk_dedup1;

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
  vector<lower = ltph, upper = htph>[n_rep] tph;
  ordered[hnc] lambda;
  vector<lower = 0, upper = 1>[hnc] deg;
  simplex[n_mk] pae;
}

transformed parameters{
  simplex[hnc] mp;
  mp = softmax(lambda);
}

model{
  tph ~ normal((ltph + htph) / 2, (htph - ltph) / 2);
  lambda ~ normal(0, 10);
  pae ~ dirichlet(ae_alpha_use);

  int index_dedup1_start = 1;
  int index_dropout_start = 1;
  int index_dropin_start = 1;
  int index_dedup2_start = 1;

  for(i in 1:n_rep){
    // Extract data for one replicate
    int n_dedup1_sum_rep1 = n_dedup1_sum_per_rep[i];
    array[2 * hnc, n_dedup1_sum_rep1] int index_dedup1_use_rep1 = index_dedup1_use[, index_dedup1_start:(index_dedup1_start + n_dedup1_sum_rep1 - 1)];
    vector[n_dedup1_sum_rep1] size_dedup1_use_rep1 = size_dedup1_use[index_dedup1_start:(index_dedup1_start + n_dedup1_sum_rep1 - 1)];
    array[n_dedup1_sum_rep1] int id_mk_dedup1_rep1 = id_mk_dedup1[index_dedup1_start:(index_dedup1_start + n_dedup1_sum_rep1 - 1)];
    vector[n_dedup1_sum_rep1] ae_dedup1_use_rep1 = n_mk * pae[id_mk_dedup1_rep1];
    vector[n_dedup1_sum_rep1] oph_dedup1_use_rep1 = oph_dedup1_use[index_dedup1_start:(index_dedup1_start + n_dedup1_sum_rep1 - 1)];
    int n_general_sum_rep1 = n_general_sum_per_rep[i];
    int n_dropout_sum_rep1 = n_dropout_sum_per_rep[i];
    int n_dropin_sum_rep1 = n_dropin_sum_per_rep[i];
    int n_dedup2_sum_rep1 = n_dedup2_sum_per_rep[i];
    array[n_dedup2_sum_rep1] int index_dedup2_use_rep1 = index_dedup2_use[index_dedup2_start:(index_dedup2_start + n_dedup2_sum_rep1 - 1)];

    matrix[2 * hnc, n_dedup1_sum_rep1] rph_dedup1_use_rep1 = calc_rph_dedup1(mp, deg, hnc, n_dedup1_sum_rep1, rph_base_use, index_dedup1_use_rep1, size_dedup1_use_rep1, min_size, ae_dedup1_use_rep1);
    matrix[2, n_dedup1_sum_rep1] gamma_param = calc_gamma_param(rph_dedup1_use_rep1, cv_square_dedup1_use, 2 * hnc, n_dedup1_sum_rep1, tph[i]);

    vector[n_dedup1_sum_rep1] lpd_dedup1_use_rep1 = rep_vector(0, n_dedup1_sum_rep1);

    if(n_general_sum_rep1 != 0){
      vector[n_general_sum_rep1] oph_general_rep1 = oph_dedup1_use_rep1[1:n_general_sum_rep1];
      vector[n_general_sum_rep1] gamma_alpha = to_vector(gamma_param[1, 1:n_general_sum_rep1]);
      vector[n_general_sum_rep1] gamma_beta = to_vector(gamma_param[2, 1:n_general_sum_rep1]);
      lpd_dedup1_use_rep1[1:n_general_sum_rep1] = calc_lpd_general(n_general_sum_rep1, oph_general_rep1, gamma_alpha, gamma_beta);
    }

    if(n_dropout_sum_rep1 != 0){
      vector[n_dropout_sum_rep1] at_dropout_rep1 = at_dropout[index_dropout_start:(index_dropout_start + n_dropout_sum_rep1 - 1)];
      vector[n_dropout_sum_rep1] gamma_alpha = to_vector(gamma_param[1, (n_general_sum_rep1 + 1):(n_general_sum_rep1 + n_dropout_sum_rep1)]);
      vector[n_dropout_sum_rep1] gamma_beta = to_vector(gamma_param[2, (n_general_sum_rep1 + 1):(n_general_sum_rep1 + n_dropout_sum_rep1)]);
      lpd_dedup1_use_rep1[(n_general_sum_rep1 + 1):(n_general_sum_rep1 + n_dropout_sum_rep1)] = calc_lpd_dropout(n_dropout_sum_rep1, at_dropout_rep1, gamma_alpha, gamma_beta);
    }

    if(n_dropin_sum_rep1 != 0){
      vector[n_dropin_sum_rep1] loglike_di_rep1 = loglike_di[index_dropin_start:(index_dropin_start + n_dropin_sum_rep1 - 1)];
      lpd_dedup1_use_rep1[(n_general_sum_rep1 + n_dropout_sum_rep1 + 1):(n_general_sum_rep1 + n_dropout_sum_rep1 + n_dropin_sum_rep1)] = loglike_di_rep1;
    }

    vector[n_dedup2_sum_rep1] lpd_dedup2_use_rep1 = lpd_dedup1_use_rep1[index_dedup2_use_rep1];

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

    index_dedup1_start += n_dedup1_sum_rep1;
    index_dropout_start += n_dropout_sum_rep1;
    index_dropin_start += n_dropin_sum_rep1;
    index_dedup2_start += n_dedup2_sum_rep1;
  }
}
