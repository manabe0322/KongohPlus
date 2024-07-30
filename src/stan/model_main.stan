functions{
  vector calc_rph_dedup1(vector mp, vector deg, int hnc, int n_dedup1_sum, vector rph_base_use, array[,] int index_dedup1_use, vector size_dedup1_use, vector ae_dedup1_use){
    vector[n_dedup1_sum] rph_dedup1 = rep_vector(0, n_dedup1_sum);
    for(i in 1:hnc){
      vector[n_dedup1_sum] rph_base_use_c1 = rph_base_use[index_dedup1_use[2 * i - 1, ]] + rph_base_use[index_dedup1_use[2 * i, ]];
      real mp_c1 = mp[i];
      real deg_c1 = deg[i];
      rph_dedup1 += rph_base_use_c1 * mp_c1 .* exp(deg_c1 * size_dedup1_use) .* ae_dedup1_use;
    }
    return rph_dedup1;
  }

  vector calc_lpd_general(int n_general_sum, vector oph_general, vector rph_general, real cv_m2, real cv_p2, real tph){
    vector[n_general_sum] lpd_general;
    vector[n_general_sum] gamma_alpha = cv_m2 * rph_general;
    real gamma_beta = 1 / (tph * cv_p2);
    for(i in 1:n_general_sum){
      lpd_general[i] = gamma_lpdf(oph_general[i] | gamma_alpha[i], gamma_beta);
    }
    return lpd_general;
  }

  vector calc_lpd_dropout(int n_dropout_sum, vector at_dropout, vector rph_dropout, real cv_m2, real cv_p2, real tph){
    vector[n_dropout_sum] lpd_dropout;
    vector[n_dropout_sum] gamma_alpha = cv_m2 * rph_dropout;
    real gamma_beta = 1 / (tph * cv_p2);
    for(i in 1:n_dropout_sum){
      lpd_dropout[i] = gamma_lcdf(at_dropout[i] | gamma_alpha[i], gamma_beta);
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
  real ltph;
  real htph;

  int n_rph_base_sum;
  vector[n_rph_base_sum] rph_base_use;
  vector[n_mk] ae_alpha_use;
  real cv;

  int n_dedup1_sum;
  array[2 * hnc, n_dedup1_sum] int index_dedup1_use;
  vector[n_dedup1_sum] size_dedup1_use;
  vector[n_dedup1_sum] oph_dedup1_use;
  array[n_dedup1_sum] int id_mk_dedup1;

  int n_general_sum;
  int n_dropout_sum;
  vector[n_dropout_sum] at_dropout;
  int n_dropin_sum;
  vector[n_dropin_sum] loglike_di;

  int n_dedup2_sum;
  array[n_mk] int n_dedup2_use;
  array[n_dedup2_sum] int index_dedup2_use;

  array[n_mk] int n_peak_use;
  array[n_mk] int n_gtc_use;
  row_vector[sum(n_gtc_use)] gtc_logprob_use;
}

transformed data{
  real tph_mean = (ltph + htph) / 2;
  real tph_sd = (htph - ltph) / 2;
  real cv_m2 = pow(cv, -2);
  real cv_p2 = pow(cv, 2);
}

parameters{
  real<lower = ltph, upper = htph> tph;
  ordered[hnc] lambda;
  vector<lower = -0.05, upper = 0>[hnc] deg;
  simplex[n_mk] pae;
}

transformed parameters{
  simplex[hnc] mp; // = lambda / sum(lambda);
  mp = softmax(lambda);
}

model{
  lambda ~ normal(0, 10);
  //target += gamma_lpdf(lambda | 1, 1);
  pae ~ dirichlet(rep_vector(1, n_mk));

  vector[n_dedup1_sum] ae_dedup1_use = n_mk * pae[id_mk_dedup1];

  vector[n_dedup1_sum] rph_dedup1_use = calc_rph_dedup1(mp, deg, hnc, n_dedup1_sum, rph_base_use, index_dedup1_use, size_dedup1_use, ae_dedup1_use);

  vector[n_dedup1_sum] lpd_dedup1_use = rep_vector(0, n_dedup1_sum);

  if(n_general_sum != 0){
    vector[n_general_sum] oph_general = oph_dedup1_use[1:n_general_sum];
    vector[n_general_sum] rph_general = rph_dedup1_use[1:n_general_sum];
    lpd_dedup1_use[1:n_general_sum] = calc_lpd_general(n_general_sum, oph_general, rph_general, cv_m2, cv_p2, tph);
  }

  if(n_dropout_sum != 0){
    vector[n_dropout_sum] rph_dropout = rph_dedup1_use[(n_general_sum + 1):(n_general_sum + n_dropout_sum)];
    lpd_dedup1_use[(n_general_sum + 1):(n_general_sum + n_dropout_sum)] = calc_lpd_dropout(n_dropout_sum, at_dropout, rph_dropout, cv_m2, cv_p2, tph);
  }

  if(n_dropin_sum != 0){
    lpd_dedup1_use[(n_general_sum + n_dropout_sum + 1):(n_general_sum + n_dropout_sum + n_dropin_sum)] = loglike_di;
  }

  vector[n_dedup2_sum] lpd_dedup2_use = lpd_dedup1_use[index_dedup2_use];
  int dedup2_start = 1;
  int gtc_start = 1;
  for(i in 1:n_mk){
    int n_peak = n_peak_use[i];
    int n_gtc = n_gtc_use[i];
    int n_dedup2 = n_dedup2_use[i];
    matrix[n_peak, n_gtc] lpd_dedup2 = to_matrix(lpd_dedup2_use[dedup2_start:(dedup2_start + n_dedup2 - 1)], n_peak, n_gtc);
    row_vector[n_gtc] gtc_logprob_mk = gtc_logprob_use[gtc_start:(gtc_start + n_gtc - 1)];
    row_vector[n_gtc] loglike_gtc = calc_loglike_gtc(n_peak, n_gtc, lpd_dedup2, gtc_logprob_mk);

    target += log_sum_exp(loglike_gtc);

    dedup2_start += n_dedup2;
    gtc_start += n_gtc;
  }
}
