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

  row_vector calc_loglike_gtc(int n_peak, int n_gtc, matrix lpd_dedup2_1st, matrix lpd_dedup2_2nd, row_vector gtc_logprob_mk){
    row_vector[n_peak] rv = rep_row_vector(1, n_peak);
    row_vector[n_gtc] loglike_gtc = rv * lpd_dedup2_1st + rv * lpd_dedup2_2nd + gtc_logprob_mk;
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

  int n_dedup1_sum_1st;
  array[2 * hnc, n_dedup1_sum_1st] int index_dedup1_use_1st;
  vector[n_dedup1_sum_1st] size_dedup1_use_1st;
  vector[n_dedup1_sum_1st] oph_dedup1_use_1st;
  array[n_dedup1_sum_1st] int id_mk_dedup1_1st;

  int n_general_sum_1st;
  int n_dropout_sum_1st;
  vector[n_dropout_sum_1st] at_dropout_1st;
  int n_dropin_sum_1st;
  vector[n_dropin_sum_1st] loglike_di_1st;

  int n_dedup1_sum_2nd;
  array[2 * hnc, n_dedup1_sum_2nd] int index_dedup1_use_2nd;
  vector[n_dedup1_sum_2nd] size_dedup1_use_2nd;
  vector[n_dedup1_sum_2nd] oph_dedup1_use_2nd;
  array[n_dedup1_sum_2nd] int id_mk_dedup1_2nd;

  int n_general_sum_2nd;
  int n_dropout_sum_2nd;
  vector[n_dropout_sum_2nd] at_dropout_2nd;
  int n_dropin_sum_2nd;
  vector[n_dropin_sum_2nd] loglike_di_2nd;

  int n_dedup2_sum;
  array[n_mk] int n_dedup2_use;
  array[n_dedup2_sum] int index_dedup2_use_1st;
  array[n_dedup2_sum] int index_dedup2_use_2nd;

  array[n_mk] int n_peak_use;
  array[n_mk] int n_gtc_use;
  row_vector[sum(n_gtc_use)] gtc_logprob_use;
}

transformed data{
  real tph_mean_1st = (ltph + htph) / 2;
  real tph_mean_2nd = (ltph + htph) / 2;
  real tph_sd_1st = (htph - ltph) / 2;
  real tph_sd_2nd = (htph - ltph) / 2;
  real cv_m2 = pow(cv, -2);
  real cv_p2 = pow(cv, 2);
}

parameters{
  real<lower = ltph, upper = htph> tph_1st;
  real<lower = ltph, upper = htph> tph_2nd;
  positive_ordered[hnc] lambda;
  vector<lower = -0.05, upper = 0>[hnc] deg;
  simplex[n_mk] pae;
}

transformed parameters{
  simplex[hnc] mp = lambda / sum(lambda);
}

model{
  tph_1st ~ normal(tph_mean_1st, tph_sd_1st);
  tph_2nd ~ normal(tph_mean_2nd, tph_sd_2nd);
  target += gamma_lpdf(lambda | rep_vector(1, hnc), 1);
  deg ~ normal(0, 0.015);
  pae ~ dirichlet(ae_alpha_use);

  vector[n_dedup1_sum_1st] ae_dedup1_use_1st = n_mk * pae[id_mk_dedup1_1st];
  vector[n_dedup1_sum_1st] rph_dedup1_use_1st = calc_rph_dedup1(mp, deg, hnc, n_dedup1_sum_1st, rph_base_use, index_dedup1_use_1st, size_dedup1_use_1st, ae_dedup1_use_1st);
  vector[n_dedup1_sum_1st] lpd_dedup1_use_1st = rep_vector(0, n_dedup1_sum_1st);
  if(n_general_sum_1st != 0){
    vector[n_general_sum_1st] oph_general_1st = oph_dedup1_use_1st[1:n_general_sum_1st];
    vector[n_general_sum_1st] rph_general_1st = rph_dedup1_use_1st[1:n_general_sum_1st];
    lpd_dedup1_use_1st[1:n_general_sum_1st] = calc_lpd_general(n_general_sum_1st, oph_general_1st, rph_general_1st, cv_m2, cv_p2, tph_1st);
  }
  if(n_dropout_sum_1st != 0){
    vector[n_dropout_sum_1st] rph_dropout_1st = rph_dedup1_use_1st[(n_general_sum_1st + 1):(n_general_sum_1st + n_dropout_sum_1st)];
    lpd_dedup1_use_1st[(n_general_sum_1st + 1):(n_general_sum_1st + n_dropout_sum_1st)] = calc_lpd_dropout(n_dropout_sum_1st, at_dropout_1st, rph_dropout_1st, cv_m2, cv_p2, tph_1st);
  }
  if(n_dropin_sum_1st != 0){
    lpd_dedup1_use_1st[(n_general_sum_1st + n_dropout_sum_1st + 1):(n_general_sum_1st + n_dropout_sum_1st + n_dropin_sum_1st)] = loglike_di_1st;
  }
  vector[n_dedup2_sum] lpd_dedup2_use_1st = lpd_dedup1_use_1st[index_dedup2_use_1st];

  vector[n_dedup1_sum_2nd] ae_dedup1_use_2nd = n_mk * pae[id_mk_dedup1_2nd];
  vector[n_dedup1_sum_2nd] rph_dedup1_use_2nd = calc_rph_dedup1(mp, deg, hnc, n_dedup1_sum_2nd, rph_base_use, index_dedup1_use_2nd, size_dedup1_use_2nd, ae_dedup1_use_2nd);
  vector[n_dedup1_sum_2nd] lpd_dedup1_use_2nd = rep_vector(0, n_dedup1_sum_2nd);
  if(n_general_sum_2nd != 0){
    vector[n_general_sum_2nd] oph_general_2nd = oph_dedup1_use_2nd[1:n_general_sum_2nd];
    vector[n_general_sum_2nd] rph_general_2nd = rph_dedup1_use_2nd[1:n_general_sum_2nd];
    lpd_dedup1_use_2nd[1:n_general_sum_2nd] = calc_lpd_general(n_general_sum_2nd, oph_general_2nd, rph_general_2nd, cv_m2, cv_p2, tph_2nd);
  }
  if(n_dropout_sum_2nd != 0){
    vector[n_dropout_sum_2nd] rph_dropout_2nd = rph_dedup1_use_2nd[(n_general_sum_2nd + 1):(n_general_sum_2nd + n_dropout_sum_2nd)];
    lpd_dedup1_use_2nd[(n_general_sum_2nd + 1):(n_general_sum_2nd + n_dropout_sum_2nd)] = calc_lpd_dropout(n_dropout_sum_2nd, at_dropout_2nd, rph_dropout_2nd, cv_m2, cv_p2, tph_2nd);
  }
  if(n_dropin_sum_2nd != 0){
    lpd_dedup1_use_2nd[(n_general_sum_2nd + n_dropout_sum_2nd + 1):(n_general_sum_2nd + n_dropout_sum_2nd + n_dropin_sum_2nd)] = loglike_di_2nd;
  }
  vector[n_dedup2_sum] lpd_dedup2_use_2nd = lpd_dedup1_use_2nd[index_dedup2_use_2nd];

  int dedup2_start = 1;
  int gtc_start = 1;
  for(i in 1:n_mk){
    int n_peak = n_peak_use[i];
    int n_gtc = n_gtc_use[i];
    int n_dedup2 = n_dedup2_use[i];
    matrix[n_peak, n_gtc] lpd_dedup2_1st = to_matrix(lpd_dedup2_use_1st[dedup2_start:(dedup2_start + n_dedup2 - 1)], n_peak, n_gtc);
    matrix[n_peak, n_gtc] lpd_dedup2_2nd = to_matrix(lpd_dedup2_use_2nd[dedup2_start:(dedup2_start + n_dedup2 - 1)], n_peak, n_gtc);
    row_vector[n_gtc] gtc_logprob_mk = gtc_logprob_use[gtc_start:(gtc_start + n_gtc - 1)];
    row_vector[n_gtc] loglike_gtc = calc_loglike_gtc(n_peak, n_gtc, lpd_dedup2_1st, lpd_dedup2_2nd, gtc_logprob_mk);

    target += log_sum_exp(loglike_gtc);

    dedup2_start += n_dedup2;
    gtc_start += n_gtc;
  }
}
