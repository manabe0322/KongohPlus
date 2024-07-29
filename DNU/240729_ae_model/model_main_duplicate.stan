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
  deg ~ normal(-0.025, 0.0125);
  pae ~ dirichlet(tph * ae_alpha_use);
}
