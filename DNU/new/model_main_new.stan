functions{
  matrix calc_rph_dedup1(real tph, vector mp, vector deg, int hnc, int n_dedup1_sum, vector rph_base_use, array[,] int index_dedup1, vector size_dedup1, real mean_size, vector ae_dedup1){
    matrix[n_dedup1_sum, 2 * hnc] rph_dedup1;
    for(i in 1:hnc){
      vector[n_dedup1_sum] rphs_1 = rph_base_use[index_dedup1[2 * i - 1, ]];
      vector[n_dedup1_sum] rphs_2 = rph_base_use[index_dedup1[2 * i, ]];
      real mp_c1 = mp[i];
      real deg_c1 = deg[i];

      rph_dedup1[, 2 * i - 1] = rphs_1 * tph * mp_c1 .* deg_c1^((size_dedup1 - mean_size) / mean_size) .* ae_dedup1;
      rph_dedup1[, 2 * i] = rphs_2 * tph * mp_c1 .* deg_c1^((size_dedup1 - mean_size) / mean_size) .* ae_dedup1;
    }
    return rph_dedup1;
  }

  real calc_sum2(row_vector rphs, real cv_al_sq, real cv_st_sq, array[] int index_1, int index_al_max){
    int n_al = size(rphs);
    real sum2 = 0;
    for(i in 1:n_al){
      real rph = rphs[i];
      if(rph > 0){
        if(index_1[i] <= index_al_max){
          sum2 += cv_al_sq * rph;
        }else{
          sum2 += cv_st_sq * rph;
        }
      }
    }
    return sum2;
  }

  vector calc_lpd_general(real cv_al_sq, real cv_st_sq, int n_mk, int n_general_sum, array[] int n_general_use, array[] int n_peak_use, array[] int n_rph_base_use, vector oph_dedup1_general, matrix rph_dedup1_general, array[,] int index_dedup1_general){
    vector[n_general_sum] lpd_general = rep_vector(0, n_general_sum);
    int n_al = cols(rph_dedup1_general);
    int index_start = 0;
    int index_al_max = 0;
    for(i in 1:n_mk){
      int n_gen = n_general_use[i];
      int n_peak = n_peak_use[i];
      index_al_max += n_peak;
      if(n_gen != 0){
        for(j in 1:n_gen){
          array[n_al] int index_1 = index_dedup1_general[, index_start + j];
          row_vector[n_al] rphs = rph_dedup1_general[index_start + j, ];
          real sum1 = sum(rphs);
          real sum2 = calc_sum2(rphs, cv_al_sq, cv_st_sq, index_1, index_al_max);
          real alpha = sum1 * sum1 / sum2;
          real beta = sum1 / sum2;
          lpd_general[index_start + j] = gamma_lpdf(oph_dedup1_general[index_start + j] | alpha, beta);
        }
        index_start += n_gen;
      }
      index_al_max -= n_peak;
      index_al_max += n_rph_base_use[i];
    }
    return lpd_general;
  }

  vector calc_lpd_dropout(real cv_al_sq, real cv_st_sq, int n_mk, int n_dropout_sum, array[] int n_dropout_use, array[] int n_peak_use, array[] int n_rph_base_use, vector ats_use, matrix rph_dedup1_dropout, array[,] int index_dedup1_dropout){
    vector[n_dropout_sum] lpd_dropout = rep_vector(0, n_dropout_sum);
    int n_al = cols(rph_dedup1_dropout);
    int index_start = 0;
    int index_al_max = 0;
    for(i in 1:n_mk){
      int n_do = n_dropout_use[i];
      int n_peak = n_peak_use[i];
      index_al_max += n_peak;
      real at = ats_use[i];
      if(n_do != 0){
        for(j in 1:n_do){
          array[n_al] int index_1 = index_dedup1_dropout[, index_start + j];
          row_vector[n_al] rphs = rph_dedup1_dropout[index_start + j, ];
          real sum1 = sum(rphs);
          real sum2 = calc_sum2(rphs, cv_al_sq, cv_st_sq, index_1, index_al_max);
          real alpha = sum1 * sum1 / sum2;
          real beta = sum1 / sum2;
          lpd_dropout[index_start + j] = gamma_lcdf(at | alpha, beta);
        }
        index_start += n_do;
      }
      index_al_max -= n_peak;
      index_al_max += n_rph_base_use[i];
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
  array[n_mk] int n_rph_base_use;
  vector[n_rph_base_sum] rph_base_use;
  vector[n_mk] ae_alpha_use;
  vector[n_mk] ats_use;
  real mean_size;
  real cv_allele;
  real cv_stutter;

  int n_dedup1_sum;

  int n_general_sum;
  array[n_mk] int n_general_use;
  array[2 * hnc, n_general_sum] int index_dedup1_general;
  vector[n_general_sum] size_dedup1_general;
  vector[n_general_sum] oph_dedup1_general;
  array[n_general_sum] int id_mk_dedup1_general;

  int n_dropout_sum;
  array[n_mk] int n_dropout_use;
  array[2 * hnc, n_dropout_sum] int index_dedup1_dropout;
  vector[n_dropout_sum] size_dedup1_dropout;
  array[n_dropout_sum] int id_mk_dedup1_dropout;

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
  real cv_al_sq = pow(cv_allele, 2);
  real cv_st_sq = pow(cv_stutter, 2);
}

parameters{
  real<lower = ltph, upper = htph> tph;
  ordered[hnc] lambda;
  vector<lower = 0, upper = 1>[hnc] deg;
  simplex[n_mk] pae;
}

transformed parameters{
  simplex[hnc] mp;
  mp = softmax(lambda);
}

model{
  tph ~ normal(tph_mean, tph_sd);
  lambda ~ normal(0, 10);
  pae ~ dirichlet(ae_alpha_use);

  vector[n_dedup1_sum] lpd_dedup1_use = rep_vector(0, n_dedup1_sum);

  // general
  if(n_general_sum != 0){
    vector[n_general_sum] ae_dedup1_general = n_mk * pae[id_mk_dedup1_general];
    matrix[n_general_sum, 2 * hnc] rph_dedup1_general = calc_rph_dedup1(tph, mp, deg, hnc, n_general_sum, rph_base_use, index_dedup1_general, size_dedup1_general, mean_size, ae_dedup1_general);
    lpd_dedup1_use[1:n_general_sum] = calc_lpd_general(cv_al_sq, cv_st_sq, n_mk, n_general_sum, n_general_use, n_peak_use, n_rph_base_use, oph_dedup1_general, rph_dedup1_general, index_dedup1_general);
  }

  // drop-out
  if(n_dropout_sum != 0){
    vector[n_dropout_sum] ae_dedup1_dropout = n_mk * pae[id_mk_dedup1_dropout];
    matrix[n_dropout_sum, 2 * hnc] rph_dedup1_dropout = calc_rph_dedup1(tph, mp, deg, hnc, n_dropout_sum, rph_base_use, index_dedup1_dropout, size_dedup1_dropout, mean_size, ae_dedup1_dropout);
    lpd_dedup1_use[(n_general_sum + 1):(n_general_sum + n_dropout_sum)] = calc_lpd_dropout(cv_al_sq, cv_st_sq, n_mk, n_dropout_sum, n_dropout_use, n_peak_use, n_rph_base_use, ats_use, rph_dedup1_dropout, index_dedup1_dropout);
  }

  // drop-in
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
