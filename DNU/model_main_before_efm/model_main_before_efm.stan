functions{
  vector calc_eph_dedup1(real tph, vector mp, vector deg, int hnc, int n_dedup1_all, vector eph_base, array[,] int index_dedup1_mcmc, vector size_dedup1_mcmc){
    vector[n_dedup1_all] eph_dedup1;
    eph_dedup1 = rep_vector(0, n_dedup1_all);
    for(i in 1:hnc){
      vector[n_dedup1_all] eph_base_c1 = eph_base[index_dedup1_mcmc[2 * i - 1, ]] + eph_base[index_dedup1_mcmc[2 * i, ]];
      real tph_c1 = tph * mp[i];
      real deg_c1 = deg[i];
      eph_dedup1 += eph_base_c1 * tph_c1 .* exp(deg_c1 * size_dedup1_mcmc);
    }
    return eph_dedup1;
  }

  real calc_lpd(real oph, real eph, real at, real beta, real di_alpha, real di_beta, real logprob_di, real di_min){
    real lpd = 0;
    if(eph > 0){
      if(oph >= at){
        lpd = gamma_lpdf(oph | eph * beta, beta);
      }else{
        lpd = gamma_lcdf(at | eph * beta, beta);
      }
    }else if(eph < 0){
      lpd = gamma_lpdf(oph - di_min | di_alpha, di_beta) + logprob_di;
    }
    return lpd;
  }

  vector calc_lpd_mk(int n_dedup_eph, vector oph_mk, vector eph_mk, real at, real beta, real di_alpha, real di_beta, real logprob_di, real di_min){
    vector[n_dedup_eph] lpd_mk;
    for(i in 1:n_dedup_eph){
      lpd_mk[i] = calc_lpd(oph_mk[i], eph_mk[i], at, beta, di_alpha, di_beta, logprob_di, di_min);
    }
    return lpd_mk;
  }

  vector calc_lpd_gtc(int n_gtc, matrix lpd_extract, vector gtc_logprob_mk){
    vector[n_gtc] lpd_gtc;
    for(i in 1:n_gtc){
      lpd_gtc[i] = sum(lpd_extract[, i]) + gtc_logprob_mk[i];
    }
    return lpd_gtc;
  }
}

data{
  int Hnc;
  int N_Mk;
  real Ltph;
  real Htph;

  int N_Eph_Base_All;
  vector[N_Eph_Base_All] Eph_Base;
  vector[N_Mk] AE_Alpha_All_Mk;

  int N_Dedup1_All;
  array[N_Mk] int N_Dedup1_Mcmc;
  array[2 * Hnc, N_Dedup1_All] int Index_Dedup1_Mcmc;
  vector[N_Dedup1_All] Size_Dedup1_Mcmc;
  vector[N_Dedup1_All] Oph_Dedup1_Mcmc;

  int N_Dedup2_All;
  array[N_Mk] int N_Dedup2_Mcmc;
  array[N_Dedup2_All] int Index_Dedup2_Mcmc;

  array[N_Mk] int N_Peak_Adopt;
  array[N_Mk] int N_Gtc_Mcmc;
  vector[sum(N_Gtc_Mcmc)] Gtc_Logprob_Mcmc;

  vector[N_Mk] At_All_Mk;
  real Beta;
  real Di_Alpha;
  real Di_Beta;
  real Logprob_Di;
  real Di_Min;
}

transformed data{
  real tph_mean = (Ltph + Htph) / 2;
  real tph_sd = (Htph - Ltph) / 2;
}

parameters{
  real<lower = Ltph, upper = Htph> tph;
  positive_ordered[Hnc] lambda;
  vector<lower = -0.05, upper = 0>[Hnc] deg;
  simplex[N_Mk] pae;
}

transformed parameters{
  simplex[Hnc] mp = lambda / sum(lambda);
}

model{
  tph ~ normal(tph_mean, tph_sd);
  target += gamma_lpdf(lambda | rep_vector(1, Hnc), 1);
  deg ~ normal(-0.025, 0.0125);
  pae ~ dirichlet(tph * AE_Alpha_All_Mk);

  vector[N_Dedup1_All] eph_dedup1 = calc_eph_dedup1(tph, mp, deg, Hnc, N_Dedup1_All, Eph_Base, Index_Dedup1_Mcmc, Size_Dedup1_Mcmc);

  int eph_start = 1;
  int dedup2_start = 1;
  int gtc_start = 1;
  for(i in 1:N_Mk){
    int n_dedup_eph = N_Dedup1_Mcmc[i];
    vector[n_dedup_eph] eph_mk = N_Mk * pae[i] * eph_dedup1[eph_start:(eph_start + n_dedup_eph - 1)];
    vector[n_dedup_eph] oph_mk = Oph_Dedup1_Mcmc[eph_start:(eph_start + n_dedup_eph - 1)];
    real at = At_All_Mk[i];

    vector[n_dedup_eph] lpd_mk = calc_lpd_mk(n_dedup_eph, oph_mk, eph_mk, at, Beta, Di_Alpha, Di_Beta, Logprob_Di, Di_Min);

    int n_dedup2 = N_Dedup2_Mcmc[i];
    int n_gtc = N_Gtc_Mcmc[i];
    int n_peak = N_Peak_Adopt[i];
    array[n_dedup2] int index_dedup2_mk = Index_Dedup2_Mcmc[dedup2_start:(dedup2_start + n_dedup2 - 1)];
    matrix[n_peak, n_gtc] lpd_extract = to_matrix(lpd_mk[index_dedup2_mk], n_peak, n_gtc);
    vector[n_gtc] gtc_logprob_mk = Gtc_Logprob_Mcmc[gtc_start:(gtc_start + n_gtc - 1)];
    vector[n_gtc] lpd_gtc = calc_lpd_gtc(n_gtc, lpd_extract, gtc_logprob_mk);

    target += log_sum_exp(lpd_gtc);

    eph_start += n_dedup_eph;
    dedup2_start += n_dedup2;
    gtc_start += n_gtc;
  }
}
