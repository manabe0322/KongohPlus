functions{
  vector calc_eph_dedup_1(real tph, vector mp, vector deg, int hnc, int n_dedup_1_all, vector eph_base, array[,] int index_dedup_1, vector size_dedup_1){
    vector[n_dedup_1_all] eph_dedup_1;
    eph_dedup_1 = rep_vector(0, n_dedup_1_all);
    for(i in 1:hnc){
      vector[n_dedup_1_all] eph_base_c1 = eph_base[index_dedup_1[2 * i - 1, ]] + eph_base[index_dedup_1[2 * i, ]];
      real tph_c1 = tph * mp[i];
      real deg_c1 = deg[i];
      eph_dedup_1 += eph_base_c1 * tph_c1 .* exp(deg_c1 * size_dedup_1);
    }
    return eph_dedup_1;
  }

  real calc_lpd(real oph, real eph, real at, real beta){
    real lpd = 0;
    if(eph != 0){
      if(oph >= at){
        lpd = gamma_lpdf(oph | eph * beta, beta);
      }else{
        lpd = gamma_lcdf(at | eph * beta, beta);
      }
    }
    return lpd;
  }

  vector calc_lpd_mk(int n_dedup_eph, vector oph_mk, vector eph_mk, real at, real beta){
    vector[n_dedup_eph] lpd_mk;
    for(i in 1:n_dedup_eph){
      lpd_mk[i] = calc_lpd(oph_mk[i], eph_mk[i], at, beta);
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
  vector[N_Mk] AE_Alpha_Mk;

  int N_Dedup_1_All;
  array[N_Mk] int N_Dedup_1_Mk;
  array[2 * Hnc, N_Dedup_1_All] int Index_Dedup_1;
  vector[N_Dedup_1_All] Size_Dedup_1;
  vector[N_Dedup_1_All] Oph_Dedup_1;

  int N_Dedup_2_All;
  array[N_Mk] int N_Dedup_2_Mk;
  array[N_Dedup_2_All] int Index_Dedup_2;

  array[N_Mk] int N_Peak_Mk;
  array[N_Mk] int N_Gtc_Mk;
  vector[sum(N_Gtc_Mk)] Gtc_Logprob_All;

  vector[N_Mk] At_Mk;
  real Beta;
}

transformed data{
  real tph_mean = (Ltph + Htph) / 2;
  real tph_sd = (Htph - Ltph) / 2;
}

parameters{
  real<lower = Ltph, upper = Htph> tph;
  simplex[Hnc] mp;
  vector<lower = -0.05, upper = 0>[Hnc] deg;
  simplex[N_Mk] ae;
}

model{
  tph ~ normal(tph_mean, tph_sd);
  mp ~ dirichlet(rep_vector(1, Hnc));
  deg ~ normal(-0.025, 0.0125);
  ae ~ dirichlet(AE_Alpha_Mk);

  vector[N_Dedup_1_All] eph_dedup_1 = calc_eph_dedup_1(tph, mp, deg, Hnc, N_Dedup_1_All, Eph_Base, Index_Dedup_1, Size_Dedup_1);

  int eph_start = 1;
  int dedup2_start = 1;
  int gtc_start = 1;
  for(i in 1:N_Mk){
    int n_dedup_eph = N_Dedup_1_Mk[i];
    vector[n_dedup_eph] eph_mk = N_Mk * ae[i] * eph_dedup_1[eph_start:(eph_start + n_dedup_eph - 1)];
    vector[n_dedup_eph] oph_mk = Oph_Dedup_1[eph_start:(eph_start + n_dedup_eph - 1)];
    real at = At_Mk[i];

    vector[n_dedup_eph] lpd_mk = calc_lpd_mk(n_dedup_eph, oph_mk, eph_mk, at, Beta);

    int n_dedup_2 = N_Dedup_2_Mk[i];
    int n_gtc = N_Gtc_Mk[i];
    int n_peak = N_Peak_Mk[i];
    array[n_dedup_2] int index_dedup_2_mk = Index_Dedup_2[dedup2_start:(dedup2_start + n_dedup_2 - 1)];
    matrix[n_peak, n_gtc] lpd_extract = to_matrix(lpd_mk[index_dedup_2_mk], n_peak, n_gtc);
    vector[n_gtc] gtc_logprob_mk = Gtc_Logprob_All[gtc_start:(gtc_start + n_gtc - 1)];
    vector[n_gtc] lpd_gtc = calc_lpd_gtc(n_gtc, lpd_extract, gtc_logprob_mk);

    target += log_sum_exp(lpd_gtc);

    eph_start += n_dedup_eph;
    dedup2_start += n_dedup_2;
    gtc_start += n_gtc;
  }
}
