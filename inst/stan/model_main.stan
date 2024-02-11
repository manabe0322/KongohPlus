functions{
  vector calc_eph_dedup_1(real tph, vector mp, vector deg, int hnc, int n_dedup_1_all, vector eph_base, array[,] int index_dedup_1, vector size_dedup_1){
    vector[n_dedup_1_all] eph_dedup_1;
    eph_dedup_1 = rep_vector(0, n_dedup_1_all);
    for(i in 1:hnc){
      real tph_c1 = tph * mp[i];
      real deg_c1 = deg[i];
      eph_dedup_1 += eph_base[index_dedup_1[, 2 * i - 1]] * tph_c1 .* exp(deg_c1 * size_dedup_1);
      eph_dedup_1 += eph_base[index_dedup_1[, 2 * i]] * tph_c1 .* exp(deg_c1 * size_dedup_1);
    }
    return eph_dedup_1;
  }

  real calc_lpd_gen(real oph, real eph, real beta){
    real lpd = gamma_lpdf(oph | eph * beta, beta);
    return lpd;
  }

  real calc_lpd_do(real eph, int at_m1, real beta){
    vector[at_m1] lpd_drop;
    for(i in 1:at_m1){
      lpd_drop[i] = gamma_lpdf(i | eph * beta, beta);
    }
    real lpd = log_sum_exp(lpd_drop) - log(at_m1);
    return lpd;
  }

  vector calc_lpd_mk(int n_dedup_eph, int n_gen, int n_do, int n_di, int n_ig, vector oph_mk, vector eph_mk, int at_m1, real beta){
    vector[n_dedup_eph] lpd_mk;
    int index_current = 0;

    if(n_gen != 0){
      for(i in 1:n_gen){
        lpd_mk[i] = calc_lpd_gen(oph_mk[i], eph_mk[i], beta);
      }
      index_current += n_gen;
    }

    if(n_do != 0){
      for(i in 1:n_do){
        index_current += 1;
        lpd_mk[index_current] = calc_lpd_do(eph_mk[index_current], at_m1, beta);
      }
    }

    if(n_di != 0){
      for(i in 1:n_di){
        index_current += 1;
        lpd_mk[index_current] = 0;
      }
    }

    if(n_ig != 0){
      for(i in 1:n_ig){
        index_current += 1;
        lpd_mk[index_current] = 0;
      }
    }

    return lpd_mk;
  }

  vector calc_lpd_gtc(int n_peak, int n_gtc, array[] int index_dedup_2_mk, vector lpd_mk, vector gtc_logprob_mk){
    vector[n_gtc] lpd_gtc;
    for(i in 1:n_gtc){
      array[n_peak] int index_dedup_gtc = index_dedup_2_mk[(n_peak * (i - 1) + 1):(n_peak * i)];
      lpd_gtc[i] = sum(lpd_mk[index_dedup_gtc]) + gtc_logprob_mk[i];
    }
    return lpd_gtc;
  }
}

data{
  int Hnc;
  int N_Cspmk;

  int N_Eph_Base_All;
  vector[N_Eph_Base_All] Eph_Base;

  int N_Dedup_1_All;
  array[N_Cspmk] int N_Dedup_1_Cspmk;
  array[N_Dedup_1_All, 2 * Hnc] int Index_Dedup_1;
  vector[N_Dedup_1_All] Size_Dedup_1;
  vector[N_Dedup_1_All] Oph_Dedup_1;

  int N_Dedup_2_All;
  array[N_Cspmk] int N_Dedup_2_Cspmk;
  array[N_Dedup_2_All] int Index_Dedup_2;

  array[N_Cspmk] int N_Peak_Cspmk;
  array[N_Cspmk] int N_Gtc_Cspmk;
  vector[sum(N_Gtc_Cspmk)] Gtc_Logprob_All;

  array[N_Cspmk] int At_M1_Cspmk;
  real Beta;

  array[4 * N_Cspmk] int N_Peak_Type_All;
}

parameters{
  real<lower = 1, upper = 50000> tph;
  simplex[Hnc] mp;
  vector<lower = -0.05, upper = 0>[Hnc] deg;
}

model{
  mp ~ dirichlet(rep_vector(1, Hnc));

  vector[N_Dedup_1_All] eph_dedup_1 = calc_eph_dedup_1(tph, mp, deg, Hnc, N_Dedup_1_All, Eph_Base, Index_Dedup_1, Size_Dedup_1);

  int eph_start = 1;
  int dedup2_start = 1;
  int gtc_start = 1;
  for(i in 1:N_Cspmk){
    int n_dedup_eph = N_Dedup_1_Cspmk[i];
    vector[n_dedup_eph] eph_mk = eph_dedup_1[eph_start:(eph_start + n_dedup_eph - 1)];
    vector[n_dedup_eph] oph_mk = Oph_Dedup_1[eph_start:(eph_start + n_dedup_eph - 1)];
    int at_m1 = At_M1_Cspmk[i];

    int n_gen = N_Peak_Type_All[4 * (i - 1) + 1];
    int n_do = N_Peak_Type_All[4 * (i - 1) + 2];
    int n_di = N_Peak_Type_All[4 * (i - 1) + 3];
    int n_ig = N_Peak_Type_All[4 * i];

    vector[n_dedup_eph] lpd_mk = calc_lpd_mk(n_dedup_eph, n_gen, n_do, n_di, n_ig, oph_mk, eph_mk, at_m1, Beta);

    int n_dedup_2 = N_Dedup_2_Cspmk[i];
    int n_gtc = N_Gtc_Cspmk[i];
    int n_peak = N_Peak_Cspmk[i];
    array[n_dedup_2] int index_dedup_2_mk = Index_Dedup_2[dedup2_start:(dedup2_start + n_dedup_2 - 1)];
    vector[n_gtc] gtc_logprob_mk = Gtc_Logprob_All[gtc_start:(gtc_start + n_gtc - 1)];
    vector[n_gtc] lpd_gtc = calc_lpd_gtc(n_peak, n_gtc, index_dedup_2_mk, lpd_mk, gtc_logprob_mk);

    target += log_sum_exp(lpd_gtc);

    eph_start += n_dedup_eph;
    dedup2_start += n_dedup_2;
    gtc_start += n_gtc;
  }
}
