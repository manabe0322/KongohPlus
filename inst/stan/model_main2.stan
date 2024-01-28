functions{
  real calc_density(real oph, real eph, int at, real beta){
    real dens = 0;
    if(eph != 0){ // other than 99
      if(oph >= at){
        dens += gamma_lpdf(oph | eph * beta, beta);
      }else{
        for(i in 1:(at - 1)){
          dens += gamma_lpdf(i | eph * beta, beta);
        }
      }
    }
    return dens;
  }
}

data{
  int Hnc;
  int N_Csp_Mk;

  array[N_Csp_Mk] int N_Eph_Base_All;
  vector[sum(N_Eph_Base_All)] Eph_Base_All;

  array[N_Csp_Mk] int N_Dedup_Eph_All;
  array[sum(N_Dedup_Eph_All), 2 * Hnc] int Index_Dedup_Eph_All;
  vector[sum(N_Dedup_Eph_All)] Size_Eph_Dedup_All;
  vector[sum(N_Dedup_Eph_All)] Oph_Dedup_All;

  array[N_Csp_Mk] int N_Index_Dedup_Gt_All;
  array[sum(N_Index_Dedup_Gt_All)] int Index_Dedup_Gt_All;

  array[N_Csp_Mk] int N_Peak_All;

  array[N_Csp_Mk] int AT_All;
  real Beta;
}

parameters{
  positive_ordered[Hnc] mp;
  vector<lower = -0.05, upper = 0>[Hnc] deg;
}

transformed parameters{
  vector[sum(N_Dedup_Eph_All)] eph_dedup;
  eph_dedup = rep_vector(0, sum(N_Dedup_Eph_All));
  for(i in 1:Hnc){
    real mp_c1 = mp[i];
    real deg_c1 = deg[i];
    eph_dedup += Eph_Base_All[Index_Dedup_Eph_All[, 2 * i - 1]] * mp_c1 .* exp(deg_c1 * Size_Eph_Dedup_All);
    eph_dedup += Eph_Base_All[Index_Dedup_Eph_All[, 2 * i]] * mp_c1 .* exp(deg_c1 * Size_Eph_Dedup_All);
  }
}

model{
  mp ~ uniform(0, 32000);
  deg ~ uniform(-0.05, 0);

  int eph_start = 1;
  int peak_start = 1;
  for(i in 1:N_Csp_Mk){
    int n_dedup_eph = N_Dedup_Eph_All[i];
    vector[n_dedup_eph] eph_mk = eph_dedup[eph_start:(eph_start + n_dedup_eph - 1)];
    vector[n_dedup_eph] oph_mk = Oph_Dedup_All[eph_start:(eph_start + n_dedup_eph - 1)];
    int at = AT_All[i];

    vector[n_dedup_eph] dens_mk;
    for(j in 1:n_dedup_eph){
      dens_mk[j] = calc_density(oph_mk[j], eph_mk[j], at, Beta);
    }

    int n_gtc = N_Index_Dedup_Gt_All[i];
    int n_peak = N_Peak_All[i];
    vector[n_gtc] dens_gtc;
    for(j in 1:n_gtc){
      array[n_peak] int index_dedup_gt = Index_Dedup_Gt_All[peak_start:(peak_start + n_peak - 1)];
      dens_gtc[j] = sum(dens_mk[index_dedup_gt]); // + log(genotype prob)
      peak_start += n_peak;
    }
    target += log_sum_exp(dens_gtc);

    eph_start += n_dedup_eph;
  }
}
