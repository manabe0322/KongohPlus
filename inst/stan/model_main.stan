functions{
  real calc_density(int n_peak, vector oph, vector eph, int AT, real Beta){
    real dens;
    real oph_one;
    real eph_one;
    dens = 0;
    for(i in 1:n_peak){
      oph_one = oph[i];
      eph_one = eph[i];
      if(oph_one >= AT){
        dens += gamma_lpdf(oph_one | eph_one * Beta, Beta);
      }else{
        for(j in 1:(AT - 1)){
          dens += gamma_lpdf(j | eph_one * Beta, Beta);
        }
      }
    }
    return dens;
  }

  real calc_likelihood(int Hnc, int N_Csp_Mk, int N_Peak_All, array[] int N_Eph_Base_All, array[] int N_Gt_Comb_All, vector Oph_All, int AT, real Beta){
    index_eph_start = 1;
    index_gtc_start = 1;
    index_peak_start = 1;
    vector[N_Csp_Mk] = loglike_all;
    for(i in 1:N_Csp_Mk){
      int n_eph_base = N_Eph_Base_All[i];
      matrix[Hnc, n_eph_base] eph_per_person_cspmk1 = eph_per_person[, index_eph_start:(index_eph_start + n_eph_base - 1)];

      int n_gt_comb = N_Gt_Comb_All[i];
      array[Hnc, n_gt_comb] int index_gt_comb_cspmk1 = Index_Gt_Comb_All[, index_gtc_start:(index_gtc_start + n_gt_comb - 1)];

      int n_peak = N_Peak_All[i];
      vector[n_peak] oph = Oph_All[index_peak_start:(index_peak_start + n_peak - 1)];

      vector[n_gt_comb] density;
      for(j in 1:n_gt_comb){
        array[Hnc] int index_gtc1 = index_gt_comb_cspmk1[, j];

        vector[n_peak] eph;
        for(k in 1:Hnc){
          int index_gt1 = index_gtc1[k];
          eph += eph_per_person_cspmk1[k, (n_peak * (index_gt1 - 1) + 1):(n_peak * index_gt1)];
        }

        // vs oph
        density[j] = calc_density(n_peak, oph, eph, AT, Beta);
      }

      loglike_all[i] = log_sum_exp(density);

      index_eph_start += n_eph_base;
      index_gtc_start += n_gt_comb;
      index_peak_start += n_peak;
    }
    return(sum(loglike_all));
  }
}

data{
  int Hnc;

  int N_Csp_Mk;
  //array[N_Csp_Mk] int N_Gt_All;
  array[N_Csp_Mk] int N_Peak_All;
  array[N_Csp_Mk] int N_Eph_Base_All; //N_Gt_All * N_Peak_All

  vector[sum(N_Eph_Base_All)] Eph_Base_All;
  vector[sum(N_Eph_Base_All)] Size_Eph_All;

  array[N_Csp_Mk] int N_Gt_Comb_All;
  array[Hnc, sum(N_Gt_Comb_All)] int Index_Gt_Comb_All;

  int N_Oph_All;
  vector[N_Oph_All] Oph_All;
  int AT;
  real Beta;
}

parameters{
  positive_ordered[Hnc] mp;
  vector<lower = -0.05, upper = 0>[Hnc] deg;
}

transformed parameters{
  matrix[Hnc, N_Eph_Base_All] eph_per_person;
  for(i in 1:Hnc){
    real mp_c1 = mp[i];
    real deg_c1 = deg[i];
    eph_per_person[i, ] = Eph_Base_All * mp_c1 .* exp(deg_c1 * Size_Eph_All);
  }
}

model{
  target += calc_likelihood(Hnc, N_Csp_Mk, N_Peak_All, N_Eph_Base_All, N_Gt_Comb_All, Oph_All, AT, Beta);
}
