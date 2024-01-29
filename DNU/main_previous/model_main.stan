functions{
  real calc_density(int n_peak, vector oph, vector eph, int at, real Beta){
    real dens;
    real oph_one;
    real eph_one;
    dens = 0;
    for(i in 1:n_peak){
      oph_one = oph[i];
      eph_one = eph[i];
      if(eph_one != 0){ // other than 99
        if(oph_one >= at){
          dens += gamma_lpdf(oph_one | eph_one * Beta, Beta);
        }else{
          for(j in 1:(at - 1)){
            dens += gamma_lpdf(j | eph_one * Beta, Beta);
          }
        }
      }
    }
    return dens;
  }

  real calc_likelihood(int Hnc, int N_Csp_Mk, array[] int N_Peak_All, array[] int N_Eph_Base_All, matrix eph_per_person, array[] int N_Gt_Comb_All, array[,] int Index_Gt_Comb_All, vector Oph_All, array[] int AT_All, real Beta){
    int index_eph_start = 1;
    int index_gtc_start = 1;
    int index_peak_start = 1;
    vector[N_Csp_Mk] loglike_all;
    loglike_all = rep_vector(0, N_Csp_Mk);
    for(i in 1:N_Csp_Mk){
      int n_eph_base = N_Eph_Base_All[i];
      matrix[n_eph_base, Hnc] eph_per_person_cspmk1 = eph_per_person[index_eph_start:(index_eph_start + n_eph_base - 1), ];

      int n_gt_comb = N_Gt_Comb_All[i];
      array[n_gt_comb, Hnc] int index_gt_comb_cspmk1 = Index_Gt_Comb_All[index_gtc_start:(index_gtc_start + n_gt_comb - 1), ];

      int n_peak = N_Peak_All[i];
      vector[n_peak] oph = Oph_All[index_peak_start:(index_peak_start + n_peak - 1)];

      vector[n_gt_comb] density;
      density = rep_vector(0, n_gt_comb);

      int at = AT_All[i];

      for(j in 1:n_gt_comb){
        array[Hnc] int index_gtc1 = index_gt_comb_cspmk1[j, ];
        vector[n_peak] eph;
        eph = rep_vector(0, n_peak);
        for(k in 1:Hnc){
          vector[n_eph_base] eph_p1 = eph_per_person_cspmk1[, k];
          int index_gt1 = index_gtc1[k];
          eph += eph_p1[(n_peak * (index_gt1 - 1) + 1):(n_peak * index_gt1)];
        }

        // vs oph
        density[j] = calc_density(n_peak, oph, eph, at, Beta);
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
  array[sum(N_Gt_Comb_All), Hnc] int Index_Gt_Comb_All;

  vector[sum(N_Peak_All)] Oph_All;
  array[N_Csp_Mk] int AT_All;
  real Beta;
}

parameters{
  positive_ordered[Hnc] mp;
  vector<lower = -0.05, upper = 0>[Hnc] deg;
}

transformed parameters{
  matrix[sum(N_Eph_Base_All), Hnc] eph_per_person;
  for(i in 1:Hnc){
    real mp_c1 = mp[i];
    real deg_c1 = deg[i];
    eph_per_person[, i] = Eph_Base_All * mp_c1 .* exp(deg_c1 * Size_Eph_All);
  }
}

model{
  mp ~ uniform(0, 32000);
  deg ~ uniform(-0.05, 0);
  target += calc_likelihood(Hnc, N_Csp_Mk, N_Peak_All, N_Eph_Base_All, eph_per_person, N_Gt_Comb_All, Index_Gt_Comb_All, Oph_All, AT_All, Beta);
}
