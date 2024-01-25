functions{
  real calc_density(int n_peak, vector oph, vector eph, int at, real beta){
    real dens;
    real oph_one;
    real eph_one;
    dens = 0;
    for(i in 1:n_peak){
      oph_one = oph[i];
      eph_one = eph[i];
      if(eph_one != 0){ // other than 99
        if(oph_one >= at){
          dens += gamma_lpdf(oph_one | eph_one * beta, beta);
        }else{
          for(j in 1:(at - 1)){
            dens += gamma_lpdf(j | eph_one * beta, beta);
          }
        }
      }
    }
    return dens;
  }

  real calc_loglike(int n_peak, int n_gt_comb, vector eph_cspmk, vector oph, int at, real beta){
    vector[n_gt_comb] density;
    density = rep_vector(0, n_gt_comb);

    int index_gtc_start = 1;
    for(i in 1:n_gt_comb){
      vector[n_peak] eph = eph_cspmk[index_gtc_start:(index_gtc_start + n_peak - 1)];
      density[i] = calc_density(n_peak, oph, eph, at, beta);
      index_gtc_start += n_peak;
    }
    real loglike = log_sum_exp(density);
    return(loglike);
  }
}

data{
  int Hnc;
  int N_Csp_Mk;
  array[N_Csp_Mk] int N_Gt_Comb_All;
  array[N_Csp_Mk] int N_Peak_All;
  array[N_Csp_Mk] int N_Eph_Base_All;

  matrix[sum(N_Eph_Base_All), Hnc] Eph_Base_All;
  vector[sum(N_Eph_Base_All)] Size_Eph_All;

  vector[sum(N_Peak_All)] Oph_All;
  array[N_Csp_Mk] int AT_All;
  real Beta;
}

parameters{
  positive_ordered[Hnc] mp;
  vector<lower = -0.05, upper = 0>[Hnc] deg;
}

transformed parameters{
  vector[sum(N_Eph_Base_All)] eph_all;
  eph_all = rep_vector(0, sum(N_Eph_Base_All));
  for(i in 1:Hnc){
    real mp_c1 = mp[i];
    real deg_c1 = deg[i];
    eph_all += Eph_Base_All[, i] * mp_c1 .* exp(deg_c1 * Size_Eph_All);
  }
}

model{
  mp ~ uniform(0, 32000);
  deg ~ uniform(-0.05, 0);

  int index_eph_start = 1;
  int index_oph_start = 1;
  for(i in 1:N_Csp_Mk){
    int n_peak = N_Peak_All[i];
    int n_gt_comb = N_Gt_Comb_All[i];
    int n_eph_base = N_Eph_Base_All[i];
    vector[n_eph_base] eph_cspmk = eph_all[index_eph_start:(index_eph_start + n_eph_base - 1)];
    vector[n_peak] oph = Oph_All[index_oph_start:(index_oph_start + n_peak - 1)];
    int at = AT_All[i];
    target += calc_loglike(n_peak, n_gt_comb, eph_cspmk, oph, at, Beta);

    index_eph_start += n_eph_base;
    index_oph_start += n_peak;
  }
}
