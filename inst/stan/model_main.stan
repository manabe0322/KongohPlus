functions{
  vector calc_eph_dedup_1(real tph, vector mp, vector deg, int Hnc, int N_Dedup_1, vector Eph_Base, array[,] int Index_Dedup_1, vector Size_Dedup_1){
    vector[N_Dedup_1] eph_dedup_1;
    eph_dedup_1 = rep_vector(0, N_Dedup_1);
    for(i in 1:Hnc){
      real tph_c1 = tph * mp[i];
      real deg_c1 = deg[i];
      eph_dedup_1 += Eph_Base[Index_Dedup_1[, 2 * i - 1]] * tph_c1 .* exp(deg_c1 * Size_Dedup_1);
      eph_dedup_1 += Eph_Base[Index_Dedup_1[, 2 * i]] * tph_c1 .* exp(deg_c1 * Size_Dedup_1);
    }
    return eph_dedup_1;
  }

  real calc_density(real oph, real eph, int at, real beta){
    real dens = 0;
    if(eph != 0){
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
  int N_Cspmk;

  int N_Eph_Base;
  vector[N_Eph_Base] Eph_Base;

  int N_Dedup_1;
  array[N_Cspmk] int N_Dedup_1_Cspmk;
  array[N_Dedup_1, 2 * Hnc] int Index_Dedup_1;
  vector[N_Dedup_1] Size_Dedup_1;
  vector[N_Dedup_1] Oph_Dedup_1;

  int N_Dedup_2;
  array[N_Dedup_2] int Index_Dedup_2;

  array[N_Cspmk] int N_Peak_Cspmk;
  array[N_Cspmk] int N_Gtc_Cspmk ;

  array[N_Cspmk] int At_Cspmk;
  real Beta;
}

parameters{
  real<lower = 1, upper = 50000> tph;
  simplex[Hnc] mp;

  //positive_ordered[Hnc] mp;
  vector<lower = -0.05, upper = 0>[Hnc] deg;
}

model{
  mp ~ dirichlet(rep_vector(1, Hnc));

  vector[N_Dedup_1] eph_dedup_1 = calc_eph_dedup_1(tph, mp, deg, Hnc, N_Dedup_1, Eph_Base, Index_Dedup_1, Size_Dedup_1);

  int eph_start = 1;
  int peak_start = 1;
  for(i in 1:N_Cspmk){
    int n_dedup_eph = N_Dedup_1_Cspmk[i];
    vector[n_dedup_eph] eph_mk = eph_dedup_1[eph_start:(eph_start + n_dedup_eph - 1)];
    vector[n_dedup_eph] oph_mk = Oph_Dedup_1[eph_start:(eph_start + n_dedup_eph - 1)];
    int at = At_Cspmk[i];

    vector[n_dedup_eph] dens_mk;
    for(j in 1:n_dedup_eph){
      dens_mk[j] = calc_density(oph_mk[j], eph_mk[j], at, Beta);
    }

    int n_gtc = N_Gtc_Cspmk[i];
    int n_peak = N_Peak_Cspmk[i];
    vector[n_gtc] dens_gtc;
    for(j in 1:n_gtc){
      array[n_peak] int index_dedup_gt = Index_Dedup_2[peak_start:(peak_start + n_peak - 1)];
      dens_gtc[j] = sum(dens_mk[index_dedup_gt]); // + log(genotype prob)
      peak_start += n_peak;
    }
    target += log_sum_exp(dens_gtc);

    eph_start += n_dedup_eph;
  }
}
