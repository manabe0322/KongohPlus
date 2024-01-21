functions{
  real calc_likelihood(int Hnc, int N_Sn_Mk, int N_Peak_All, array[] int N_Eph_Base_All, array[] int N_Gt_Comb_All){
    index_eph_start = 1;
    index_gtc_start = 1;
    for(i in 1:N_Sn_Mk){
      int N_Eph_Base = N_Eph_Base_All[i];
      matrix[Hnc, N_Eph_Base] eph_per_person_snmk1 = eph_per_person[, index_eph_start:(index_eph_start + N_Eph_Base - 1)];

      int N_Gt_Comb = N_Gt_Comb_All[i];
      array[Hnc, N_Gt_Comb] int index_gt_comb_snmk1 = Index_Gt_Comb_All[, index_gtc_start:(index_gtc_start + N_Gt_Comb - 1)];

      int N_Peak = N_Peak_All[i];
      for(j in 1:N_Gt_Comb){
        array[Hnc] int index_gtc1 = index_gt_comb_snmk1[, j];

        vector[N_Peak] eph;
        for(k in 1:Hnc){
          int index_gt1 = index_gtc1[k];
          eph += eph_per_person_snmk1[k, (N_Peak * (index_gt1 - 1) + 1):(N_Peak * index_gt1)];
        }

        // vs oph
      }
    }
  }
}

data{
  int Hnc;

  int N_Sn_Mk;
  //array[N_Sn_Mk] int N_Gt_All;
  array[N_Sn_Mk] int N_Peak_All;
  array[N_Sn_Mk] int N_Eph_Base_All; //N_Gt_All * N_Peak_All

  vector[sum(N_Eph_Base_All)] Eph_Base_All;
  vector[sum(N_Eph_Base_All)] Size_Eph_All;

  array[N_Sn_Mk] int N_Gt_Comb_All;
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

}
