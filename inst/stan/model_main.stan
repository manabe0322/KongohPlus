data{
  int Hnc;
  int AT;
  int N_Sn;
  int N_Mk;
  array[N_Mk] int N_Gt_All;
  array[N_Mk] int N_Peak_All;
  int N_Eph_Base_All;
  vector[N_Eph_Base_All] Eph_Base_All;
  vector[N_Eph_Base_All] Size_Eph_All;
  int N_Oph_All;
  vector[N_Oph_All] Oph_All;
  real Beta;
}

parameters{
  positive_ordered[Hnc] mp;
  vector<lower = -0.05, upper = 0>[Hnc] deg;
}

transformed parameters{
  matrix[N_Eph_Base_All, Hnc] eph_per_contributor;
  for(i in 1:Hnc){
    real mp_c1 = mp[i];
    real deg_c1 = deg[i];
    eph_per_contributor[, i] = Eph_Base_All * mp_c1 .* exp(deg_c1 * Size_Eph_All);
  }
}

model{

}
