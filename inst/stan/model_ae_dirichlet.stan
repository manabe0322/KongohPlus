data{
  int N;
  int N_Mk;
  vector[N] AE;
}

parameters{
  vector<lower = 0>[N_Mk] alpha;
}

model{
  for(i in 1:N){
    mp ~ dirichlet();

    AE[i] ~ lognormal(log(myu), sqrt(sigma / Height[i]));
  }
}
