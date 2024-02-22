data{
  int N;
  int N_Mk;
  matrix[N_Mk, N] AE;
}

parameters{
  vector<lower = 0>[N_Mk] alpha;
}

model{
  for(i in 1:N){
    AE[, i] ~ dirichlet(alpha);
  }
}
