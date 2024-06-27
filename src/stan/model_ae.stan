data{
  int N;
  int N_Mk;
  matrix[N_Mk, N] PAE;
  vector[N] Tph;
}

parameters{
  vector<lower = 0>[N_Mk] alpha;
}

model{
  for(i in 1:N){
    PAE[, i] ~ dirichlet(Tph[i] * alpha);
  }
}
