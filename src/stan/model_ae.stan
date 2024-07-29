data{
  int N;
  int N_Mk;
  matrix[N_Mk, N] PAE;
}

parameters{
  vector<lower = 0>[N_Mk] alpha;
}

model{
  for(i in 1:N){
    PAE[, i] ~ dirichlet(alpha);
  }
}
