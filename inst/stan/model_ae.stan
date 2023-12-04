data{
  int N;
  vector[N] AE;
  vector[N] Height;
}

parameters{
  real myu;
  real<lower = 0> sigma;
}

model{
  for(i in 1:N){
    AE[i] ~ lognormal(myu, sqrt(sigma / Height[i]));
  }
}
