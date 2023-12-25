data{
  int N;
  vector[N] AE;
  vector[N] Height;
}

parameters{
  real<lower = 0> myu;
  real<lower = 0> sigma;
}

model{
  for(i in 1:N){
    AE[i] ~ lognormal(log(myu), sqrt(sigma / Height[i]));
  }
}
