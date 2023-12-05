data{
  int N;
  vector[N] SR;
  vector[N] Height;
}

parameters{
  real myu;
  real<lower = 0> sigma;
}

model{
  for(i in 1:N){
    SR[i] ~ lognormal(myu, sqrt(sigma / Height[i]));
  }
}
