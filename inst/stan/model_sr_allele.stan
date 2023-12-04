data{
  int N;
  vector[N] SR;
  vector[N] Allele;
  vector[N] Height;
}

parameters{
  real<lower = 0> slope;
  real intercept;
  real<lower = 0> sigma;
}

model{
  for(i in 1:N){
    SR[i] ~ lognormal(log(slope * Allele[i] + intercept), sqrt(sigma / Height[i]));
  }
}
