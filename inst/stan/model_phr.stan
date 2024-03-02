data{
  int N;
  vector[N] PHR;
  vector[N] Mean_Ph;
}

parameters{
  real<lower = 0> sigma;
}

model{
  for(i in 1:N){
    PHR[i] ~ lognormal(0, sqrt(sigma / Mean_Ph[i]));
  }
}
