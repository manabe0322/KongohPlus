data{
  int N;
  vector[N] OPH;
  vector[N] RPH;
  vector[N] TPH;
}

parameters{
  real<lower = 0, upper = 10> cv;
}

model{
  for(i in 1:N){
    real sigma = cv / sqrt(TPH[i]);
    OPH[i] ~ gamma(pow(sigma, -2) * RPH[i], 1 / (TPH[i] * pow(sigma, 2)));
  }
}
