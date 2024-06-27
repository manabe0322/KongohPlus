data{
  int N;
  vector[N] OPH;
  vector[N] RPH;
  vector[N] TPH;
}

parameters{
  real<lower = 0, upper = 1> cv;
}

model{
  for(i in 1:N){
    OPH[i] ~ gamma(pow(cv, -2) * RPH[i], 1 / (TPH[i] * pow(cv, 2)));
  }
}
