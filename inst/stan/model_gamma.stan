data{
  int N;
  vector[N] OPH;
  vector[N] EPH;
}

parameters{
  real<lower = 0> beta;
}

model{
  for(i in 1:N){
    OPH[i] ~ gamma(beta * EPH[i], beta);
  }
}
