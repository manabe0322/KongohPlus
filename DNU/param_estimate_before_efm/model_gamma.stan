data{
  int N;
  vector[N] OPH;
  vector[N] EPH;
  vector[N] AT;
}

parameters{
  real<lower = 0, upper = 0.1> beta;
}

model{
  for(i in 1:N){
    if(OPH[i] < AT[i]){
      target += gamma_lcdf(AT[i] | EPH[i] * beta, beta);
    }else{
      OPH[i] ~ gamma(EPH[i] * beta, beta);
    }
  }
}
