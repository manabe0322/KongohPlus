data{
  int N;
  vector[N] OPH;
  vector[N] RPH;
  vector[N] TPH;
  vector[N] AT;
}

parameters{
  real<lower = 0, upper = 30> cv;
}

model{
  for(i in 1:N){
    real sigma = cv / sqrt(TPH[i]);
    real gamma_alpha = pow(sigma, -2) * RPH[i];
    real gamma_beta = 1 / (TPH[i] * pow(sigma, 2));
    if(OPH[i] == 0){
      target += gamma_lcdf(AT[i] | gamma_alpha, gamma_beta);
    }else{
      target += gamma_lpdf(OPH[i] | gamma_alpha, gamma_beta);
    }
  }
}
