data{
  int N;
  matrix[10, N] RPH;
  vector[N] TPH;
  vector[N] OPH;
  vector[N] AT;
}

parameters{
  real<lower = 0> cv_allele;
  real<lower = 0> cv_stutter;
}

model{
  for(i in 1:N){
    vector[10] tmp1 = RPH[, i] * TPH[i];
    vector[10] tmp2 = tmp1;
    for(j in 1:2){
      tmp2[j] *= pow(cv_allele, 2);
    }
    for(j in 3:10){
      tmp2[j] *= pow(cv_stutter, 2);
    }

    real tmp1_sum = sum(tmp1);
    real tmp2_sum = sum(tmp2);
    real gamma_alpha = tmp1_sum^2 / tmp2_sum;
    real gamma_beta = tmp1_sum / tmp2_sum;

    if(OPH[i] == 0){
      target += gamma_lcdf(AT[i] | gamma_alpha, gamma_beta);
    }else{
      target += gamma_lpdf(OPH[i] | gamma_alpha, gamma_beta);
    }
  }
}
