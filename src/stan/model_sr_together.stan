data{
  int N;
  int NMk;
  vector[N] SR;
  vector[N] Height;
  array[N] int IndexMk;
}

parameters{
  vector<lower = 0>[NMk] myu;
  real<lower = 0> sigma;
}

model{
  for(i in 1:N){
    SR[i] ~ lognormal(log(myu[IndexMk[i]]), sqrt(sigma / Height[i]));
  }
}
