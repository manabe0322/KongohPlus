data{
  int N;
  int NMk;
  vector[N] SR;
  vector[N] Height;
  array[N] int IndexMk;
}

parameters{
  vector[NMk] myu;
  real<lower = 0> sigma;
}

model{
  for(i in 1:N){
    SR[i] ~ lognormal(myu[IndexMk[i]], sqrt(sigma / Height[i]));
  }
}
