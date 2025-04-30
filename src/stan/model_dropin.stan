data{
  int N;
  vector[N] Height;
}

parameters{
  real<lower = 0> shape_di;
  real<lower = 0> rate_di;
}

model{
  Height ~ gamma(shape_di, rate_di);
}
