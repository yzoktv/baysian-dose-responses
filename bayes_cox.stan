data {
  int<lower=1> N;
  int<lower=1> H;
  matrix[N, H] B;
  array[N] int<lower=0, upper=1> status;
  vector[N] time;

  int<lower=1> J;
  vector[J+1] time_knots;
  vector[H] dose_knots;
}

parameters {
  vector[H] beta_raw;                 // Non-centered version of beta
  vector<lower=0>[J] lambda;
  real<lower=0, upper=1> rho;
  real<lower=0> tau;
}

transformed parameters {
  vector[H] beta;
  beta[1] = beta_raw[1] * tau;
  for (h in 2:H) {
    real delta = dose_knots[h] - dose_knots[h - 1];
    beta[h] = pow(rho, delta) * beta[h - 1] + beta_raw[h] * tau;
  }
}

model {
  // Priors
  beta_raw ~ normal(0, 1);             // Standard normal for non-centered param
  rho ~ beta(4, 4);
  tau ~ normal(0, 1) T[0, ];
  lambda ~ gamma(2, 2);

  // Likelihood
  for (i in 1:N) {
    int t_idx = J;
    while (t_idx > 1 && time[i] <= time_knots[t_idx]) {
      t_idx -= 1;
    }

    real linpred = dot_product(B[i], beta);
    real log_haz = log(lambda[t_idx]) + linpred;
    real cumhaz = 0;

    for (j in 1:t_idx) {
      real t_start = time_knots[j];
      real t_end = (j == t_idx) ? time[i] : time_knots[j + 1];
      cumhaz += lambda[j] * (t_end - t_start) * exp(linpred);
    }

    target += status[i] * log_haz - cumhaz;
  }
}
