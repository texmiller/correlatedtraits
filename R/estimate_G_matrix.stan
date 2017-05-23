data {
  int<lower=0> N;              // number of observations
  int<lower=0> NS;             // number of sires
  int<lower=0> ND;             // number of dams
  matrix[3] zeros[NS];         // matrix of zeros for mean sire effect
  matrix[ND, NS] DS;           // dam-sire model matrix
  vector<lower=0>[N] dis;      // observed dispersal distances
  vector<lower=0>[N] off;      // observed number of offspring

}

parameters {
  cov_matrix[3] S;             // sire variance-covariance matrix
  cov_matrix[3] D;             // dam variance-covariance matrix
  real<lower=0> phi_dis;       // dispersal overdispersion parameter
  real<lower=0> phi_off;       // offspring overdispersion parameter
  real beta_mu;                // intercept of mean dispersal distance
  real beta_lambda;            // intercept of low-density per-capita growth rate
  real beta_alpha;             // intercept of Beverton-Holt parameter
}

transformed parameters {
  vector<lower=0>[3] d_S[NS];  // sire effect
  vector<lower=0>[3] d_DS[ND]; // dam (nested within sire) effect
  real<lower=0> mu_DS;         // estimated mean dispersal distance
  real<lower=0> lambda_DS;     // estimated low-density per-capita growth rate
  real<lower=0> alpha_DS;      // estimated Beverton-Holt parameter
  real<lower=0> lambda;        // expected density-dependent growth rate

  // dam effects are nested within sire effects
  d_S  ~ multi_normal(zeros, S)
  d_DS ~ multi_normal(DS * d_S, D)

  // expected trait values for individuals with dam D and sire S
  mu_DS      = exp(d_DS[, 1] + beta_mu)
  lambda_DS  = exp(d_DS[, 2] + beta_lambda)
  alpha_DS   = exp(d_DS[, 3] + beta_alpha)

  // expected per-capita growth rate
  lambda = (lambda_DS * N) / (1 + (N / alpha_DS))
}

model {
  // observed dispersal distance
  dis ~ NegBinomial2(mu, phi_dis)

  // observed number of offspring
  off ~ NegBinomial2(lambda, phi_off)
}

generated quantities{
  cov_matrix[3] G;

  // calculate G-matrix
  G = 2 * (S + D)
}
