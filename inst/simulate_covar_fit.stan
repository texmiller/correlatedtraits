data {
  int<lower=0> N;              // number of observations
  int<lower=0> G;              // number of groups
  int<lower=0> family[G * N/4];              // number of observations
  vector[N/4 * G] lo_den_off;        // estimated female density
  vector[N/4 * G] hi_den_off;        // estimated female density
  row_vector[2] zeros[G];      // matrix of zeros for mean residual
}

parameters {
  #cov_matrix[2] S;             // variance-covariance matrix
  #row_vector[2] r[G];          // residuals for each group
  vector[G] a;          // residuals for each group
  vector[G] b;          // residuals for each group
  real beta_A;            // intercept of low-density per-capita growth rate
  real beta_B;             // intercept of Beverton-Holt parameter
  real<lower=0> sigmaa;
  real<lower=0> sigmab;
  real<lower=0> sda;
  real<lower=0> sdb;
}

transformed parameters {
  vector[N/4 * G] A;         // expected density-dependent growth rate
  vector[N/4 * G] B;          // expected density-dependent growth rate
  #row_vector[2] mus[G];
  vector[G] mua;
  vector[G] mub;

  for (i in 1:(N/4 * G)) {
    // expected trait values for individuals with dam D and sire S
    A[i] = a[family[i]];
    B[i] = b[family[i]];
  }

  for (i in 1:G) {
    #mus[i] = [beta_A, beta_B];
    mua[i]  = beta_A;
    mub[i]  = beta_B;
  }

}

model {

  // residuals drawn from variance-covariance matrix S
  #r ~ multi_normal(mus, S);

  a ~ normal(mua, sigmaa);
  b ~ normal(mub, sigmab);


  // observed number of offspring
  lo_den_off ~ normal(A, sda);
  hi_den_off ~ normal(B, sdb);
}

generated quantities{
  real est_A;            // intercept of low-density per-capita growth rate
  real est_B;             // intercept of Beverton-Holt parameter

  // calculate other variables of interest
  est_A = beta_A;
  est_B = beta_B;
}
