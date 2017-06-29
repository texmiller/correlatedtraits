// Large portions of this code were originally developed by Diogo Melo. His
// original code can be found here:
// https://github.com/diogro/stanAnimal/blob/master/animalModel.stan

functions {
  // cholesky kronecker product
  vector chol_kron_prod(matrix LA, vector sigma, matrix L_Omega, vector bv) {
    matrix[rows(L_Omega), cols(L_Omega)] LG;
    vector[num_elements(bv)] a;
    int ii;
    int jj;
    int col_LG;

    LG     = diag_pre_multiply(sigma, L_Omega);
    a      = rep_vector(0, num_elements(bv));
    col_LG = cols(LG);

    for(iA in 1:cols(LA)){
      for(jA in 1:iA){
        // don't calculate products between unrelated individuals
        if(LA[iA, jA] > 1e-10){
          for(iG in 1:col_LG){
            for(jG in 1:iG){
              ii    = (col_LG * (iA - 1)) + iG;
              jj    = (col_LG * (jA - 1)) + jG;
              a[ii] = a[ii] + LA[iA, jA] * LG[iG, jG] * bv[jj];
            }
          }
        }
      }
    }
    return a;
  }
}

data {
  int<lower=1> M;      // number of traits
  int<lower=0> N;      // number of individuals in pedigree
  int<lower=0> ND;     // number of dispersal observations
  int<lower=0> NH;     // number of offspring observations
  int<lower=0> D [ND]; // observed dispersal distances
  int<lower=0> H [NH]; // observed number of offspring produced
  int<lower=0> iD[ND]; // vector to identify observed dispersal values
  int<lower=0> iH[NH]; // vector to identify observed low-density growth values
  vector  [NH] den;    // population densities
  matrix [N,N] A;      // relationship matrix
}

transformed data{
  matrix[N, N] LA;
  LA = cholesky_decompose(A); // cholesky-decomposed relationship matrix
}

parameters {
  real  <lower=0>          b;         // Bev-Holt interaction parameter
  vector<lower=0>      [M] beta;      // fixed effects
  vector<lower=0>      [M] sigma_G;   // stddevs (for G matrix)
  vector<lower=0>      [M] sigma_R;   // standard deviations (for R matrix)
  cholesky_factor_corr [M] L_Omega_G; // Cholesky-decomposed corr matrix, G
  cholesky_factor_corr [M] L_Omega_R; // Cholesky-decomposed corr matrix, R
  matrix             [M,N] z;         // vector of standard normal deviates
  vector             [N*M] bv;        // breeding values
}

transformed parameters {
  matrix [N,M] a;     // sampled additive genetic deviances
  matrix [N,M] d;     // sampled non-additive genetic deviances
  matrix [N,M] apd;   // sum of a and d (a plus d)
  vector [ND]  mu_D;  // expected D values for observed individuals
  vector [NH]  mu_H;  // expected number of offspring for observed individuals
  vector [NH]  mu_r;  // expected r values for observed individuals

  // additive and non-additive genetic deviates
  a   = to_matrix(chol_kron_prod(LA, sigma_G, L_Omega_G, bv), N, M, 0);
  d   = (diag_pre_multiply(sigma_R, L_Omega_R) * z)';
  apd = a + d;

  // expected trait values
  mu_D = to_vector(apd[iD, 1]) + beta[1];
  mu_r = to_vector(apd[iH, 2]) + beta[2];
  mu_H = exp(mu_r) ./ (1 + b * den);
}

model {
  // observations
  D ~ poisson_log(mu_D);
  H ~ poisson(mu_H);

  // (unscaled) additive and non-additive genetic deviates
  bv           ~ normal(0, 1);
  to_vector(z) ~ normal(0, 1);
  // ^^^ this enables a non-centered parameterization ("the Matt trick""):
  // https://github.com/hmods/notes/blob/master/ch7/chapter7.Rmd

  // priors
  L_Omega_G    ~ lkj_corr_cholesky(2);
  L_Omega_R    ~ lkj_corr_cholesky(2);
  sigma_R      ~ normal(0, 2);
  sigma_G      ~ normal(0, 2);
  beta         ~ normal(0, 1);
}

generated quantities {
  cov_matrix [M]      G;
  cov_matrix [M]      E;
  cov_matrix [M]      P;
  vector     [M]     h2;
  corr_matrix[M]  corrG;
  corr_matrix[M]  corrE;
  vector     [ND] sim_D;
  vector     [NH] sim_H;

  G  = multiply_lower_tri_self_transpose(diag_pre_multiply(sigma_G, L_Omega_G));
  E  = multiply_lower_tri_self_transpose(diag_pre_multiply(sigma_R, L_Omega_R));
  P  = G + E;
  h2 = diagonal(G) ./ diagonal(P);

  corrG = multiply_lower_tri_self_transpose(L_Omega_G);
  corrE = multiply_lower_tri_self_transpose(L_Omega_R);

  // for some reason poisson_log_rng uses values that are too large sometimes,
  // so we need to constrain them somewhat.
  for(i in 1:ND) sim_D[i] = poisson_log_rng(mu_D[i] < 10  ? mu_D[i] :     10);
  for(i in 1:NH) sim_H[i] = poisson_rng(mu_H[i] < exp(10) ? mu_H[i] : exp(10));
}
