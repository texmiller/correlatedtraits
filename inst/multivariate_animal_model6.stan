functions {
  matrix as_matrix(vector X, int N, int M) {
    matrix[N, M] Y;

    for (i in 1:N) {
      Y[i] = to_row_vector(X[((i - 1) * M + 1):(i * M)]);
    }
    return Y;
  }

  vector chol_kronecker_product(matrix LA, matrix LG, vector a) {
    vector[num_elements(a)] new_a;

    new_a = rep_vector(0, num_elements(a));

    for(iA in 1:cols(LA)){
      for(jA in 1:iA){

        // avoid calculating products between unrelated individuals
        if(LA[iA, jA] > 1e-10){
          for(iG in 1:cols(LG)){
            for(jG in 1:iG){
              new_a[(cols(LG) * (iA - 1)) + iG] = new_a[(cols(LG) * (iA - 1)) + iG] + LA[iA, jA] * LG[iG, jG] * a[(cols(LG) * (jA - 1)) + jG];
            }
          }
        }
      }
    }
    return new_a;
  }

}

data {
  int<lower=1> M;      // number of traits
  int<lower=1> J;      // number of fixed effects
  int<lower=0> N;      // number of individuals in pedigree
  int<lower=0> ND;     // number of dispersal observations
  int<lower=0> Nr;     // number of low-density growth observations
  int<lower=0> NK;     // number of high-density growth observations
  int<lower=0> D[ND];  // observed dispersal distance
  int<lower=0> r[Nr];  // observed number of offspring produced
  int<lower=0> K[NK];  // observed number of (high density) offspring
  int<lower=0> iD[ND]; // vector to identify observed dispersal values
  int<lower=0> ir[Nr]; // vector to identify observed low-density growth values
  matrix[N, N] A;      // relationship matrix
}

transformed data{
  matrix[N, N] LA;
  LA = cholesky_decompose(A);
}

parameters {
  vector[N * M] a_tilde;          // breeding values
  vector[M]    beta;              // fixed effects
  vector[M]     Y[N];             // expected D and r values
  real         mu_K;              // carrying capacity mean (expected) value
  vector<lower=0, upper=10>[3] s; // negative binomial scale parameters

  # G matrix
  cholesky_factor_corr[M] L_Omega_G;
  vector<lower=0>     [M] L_sigma_G;

  // R matrix
  cholesky_factor_corr[M] L_Omega_R;
  vector<lower=0>     [M] L_sigma_R;
}

transformed parameters {
  matrix[M, M] L_Sigma_R;
  matrix[N, M] a;
  vector[ND]   mu_D;
  vector[Nr]   mu_r;
  vector[M]    mu_Dr[N];

  a = as_matrix(chol_kronecker_product(LA,
                                       diag_pre_multiply(L_sigma_G, L_Omega_G),
                                       a_tilde),
                                       N, M);
  L_Sigma_R = diag_pre_multiply(L_sigma_R, L_Omega_R);

  // // all mean trait values
  mu_D = to_vector(Y[iD, 1]);
  mu_r = to_vector(Y[ir, 2]);

  for(i in 1:N){
    mu_Dr[i, 1] = a[i, 1] + beta[1];
    mu_Dr[i, 2] = a[i, 2] + beta[2];
  }

}

model {

  // priors
  beta       ~ normal(0, 1);
  a_tilde    ~ normal(0, 1);
  L_sigma_R  ~ cauchy(0, 5);
  L_sigma_G  ~ cauchy(0, 5);
  L_Omega_G  ~ lkj_corr_cholesky(4);
  L_Omega_R  ~ lkj_corr_cholesky(4);

  // draw expected trait values
  Y ~ multi_normal_cholesky(mu_Dr, L_Sigma_R);

  // draw observations
  D ~ neg_binomial_2_log(mu_D, s[1]);
  r ~ neg_binomial_2_log(mu_r, s[2]);
  K ~ neg_binomial_2_log(mu_K, s[3]);
}

generated quantities {
  cov_matrix [M]     P;
  cov_matrix [M]     G;
  cov_matrix [M]     E;
  corr_matrix[M] corrG;
  corr_matrix[M] corrE;
  vector[2]          h;

  G = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_G, L_Omega_G));
  E = multiply_lower_tri_self_transpose(diag_pre_multiply(L_sigma_R, L_Omega_R));
  P = G + E;

  h[1] = G[1, 1] / P[1, 1];
  h[2] = G[2, 2] / P[2, 2];

  corrG = multiply_lower_tri_self_transpose(L_Omega_G);
  corrE = multiply_lower_tri_self_transpose(L_Omega_R);
}
