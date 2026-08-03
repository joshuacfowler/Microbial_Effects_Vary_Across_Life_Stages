functions {
 /* compute correlated group-level effects
  * Args:
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns:
  *   matrix of scaled group-level effects
  */
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }
  
  
  // parallelized likelihood truncation calculation function
  real partial_log_lik(array[] int seq, int start, int end,
                       vector Y, vector mu, vector total_sigma,
                       real lb, real ub) {
    int N_slice = end - start + 1;
    vector[N_slice] Y_slice     = Y[start:end];
    vector[N_slice] mu_slice    = mu[start:end];
    vector[N_slice] sigma_slice = total_sigma[start:end];
    // 1. Core Likelihood (Vectorized for this slice)
    real lp = normal_lpdf(Y_slice | mu_slice, sigma_slice);
    // 2. Truncation Correction (MUST loop element-wise)
    for (i in 1:N_slice) {
      lp -= log_diff_exp(normal_lcdf(ub | mu_slice[i], sigma_slice[i]), 
                         normal_lcdf(lb | mu_slice[i], sigma_slice[i]));
    }
    return lp;
  }

}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  vector<lower=0>[N] se;  // known sampling error
  real lb;  // lower truncation bounds;
  real ub;  // upper truncation bounds
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  array[N] int<lower=1> J_2;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_1;
  vector[N] Z_2_2;
  vector[N] Z_2_3;
  vector[N] Z_2_4;
  int<lower=1> NC_2;  // number of group-level correlations
  // data for group-level effects of ID 3
  int<lower=1> N_3;  // number of grouping levels
  int<lower=1> M_3;  // number of coefficients per level
  array[N] int<lower=1> J_3;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_3_1;
  vector[N] Z_3_2;
  vector[N] Z_3_3;
  vector[N] Z_3_4;
  int<lower=1> NC_3;  // number of group-level correlations
  // data for group-level effects of ID 4
  int<lower=1> N_4;  // number of grouping levels
  int<lower=1> M_4;  // number of coefficients per level
  array[N] int<lower=1> J_4;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_4_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  vector<lower=0>[N] se2 = square(se);
  array[N] int sequence = linspaced_int_array(N, 1, N);
}
parameters {
  vector[K] b;  // regression coefficients
  real<lower=0> sigma;  // dispersion parameter
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  array[M_1] vector[N_1] z_1;  // standardized group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  matrix[M_2, N_2] z_2;  // standardized group-level effects
  cholesky_factor_corr[M_2] L_2;  // cholesky factor of correlation matrix
  vector<lower=0>[M_3] sd_3;  // group-level standard deviations
  matrix[M_3, N_3] z_3;  // standardized group-level effects
  cholesky_factor_corr[M_3] L_3;  // cholesky factor of correlation matrix
  vector<lower=0>[M_4] sd_4;  // group-level standard deviations
  array[M_4] vector[N_4] z_4;  // standardized group-level effects
}
transformed parameters {
  vector[N_1] r_1_1;  // actual group-level effects
  matrix[N_2, M_2] r_2;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_2] r_2_1;
  vector[N_2] r_2_2;
  vector[N_2] r_2_3;
  vector[N_2] r_2_4;
  matrix[N_3, M_3] r_3;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_3] r_3_1;
  vector[N_3] r_3_2;
  vector[N_3] r_3_3;
  vector[N_3] r_3_4;
  vector[N_4] r_4_1;  // actual group-level effects
  // prior contributions to the log posterior
  real lprior = 0;
  r_1_1 = (sd_1[1] * (z_1[1]));
  // compute actual group-level effects
  r_2 = scale_r_cor(z_2, sd_2, L_2);
  r_2_1 = r_2[, 1];
  r_2_2 = r_2[, 2];
  r_2_3 = r_2[, 3];
  r_2_4 = r_2[, 4];
  // compute actual group-level effects
  r_3 = scale_r_cor(z_3, sd_3, L_3);
  r_3_1 = r_3[, 1];
  r_3_2 = r_3[, 2];
  r_3_3 = r_3[, 3];
  r_3_4 = r_3[, 4];
  r_4_1 = (sd_4[1] * (z_4[1]));
  lprior += normal_lpdf(b | 0,.2);
  lprior += normal_lpdf(sigma | 0,.5)
    - 1 * normal_lccdf(0 | 0,.5);
  lprior += normal_lpdf(sd_1 | 0,.5)
    - 1 * normal_lccdf(0 | 0,.5);
  lprior += normal_lpdf(sd_2 | 0,.5)
    - 4 * normal_lccdf(0 | 0,.5);
  lprior += lkj_corr_cholesky_lpdf(L_2 | 2.0);
  lprior += normal_lpdf(sd_3 | 0,.5)
    - 4 * normal_lccdf(0 | 0,.5);
  lprior += lkj_corr_cholesky_lpdf(L_3 | 2.0);
  lprior += normal_lpdf(sd_4 | 0,.5)
    - 1 * normal_lccdf(0 | 0,.5);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    // trying to add the population levels effects in a vectorized way
    mu += X * b 
            + r_1_1[J_1] .* Z_1_1 
            + r_2_1[J_2] .* Z_2_1 
            + r_2_2[J_2] .* Z_2_2 
            + r_2_3[J_2] .* Z_2_3 
            + r_2_4[J_2] .* Z_2_4 
            + r_3_1[J_3] .* Z_3_1 
            + r_3_2[J_3] .* Z_3_2 
            + r_3_3[J_3] .* Z_3_3 
            + r_3_4[J_3] .* Z_3_4 
            + r_4_1[J_4] .* Z_4_1;
    // for (n in 1:N) {
    //   // add more terms to the linear predictor
    //   mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n] + r_2_2[J_2[n]] * Z_2_2[n] + r_2_3[J_2[n]] * Z_2_3[n] + r_2_4[J_2[n]] * Z_2_4[n] + r_3_1[J_3[n]] * Z_3_1[n] + r_3_2[J_3[n]] * Z_3_2[n] + r_3_3[J_3[n]] * Z_3_3[n] + r_3_4[J_3[n]] * Z_3_4[n] + r_4_1[J_4[n]] * Z_4_1[n];
    // }
    // variance and sampling error term as a vector
    vector[N] total_sigma = sqrt(square(sigma) + se2);

    // Likelihood model with reduce_sum for multithreading
    int grainsize = 1; 
    target += reduce_sum(partial_log_lik, sequence, grainsize,
                         Y, mu, total_sigma, lb, ub);
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(z_1[1]);
  target += std_normal_lpdf(to_vector(z_2));
  target += std_normal_lpdf(to_vector(z_3));
  target += std_normal_lpdf(z_4[1]);
}
generated quantities {
  // compute group-level correlations
  corr_matrix[M_2] Cor_2 = multiply_lower_tri_self_transpose(L_2);
  vector<lower=-1,upper=1>[NC_2] cor_2;
  // compute group-level correlations
  corr_matrix[M_3] Cor_3 = multiply_lower_tri_self_transpose(L_3);
  vector<lower=-1,upper=1>[NC_3] cor_3;
  // extract upper diagonal of correlation matrix
  for (k in 1:M_2) {
    for (j in 1:(k - 1)) {
      cor_2[choose(k - 1, 2) + j] = Cor_2[j, k];
    }
  }
  // extract upper diagonal of correlation matrix
  for (k in 1:M_3) {
    for (j in 1:(k - 1)) {
      cor_3[choose(k - 1, 2) + j] = Cor_3[j, k];
    }
  }
}
