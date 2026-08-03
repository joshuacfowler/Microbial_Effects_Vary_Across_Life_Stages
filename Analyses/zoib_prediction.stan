functions {
  //  * Returns:
  //   matrix of scaled group-level effects
  //
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }
  
}
data {
// --- 1. TRAINING SIZES (Required to read your old CSV columns safely) ---
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
  
  // --- 2. NEW PREDICTION DATA FRAME ---
  int<lower=1> N_new;                       // Number of rows in prediction df
  matrix[N_new, K] X_new;                   // New design matrix
  array[N_new] int<lower=1> J_1_new;        // group IDs for new rows
  array[N_new] int<lower=1> J_2_new;        // group IDs for new rows
  array[N_new] int<lower=1> J_3_new;        // group IDs for new rows
  array[N_new] int<lower=1> J_4_new;        // group IDs for new rows
  
  vector[N_new] Z_1_1_new;
  
  vector[N_new] Z_2_1_new;
  vector[N_new] Z_2_2_new;
  vector[N_new] Z_2_3_new;
  vector[N_new] Z_2_4_new;
  
  vector[N_new] Z_3_1_new;
  vector[N_new] Z_3_2_new;
  vector[N_new] Z_3_3_new;
  vector[N_new] Z_3_4_new;
  
  vector[N_new] Z_4_1_new;
  
  // Programmatic re_formula switch (0 = fixed only, 1 = include RE)
  int<lower=0, upper=1> fitted;
  array[4] int<lower=0, upper=1> inc_re;
}
transformed data {
    // vector<lower=0>[N] se2 = square(se);
  
  //  set up sequence for multi-threading
  // array[N] int sequence = linspaced_int_array(N, 1, N);

  // implementing a zero-one-inflated beta distribution by transforming range of RII to 0,1
  vector[N] Y_new = (Y + 1.0) / 2.0;
  vector<lower=0>[N] se_new = se / 2.0;
  
  // counts of exact 0s and 1s
  int<lower=0> N_zero = 0;
  int<lower=0> N_one = 0;
  int<lower=0> N_cont = 0;
  for (n in 1:N) {
    if (Y_new[n] == 0.0) N_zero += 1;
    else if (Y_new[n] == 1.0) N_one += 1;
    else N_cont += 1;
  }
  
  // Store indices and values of continuous data
  array[N_zero] int id_zero;
  array[N_one] int id_one;
  array[N_cont] int id_cont;
  
  vector[N_cont] Y_cont;
  vector[N_cont] se_cont;

  int pos_z = 1;
  int pos_o = 1;
  int pos_c = 1;
  for (n in 1:N) {
      if (Y_new[n] == 0.0) {
        id_zero[pos_z] = n;
        pos_z += 1;
      } else if (Y_new[n] == 1.0) {
        id_one[pos_o] = n;
        pos_o += 1;
      } else {
        id_cont[pos_c] = n;
        Y_cont[pos_c] = Y_new[n];
        se_cont[pos_c] = se_new[n];
        pos_c += 1;
      }
    }
}

parameters {
  // incorporating measurement error using a latent true scores for the continuous measurements
  vector<lower=0.0, upper=1.0>[N_cont] Y_cont_true;

  vector[K] b;  // regression coefficients
  vector[K] b_phi;  // dispersion intercept coefficients
  real<lower=0> sd_phi_2;
  vector[N_2] z_phi_2; // standardized group-level effects (varying slopes and intercepts)

  

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
  
   
  real theta_zero; // fixed effects on zero-one inflation (intercepts)
  real theta_one;
    
  // allowing inflation parameters to vary with the random effects for host family
  // real<lower=0> sd_zero_2;
  // real<lower=0> sd_one_2;
  
  // vector[N_2] z_zero_2;// standardized group-level effects (varying family level intercepts only)
  // vector[N_2] z_one_2;

  // 3D simplex to model 0,1,continuous boundary
  // simplex[3] theta;
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
  
  // random effects on 01 inflation
  // vector[N_2] r_zero_2 = sd_zero_2 * z_zero_2;
  // vector[N_2] r_one_2  = sd_one_2  * z_one_2;
  
    // random effects on phi inflation
      vector[N_2] r_phi_2 = sd_phi_2 * z_phi_2;

  // matrix[N_2, M_2] r_phi_2 = transpose(diag_pre_multiply(sd_phi_2, z_phi_2));  // actual group-level effects
  // vector[N_2] r_phi_2_1 = r_phi_2[, 1];
  // vector[N_2] r_phi_2_2 = r_phi_2[, 2];
  // vector[N_2] r_phi_2_3 = r_phi_2[, 3];
  // vector[N_2] r_phi_2_4 = r_phi_2[, 4];

  // using vectors speeds up indexing in loops

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
  lprior += normal_lpdf(b | 0,1);
  lprior += normal_lpdf(b_phi | 3.8, 1);
    lprior += normal_lpdf(sd_phi_2 | 0, 1.0)
    - 1 * normal_lccdf(0 | 0,1);
  lprior += normal_lpdf(theta_zero | -4, 1);
  lprior += normal_lpdf(theta_one | -4, 1);
  // lprior += normal_lpdf(sd_zero_2 | 0, 1.0)
  //   - 1 * normal_lccdf(0 | 0,1);
  // lprior += normal_lpdf(sd_one_2  | 0, 1.0)
  //   - 1 * normal_lccdf(0 | 0,1);
  // lprior += dirichlet_lpdf(theta | rep_vector(1.0, 3));
  lprior += normal_lpdf(sd_1 | 0,1)
    - 1 * normal_lccdf(0 | 0,1);
  lprior += normal_lpdf(sd_2 | 0,1)
    - 4 * normal_lccdf(0 | 0,1);
  lprior += lkj_corr_cholesky_lpdf(L_2 | 2.0);
  lprior += normal_lpdf(sd_3 | 0,1)
    - 4 * normal_lccdf(0 | 0,1);
  lprior += lkj_corr_cholesky_lpdf(L_3 | 2.0);
  lprior += normal_lpdf(sd_4 | 0,1)
    - 1 * normal_lccdf(0 | 0,1);
}
model {
}

generated quantities {
  vector[N_new] y_fitted;
  if(fitted == 1){
    // note that N and X here are different values from above because we fit in our prediction data frame

    
    // simplex theta_i
    vector[3] theta = softmax([theta_zero, theta_one, 0.0]');
    for (n_i in 1:N_new) {
    // reconstruct linear predictor with indicator to turn off Random effects
    real mu_lin = X_new[n_i] * b 
                + (r_1_1[J_1_new[n_i]]) * inc_re[1]
                + (r_2_1[J_2_new[n_i]] * Z_2_1_new[n_i] + r_2_2[J_2_new[n_i]] * Z_2_2_new[n_i] 
                + r_2_3[J_2_new[n_i]]  * Z_2_3_new[n_i] + r_2_4[J_2_new[n_i]] * Z_2_4_new[n_i]) * inc_re[2]
                + (r_3_1[J_3_new[n_i]]  * Z_3_1_new[n_i] + r_3_2[J_3_new[n_i]] * Z_3_2_new[n_i] 
                + r_3_3[J_3_new[n_i]]  * Z_3_3_new[n_i] + r_3_4[J_3_new[n_i]] * Z_3_4_new[n_i]) * inc_re[3]
                + (r_4_1[J_4_new[n_i]]) * inc_re[4];
    real mu_fitted = inv_logit(mu_lin);
    
    // Calculate expected value 
    real mu_cont_scaled = (mu_fitted * 2.0) - 1.0;
    
    y_fitted[n_i] = theta[1] * (-1.0) 
                  + theta[2] * (1.0) 
                  + theta[3] * mu_cont_scaled;
    }
  }
    
}
