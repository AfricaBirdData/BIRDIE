
data {
  // number of observations
  int N;
  // number of observed states
  int D;
  // number of states
  int M;
  // observed data
  vector[D] y[N];
  // season data
  vector[N] summer;
  vector[N] winter;
  // observation matrix
  matrix[D, M] H;
  // time intervals
  vector[N-1] dt;
  // Priors for initial values
  real pmu_mu0;
  real psig_mu0;
  real psig_beta0;
  real psig_lambda0;
  // Priors for rate parameters
  real<lower = 0> psig_alpha;
  real<lower = 0> psig_beta;
  real<lower = 0> psig_lambda;
  // Priors for observation error
  real<lower = 0> psig_s;
  real<lower = 0> psig_w;
  real<lower = 0> psig_o;
}

parameters {
  // Initial states
  real mu0;
  real beta0;
  real lambda0;

  // Standard deviation of rate parameters
  real<lower = 0> sig_alpha;
  real<lower = 0> sig_beta;
  real<lower = 0> sig_lambda;

  // Standard deviation of observation error
  real<lower = 0> sig_s;
  real<lower = 0> sig_w;
  real<lower = 0> sig_o;
}

transformed parameters {
  // observation errors
  vector[N] R;
  // Predicted states
  vector[M] a[N];
  // Predicted covariance
  matrix[M, M] P[N];
  // Estimated states
  vector[M] ahat[N];
  // Estimated covariances
  matrix[M, M] Phat[N];
  // Transition matrices
  matrix[M,M] At[N-1];
  matrix[M,M] Qt[N-1];
  // Likelihood
  vector[N] llik_obs = rep_vector(0, N);

  // variances
  real<lower = 0> var_alpha = pow(sig_alpha,2);
  real<lower = 0> var_beta = pow(sig_beta,2);
  real<lower = 0> var_lambda = pow(sig_lambda,2);
  real<lower = 0> var_s = pow(sig_s,2);
  real<lower = 0> var_w = pow(sig_w,2);
  real<lower = 0> var_o = pow(sig_o,2);

  // 1st observation
  vector[3] aa1 = [mu0, beta0, lambda0]';
  vector[3] pp1 = [1, 1, 1]';
  a[1] = aa1;
  P[1] = diag_matrix(pp1);
  R[1] = var_s * summer[1] + var_w * winter[1] + var_o - var_o * (summer[1] + winter[1]);
  llik_obs[1] += normal_lpdf(y[1]| mu0+beta0, sqrt(R[1]));

  // precompute all transition matrices
  for(i in 1:(N-1)){
    // Transition matrices
    At[i] = [ [1, summer[i], winter[i]], [0, 1, 0], [0, 0, 1] ];
    // Covariance matrices
    Qt[i] = [ [var_alpha, 0, 0], [0, var_beta, 0], [0, 0, var_lambda] ] * dt[i];
    // observation error
    R[i+1] = var_s * summer[i+1] + var_w * winter[i+1] + var_o - var_o * (summer[i+1] + winter[i+1]);
  }

  // iterative estimation
  for (i in 1:(N-1)) {
    vector[D] v;
    matrix[D, D] S;
    matrix[D, D] Sinv;
    matrix[M, D] K;
    matrix[M, M] L;

    // Prediction
    ahat[i] = At[i] * a[i];

    Phat[i] = At[i] * P[i] * At[i]' + Qt[i];

    // Update
    v = y[i+1] - H * ahat[i];
    S = H * Phat[i] * H' + R[i+1];
    K = mdivide_right_spd(Phat[i] * H', S);

    a[i + 1] = ahat[i] + K * v;
    P[i + 1] = Phat[i] - K * S * K';

    // manual update of multivariate normal
    for(d in 1:D){
      llik_obs[i] += normal_lpdf(v[d]| 0, sqrt(S[d,d]));
    }

  }

}

model {
  // Priors for initial states
  mu0 ~ normal(pmu_mu0, psig_mu0);
  beta0 ~ normal(0, psig_beta0);
  lambda0 ~ normal(0, psig_lambda0);

  // Priors for standard deviation of rate parameters
  sig_alpha ~ cauchy(0, psig_alpha);
  sig_beta ~ cauchy(0, psig_beta);
  sig_lambda ~ cauchy(0, psig_lambda);

  // Priors for standard deviation of observation error
  sig_s ~ cauchy(0, psig_s);
  sig_w ~ cauchy(0, psig_w);
  sig_o ~ cauchy(0, psig_o);

  // Likelihood
  target += sum(llik_obs);
}

generated quantities{

  // predicted states
  vector[M] pred[N];

  for(i in 1:N){
    pred[i] = multi_normal_rng(a[i], P[i]);
  }

}

// generated quantities{
//   // Smoothed states
//   vector[M] a_s[N];
//   matrix[M,M] P_s[N];
//   // inverse estimated covariances
//   matrix[M, M] invPhat[N];
//   matrix[M,M] G[N];
//
//   // Initial prediction
//   a_s[N] = a[N];
//   P_s[N] = P[N];
//
//   for(k in 2:N){
//     int i = N - k + 1;
//
//     invPhat[i+1] = inverse(Phat[i+1]);
//
//     // Kalman gain
//     G[i] = mdivide_right_spd(P[i] * At[i]', Phat[i+1]);
//
//     a_s[i] = a[i] + G[i] * (a_s[i+1] - ahat[i+1]);
//     P_s[i] = P[i] + G[i] * (P_s[i+1] - Phat[i+1]) * G[i]';
//   }
//
// }
