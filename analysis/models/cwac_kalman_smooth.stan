
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
  vector[N] dt;
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
  real<lower = 0> sig_alpha0;
  real<lower = 0> sig_beta0;
  real<lower = 0> sig_lambda0;

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
  vector[M] a[N + 1];
  // Predicted covariance
  matrix[M, M] P[N + 1];
  // Estimated states
  vector[M] ahat[N];
  // Estimated covariances
  matrix[M, M] Phat[N];
  // Transition matrices
  matrix[M,M] At[N];
  matrix[M,M] Qt[N];
  // Likelihood
  vector[N] llik_obs = rep_vector(0, N);

  // 1st observation
  vector[3] aa1 = [mu0, beta0, lambda0]';
  vector[3] pp1 = [sig_alpha0, sig_beta0, sig_lambda0]';
  a[1] = aa1;
  P[1] = diag_matrix(pp1);

  // precompute all transition matrices
  for(i in 1:N){
    // Transition matrices
    At[i] = [ [1, summer[i], winter[i]], [0, 1, 0], [0, 0, 1] ];
    // Covariance matrices
    Qt[i] = [ [sig_alpha, 0, 0], [0, sig_beta, 0], [0, 0, sig_lambda] ] * sqrt(dt[i]);
    // observation error
    R[i] = sig_s * summer[i] + sig_w * winter[i] + sig_o - sig_o * (summer[i] + winter[i]);
  }

  // iterative estimation
  for (i in 1:N) {
    matrix[M,M] Bt;
    vector[D] v;
    matrix[D, D] S;
    // vector[D] Sdiag;
    matrix[D, D] Sinv;
    matrix[M, D] K;
    matrix[M, M] L;

    // Prediction
    ahat[i] = At[i] * a[i];

    Phat[i] = At[i] * P[i] * At[i]' + Qt[i];

    // Update
    v = y[i] - H * ahat[i];
    S = H * Phat[i] * H' + R[i];
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
  mu0 ~ normal(0, 5);
  beta0 ~ normal(0, 5);
  lambda0 ~ normal(0, 5);
  sig_alpha0 ~ cauchy(0, 5);
  sig_beta0 ~ cauchy(0, 5);
  sig_lambda0 ~ cauchy(0, 5);

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
