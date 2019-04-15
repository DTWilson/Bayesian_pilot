data {
  int<lower=0> K;  // # of clusters
    real<lower=0> M[K]; // # array of cluster sizes
      int<lower=0> N;  // total sample size
      int<lower=0> clus[N]; // cluster index allocations
      vector[N] y;  // responses
      int<lower=0,upper=1> a_p[N]; // adherences - patient level
      int<lower=0,upper=1> a_c[K/2]; // adherences - cluster level
      int<lower=0,upper=1> trt[N]; // treatment allocation
      int<lower=0> lost;  // number missed in follow-up
}
parameters {
  real<lower=0,upper=1> p_f; // prob of being followed-up
  real<lower=0,upper=1> p_a;  // prob of adhering to treatment
  real<lower=0> c_m; // Mean cluster size
  real<lower=0> c_sigsq; //cluster size variance
  real d; // Mean treatment effect
  real<lower=0,upper=1> d_rho;  // ICC for treatment effect
  real<lower=0> d_sigsq_w;  // whithin cluster variance for treatment effect
  
  vector[K] u;  // random effects
}
transformed parameters {
  real<lower=0> c_sig;    
  real<lower=0> d_sig_w;
  real<lower=0> d_sig_b;
  
  c_sig = sqrt(c_sigsq);
  d_sig_w = sqrt(d_sigsq_w);
  d_sig_b = sqrt(d_rho*d_sigsq_w/(1-d_rho));
}
model {
  // Priors
  //p_f ~ beta(22.4, 9.6);
  p_f ~ beta(1, 1);
  //p_a ~ beta(28.8, 3.2);
  p_a ~ beta(1, 1);
  c_sigsq ~ inv_gamma(20, 39);
  //c_sigsq ~ inv_gamma(1, 2);
  //c_m ~ normal(10, sqrt(c_sigsq/6));
  c_m ~ normal(10, 10);
  
  d_rho ~ beta(1.6, 30.4);
  //d_rho ~ beta(1, 1);
  d_sigsq_w ~ inv_gamma(50, 45);
  //d_sigsq_w ~ inv_gamma(1, 2);
  //d ~ normal(0.2, 0.25);
  d ~ normal(0.2, 1);
  
  u ~ normal(0, d_sig_b);
  
  // Follow-up
  N ~ binomial(N+lost, p_f);
  
  // Cluster size
  M ~ normal(c_m, c_sig);
  
  // Adherance
  a_c ~ bernoulli(p_a);
  
  for(i in 1:N){
    // Efficacy
    y[i] ~ normal(trt[i]*a_p[i]*d + u[clus[i]], d_sig_w);
  }
}

generated quantities {
  real<lower=0,upper=1> pr;
  real<lower=0,upper=1> pa;
  real<lower=0,upper=1> pg;
  
  pg = !(p_f < 0.6666666 || (22-15*p_f > c_m) || p_a < 0.6 || (1.05714286-0.5714286*d > p_a));
  pr = (p_f < 0.6 || (20-15*p_f > c_m) || p_a < 0.5 || (0.9571429-0.5714286*d > p_a));
  pa = !(pr == 1 || pg == 1);
  
  //pg = !(c_m < 10 || p_f < 0.7 || (1.75 -0.1*c_m > p_f) || (d < 0.2) || (p_a < 0.9));
  //pr = (c_m < 9 || p_f < 0.6 || (1.6 - 0.1*c_m > p_f) || (d < 0.0) || (p_a < 0.8));
  //pa = !(pr == 1 || pg == 1);
  
}