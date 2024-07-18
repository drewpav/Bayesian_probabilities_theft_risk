functions {
    real icar_normal_lpdf(vector phi, int N, array[] int node1, array[] int node2) {
        return -0.5 * dot_self(phi[node1] - phi[node2]);
    }
}
data {
    int<lower=0> N;
    int<lower=0> k;
    int<lower=0> N_edges;
    array[N_edges] int<lower=1, upper=N> node1;
    array[N_edges] int<lower=1, upper=N> node2;
    array[N] int<lower=0> Y;                                        // the dependent variable here is counts of theft
    matrix[N, k] X;                                                // the independent variables here are the IMD deprivation score and house prices
    vector<lower=0>[N] E;                                           // estimated number of expected counts of theft
}
transformed data {
    vector[N] log_offset = log(E);                                  // the expected incidents of theft are used as an offset and added to the regression model
}
parameters {
    real alpha;                                                     // defining the intercept (overall risk in population)
    vector[k] beta;                                                   // defining the coefficient for both the independent variable
    real<lower=0> sigma;                                            // defining the overall standard deviation producted with spatial effect smoothing term phi
    vector[N] phi;                                                  // spatial effect smoothing term or spatial ICAR component of the model 
}
model {
    phi ~ icar_normal(N, node1, node2);                             // prior for the spatial random effects
    Y ~ poisson_log(log_offset + alpha + X*beta + phi*sigma);       // likelihood function i.e., spatial ICAR model using Possion distribution
    alpha ~ normal(0.0, 1.0);                                       // prior for intercept   (weak/uninformative prior)
    beta ~ normal(0.0, 1.0);                                        // prior for coefficient (weak/uninformative prior)
    sigma ~ normal(0.0, 1.0);                                       // prior for SD          (weak/uninformative prior)
    sum(phi) ~ normal(0, 0.001*N);
}
generated quantities {
    vector[N] eta = alpha + X*beta + phi*sigma;                     
    vector[N] mu = exp(eta);                                        
}
