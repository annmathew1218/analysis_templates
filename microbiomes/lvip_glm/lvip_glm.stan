#include functions.stan
data {
    int<lower=2> NT;                        // number of microbial tips
    int<lower=1> NI;                        // number of microbial internal nodes
    vector<lower=0>[NT + NI] time;      // time between adjacent microbial tips/nodes. standardize so mean time is 1
    array[NT + NI] int self;                // index for each microbial tip/node in a single vector
    array[NT + NI] int ancestor;            // index for each microbial tip/node's ancestor in a single vector
    int NS;                                 // number of samples
    int NB_s;                               // number of factor levels
    int NSB;                                // number of factors
    array[NB_s] int idx;                    // mapping of sigmas to factor levels
    matrix[NS,NB_s] X_s;                    // model matrix for samples (e.g. tissue compartments, duplicates, sequencing depth, etc.). must include intercept
    array[NT,NS] int count;                 // observations
    real inv_log_max_contam;                       // prior expectation of contamination rate
    real<lower=0> shape_gnorm;                          // strength of prior pulling contamination toward zero

}
transformed data {
    int NN = NT + NI;
    vector<lower=0>[NN] time_sqrt = sqrt(time);          // time between adjacent microbial tips/nodes. standardize so mean time is 1
    vector[NN] time_log = log(time);
    array[NN] int self_i;
    array[NN] int ancestor_i;
    array[NI] int self_i2;
    for(m in 1:NN){
        self_i[m] = self[m] - NT;
        ancestor_i[m] = ancestor[m] - NT;
        if(self_i[m] > 0) {
            self_i2[self_i[m]] = m;
        }
    }
}
parameters {
    real theta_prevalence;           // traits effects on default logit probability of association
    real theta_abundance;            // traits effects on default strength of association
    vector<lower=0>[NSB] sd_prevalence; // variance of sample effects
    vector<lower=0>[NSB+1] sd_abundance;  // variance of sample effects
    vector<lower=0>[NSB] sigma_prevalence; // rate of evolution of sample effects
    vector<lower=0>[NSB] sigma_abundance;  //  rate of evolution of sample effects
    matrix[NB_s,NN] delta_prevalence;
    matrix[NB_s,NN] delta_abundance;
    matrix[NS,NT] abundance_observed;
    vector[NS] multinomial_nuisance;
    real<upper=0> inv_log_less_contamination; // smaller = less average contamination
    real<lower=0> contaminant_overdisp;            // dispersion parameter for amount of contamination in true negative count observations
}

transformed parameters {
    real log_less_contamination = inv(inv_log_less_contamination);
    vector[NSB] alpha_prevalence_log = 2*(sigma_prevalence - sd_prevalence) - log2(); // OU alpha_prevalence is function of total variance for each factor () and the rate of evolution at each branch ()
    vector[NSB] alpha_abundance_log = 2*(sigma_abundance - sd_abundance[1:NSB]) - log2(); // OU alpha_abundance is function of total variance for each factor and the rate of evolution at each branch
    vector[NSB] alpha_prevalence = exp(alpha_prevalence_log);
    vector[NSB] alpha_abundance = exp(alpha_abundance_log);
    matrix[NB_s,NN] beta_prevalence;
    matrix[NB_s,NN] beta_abundance;
    beta_prevalence[,self[1]] = theta_prevalence + sigma_prevalence[idx] * time_sqrt[self[1]] .* delta_prevalence[,self[1]];
    beta_abundance[,self[1]] = theta_abundance + sigma_abundance[idx] * time_sqrt[self[1]] .* delta_abundance[,self[1]];
    for(m in 2:NN) {
       vector[NSB] prop_prevalence = exp(-exp(alpha_prevalence_log + time_log[self[m]])); // OU process takes time-and-sigma-weighted mean state
       vector[NSB] prop_abundance = exp(-exp(alpha_abundance_log + time_log[self[m]])); // OU process takes time-and-sigma-weighted mean state
       beta_prevalence[,self[m]]
           = prop_prevalence[idx] .* beta_prevalence[,ancestor[m]]
             + (1 - prop_prevalence[idx]) * theta_prevalence
             + sigma_prevalence[idx] * time_sqrt[self[m]] .* delta_prevalence[,self[m]];
       beta_abundance[,self[m]]
           = prop_abundance[idx] .* beta_abundance[,ancestor[m]]
             + (1 - prop_abundance[idx]) * theta_abundance
             + sigma_abundance[idx] * time_sqrt[self[m]] .* delta_abundance[,self[m]];
    }
}
model {
    matrix[NS,NN] prevalence = X_s * beta_prevalence;
    matrix[NS,NN] abundance_predicted = X_s * beta_abundance;
    matrix[NS,NN] abundance_contam
      = rep_matrix(beta_abundance[1,]
                   + log_inv_logit(beta_prevalence[1,])
                   + log_less_contamination,
                   NS);
    target += normal_lpdf(theta_prevalence | 0,5);
    target += normal_lpdf(theta_abundance | 0,5);
    target += normal_lpdf(sd_prevalence | 0,2.5);
    target += normal_lpdf(sd_abundance | 0,2.5);
    target += lognormal_lpdf(sigma_prevalence | log(sd_prevalence), 0.5);
    target += lognormal_lpdf(sigma_abundance | log(sd_abundance[1:NSB]), 0.5);
    target += std_normal_lpdf(to_vector(delta_prevalence)); // fixed rate of evolution
    target += std_normal_lpdf(to_vector(delta_abundance)); // fixed rate of evolution
    target += generalized_std_normal_1_lpdf(inv_log_less_contamination / inv_log_max_contam | shape_gnorm);   // shrink amount of contamination in 'true zeros' toward zero
    target += lognormal_lpdf(contaminant_overdisp | 0, 0.1);                                               // shrink overdispersion of contaminant counts in 'true zeros' toward zero
    for(m in 1:NT) {
        for(s in 1:NS) {
            target += log_sum_exp(log1m_inv_logit(prevalence[s,self[m]])
                                  + normal_lpdf(abundance_observed[s,m] |
                                                abundance_contam[s,self[m]],
                                                contaminant_overdisp * sd_abundance[NSB+1]), //estimated abundance if true negative
                                  log_inv_logit(prevalence[s,self[m]])
                                  + normal_lpdf(abundance_observed[s,m] |
                                                log_sum_exp(abundance_contam[s,self[m]], abundance_predicted[s,self[m]]),
                                                sd_abundance[NSB+1])); //estimated abundance if true positive
        }
    }
    target += poisson_log_lpmf(to_array_1d(count) |
                               to_vector(abundance_observed + rep_matrix(multinomial_nuisance, NT)));
}
