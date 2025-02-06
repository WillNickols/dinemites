data {
    int<lower=0> N; // Number of person-day-alleles
    int<lower=0> K_infection_general;
    int<lower=0> K_infection_persistent; // Number of persistent covariates for infections
    int<lower=1> K_alleles_persistent; // Number of persistent covariates for alleles
    int<lower=1> J; // Number of alleles

    matrix[N, K_infection_general] X_infection_general;
    matrix[N, K_infection_persistent] X_infection_persistent;
    matrix[N, K_alleles_persistent] X_alleles_persistent;

    array[N] int<lower=1,upper=J> group; // Alleles
    array[N] int<lower=0> t_gap; // Time gap

    array[N] int<lower=0, upper=1> z; // Infection present
    array[N] int<lower=0, upper=1> y; // Allele present
    array[N] int<lower=0> w; // Weights
}

parameters {
    real mu_alpha_alleles_new;
    real<lower=0> sigma_alpha_alleles_new;
    real<lower=0> a_drop_out;
    real<lower=0> b_drop_out;
    vector[K_alleles_persistent] mu_alleles_old;
    vector<lower=0>[K_alleles_persistent] sigma_alleles_old;

    real alpha_infection;
    vector[K_infection_general] beta_infection_general;
    vector[K_infection_persistent] beta_infection_persistent;

    vector[J] alpha_alleles_new;
    vector<lower=0, upper=1>[J] drop_out;
    matrix[J, K_alleles_persistent] beta_alleles_old;
}

model {
    mu_alpha_alleles_new ~ normal(-3, 10);
    sigma_alpha_alleles_new ~ cauchy(0, 10);
    a_drop_out ~ exponential(1);
    b_drop_out ~ exponential(1);
    mu_alleles_old ~ normal(0, 10);
    sigma_alleles_old ~ cauchy(0, 10);

    alpha_infection ~ normal(-3, 10);
    beta_infection_general ~ normal(0, 10);
    beta_infection_persistent ~ normal(0, 10);

    alpha_alleles_new ~ normal(mu_alpha_alleles_new, sigma_alpha_alleles_new);
    drop_out ~ beta(a_drop_out, b_drop_out);

    for (j in 1:J) {
        for (k in 1:K_alleles_persistent) {
            beta_alleles_old[j, k] ~ normal(mu_alleles_old[k], sigma_alleles_old[k]);
        }
    }

    real intermediate_likelihood = 0;
    for (n in 1:N) {
        // Overall infection piece
        real per_day_rate = exp(alpha_infection + dot_product(X_infection_general[n], beta_infection_general));
        real p_new_infection = 1 - exp(-t_gap[n] * per_day_rate);
        real p_old_infection;
        if (sum(abs(X_infection_persistent[n])) > 0) {
            p_old_infection = inv_logit(dot_product(X_infection_persistent[n], beta_infection_persistent));
        } else {
            p_old_infection = 0;
        }
        real p_any_infection = p_new_infection + p_old_infection - p_new_infection * p_old_infection;

        intermediate_likelihood += w[n] * bernoulli_logit_lpmf(z[n] | p_any_infection);

        if (z[n] == 1) {
            // Allele-related probabilities
            real p_allele_if_new = inv_logit(alpha_alleles_new[group[n]]);
            real p_allele_if_old;

            if (sum(abs(X_alleles_persistent[n])) > 0) {
                p_allele_if_old = inv_logit(dot_product(X_alleles_persistent[n], beta_alleles_old[group[n]]));
            } else {
                p_allele_if_old = inv_logit(alpha_alleles_new[group[n]]) * drop_out[group[n]];
            }

            // Final allele infection probabilities
            real p_allele_new = p_allele_if_new * p_new_infection;
            real p_allele_old = p_allele_if_old * p_old_infection;

            // Likelihood for y[n]
            target += w[n] * bernoulli_lpmf(y[n] | (p_allele_new + p_allele_old - p_allele_new * p_allele_old) / p_any_infection);
        }
    }
    intermediate_likelihood = intermediate_likelihood / J;
    target += intermediate_likelihood;
}

