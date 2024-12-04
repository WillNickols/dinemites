data {
    int<lower=0> N; // Number of person-day-alleles
    int<lower=1> J; // Number of alleles
    int<lower=0> n_times; // Number of days
    int<lower=1> C; // Number of clusters
    
    array[N] int<lower=1,upper=J> alleles; // Alleles
    array[N] int<lower=0,upper=1> y; // Allele present
    array[N] int<lower=1,upper=n_times> times; // Times
    
    matrix[J, C] alleles_in_clusters;
}

parameters {
    matrix<lower=0, upper=1>[n_times, C] prob_cluster_in_time;
    matrix<lower=0, upper=1>[J, C] prob_allele_in_cluster;
}

model {
    for (n_time in 1:n_times) {
        for (c in 1:C) {
            prob_cluster_in_time[n_time, c] ~ beta(1.1,1.1);
        }
    }
    
    for (j in 1:J) {
        for (c in 1:C) {
            prob_allele_in_cluster[j, c] ~ beta(1.1,1.1);
        }
    }

    for (n in 1:N) {
        real prob_current = 1;
        // Iterate through clusters the allele could be in
        for (c in 1:C) {
            if (alleles_in_clusters[alleles[n], c] == 1) {
                prob_current = prob_current * (1 - prob_cluster_in_time[times[n], c] * prob_allele_in_cluster[alleles[n], c]);
            }
        }
        prob_current = 1 - prob_current;
        y[n] ~ bernoulli(prob_current);
    }
}
