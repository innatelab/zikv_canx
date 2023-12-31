data {
  int<lower=1> Nmsruns;
  int<lower=1> NeffectPairs;
  real<lower=0, upper=1> outlierProb;

  vector[NeffectPairs] effect_x;
  vector[NeffectPairs] effect_y;
  int<lower=1, upper=Nmsruns> msrun_x[NeffectPairs];
  int<lower=1, upper=Nmsruns> msrun_y[NeffectPairs];

  real<lower=0> sigma_t_offset;
  real<lower=0> scale_sigma_df;
  real<lower=0> data_sigma_df;
}

transformed data {
    int<lower=0> NmsrunPairs = Nmsruns * Nmsruns;
    int<lower=1, upper=Nmsruns> msrunpair2x[NmsrunPairs];
    int<lower=1, upper=Nmsruns> msrunpair2y[NmsrunPairs];
    int<lower=1, upper=NmsrunPairs> effectpair2msrunpair[NeffectPairs];

    for (y in 1:Nmsruns) {
      for (x in 1:Nmsruns) {
        msrunpair2x[(y - 1) * Nmsruns + x] = x;
        msrunpair2y[(y - 1) * Nmsruns + x] = y;
      }
    }
    for (i in 1:NeffectPairs) {
      effectpair2msrunpair[i] = (msrun_y[i] - 1) * Nmsruns + msrun_x[i];
    }
}

parameters {
    real<lower=0> data_sigma_a;
    real<lower=0> data_sigma_t;
    real<lower=0> scale_sigma_a;
    real<lower=0> scale_sigma_t;
    vector<lower=-5,upper=5>[Nmsruns] scale_log2;
}

transformed parameters {
    real<lower=0> data_sigma = data_sigma_a * sqrt(data_sigma_t);
    real<lower=0> scale_sigma = scale_sigma_a * sqrt(scale_sigma_t);
}

model {
    vector[NmsrunPairs] scale_ratio_xy = exp2(scale_log2[msrunpair2x] - scale_log2[msrunpair2y]);

    data_sigma_t - sigma_t_offset ~ inv_gamma(0.5 * data_sigma_df, 0.5 * data_sigma_df);
    data_sigma_a ~ std_normal();

    scale_sigma_t - sigma_t_offset ~ inv_gamma(0.5 * scale_sigma_df, 0.5 * scale_sigma_df);
    scale_sigma_a ~ std_normal();

    scale_log2 ~ normal(0.0, scale_sigma);

    for (i in 1:NeffectPairs) {
      real delta = effect_x[i] - scale_ratio_xy[effectpair2msrunpair[i]] * effect_y[i];
      target += log_mix(outlierProb,
                        cauchy_lpdf(delta | 0, data_sigma),
                        normal_lpdf(delta | 0, data_sigma));
    }
}

generated quantities {
    vector<lower=0>[Nmsruns] scale_mult = exp2(scale_log2);
}
