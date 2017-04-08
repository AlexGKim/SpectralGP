functions {
  int[] mastertosnspec(int m, int[] cumsum){
    int i;
    int ans[2];
    i=1;
    while (m > cumsum[i]){
      i=i+1;
    }
    ans[1]= i;
    if (i==1) {
      ans[2] = m;
    }
    else{
      ans[2]= m-cumsum[i-1];
    }
    return ans;
  }
}

data {
  int nsne;
  int nbands;
  int nphases[nsne];

  # int maxphases;
  int ntotspec;
  int ind4400;

  vector[ntotspec] flux[nbands];
  vector[ntotspec] sig[nbands];

  vector[ntotspec] phase;
}

transformed data {
  int cumsum[nsne];
  vector[ntotspec] zerovec1;
  vector[ntotspec+1] zerovec2;

  cumsum[1] = nphases[1];
  for (i in 2:nsne) cumsum[i] = cumsum[i-1]+nphases[i];

  for (i in 1:ntotspec) zerovec1[i]=0;
  for (i in 1:ntotspec+1) zerovec2[i]=0;
}

parameters {
  real<lower=0> c_eta_sq[nbands];
  real<lower=0> c_inv_rho_sq[nbands];
  real<lower=0> c_sigma_sq[nbands];

  real <lower=0> Delta_scale;
  simplex[nsne] Delta_unit;

  real t_max[nsne];
  vector[ntotspec] cfn[nbands];

  real<lower=0> dust1_0;
  vector[nbands-1] dust1_rest;

  real<lower=0> dust2_0;
  vector[nbands-1] dust2_rest;

  simplex[nsne] d1_unit;
  simplex[nsne] d2_unit;
}

transformed parameters {
  vector[nsne] Delta;
  real<lower=0> c_rho_sq[nbands];

  vector[nbands] dust1;
  vector[nbands] dust2;

  vector[nsne] d1;
  vector[nsne] d2;

  Delta = Delta_scale*(Delta_unit-1./nsne);
  for (i in 1:nbands) c_rho_sq [i]= inv(c_inv_rho_sq[i]);

  dust1[1]=dust1_0;
  dust2[1]=dust2_0;
  for (i in 2:nbands){
    dust1[i] = dust1_rest[i-1];
    dust2[i] = dust2_rest[i-1];
  }
  d1=(d1_unit-1./nsne);
  d2=(d2_unit-1./nsne);
}

model {

  int dim;
  int in1[2];
  int in2[2];
  vector[ntotspec] D_;
  vector[ntotspec] f_model;

  vector[ntotspec] d1_;
  vector[ntotspec] d2_;

  # convenience vectors for SN-based values
  for (index in 1:ntotspec){
    in1=mastertosnspec(index, cumsum);
    D_[index] = Delta[in1[1]];
    d1_[index] = d1[in1[1]];
    d2_[index] = d2[in1[1]];
  }

  #treat each band indepedently
  for (k1 in 1:nbands){

    # the dimensionality of the matrix
    if (k1==ind4400){
      dim = ntotspec+1;
    }
    else {
      dim = ntotspec;
    }

    {
      # fill the covariance matrix
      matrix[dim, dim] Sigma;
      for (i1 in 1:ntotspec){
        in1 = mastertosnspec(i1, cumsum);
        Sigma[i1,i1]= c_eta_sq[k1] + c_sigma_sq[k1];

        for (i2 in i1+1:ntotspec){
          in2 = mastertosnspec(i2, cumsum);
          Sigma[i1,i2]= c_eta_sq[k1] * exp(-c_rho_sq[k1] * pow(phase[i1]-t_max[in1[1]] - phase[i2]+t_max[in2[1]],2));
          Sigma[i2,i1]=Sigma[i1,i2];
        }
        #if it is 4400 add in extra row/column
        if (k1==ind4400){
          Sigma[i1,ntotspec+1]= 2 * c_rho_sq[k1] * (phase[i1]-t_max[in1[1]]) * c_eta_sq[k1] * exp(-c_rho_sq[k1] * pow(phase[i1]-t_max[in1[1]],2));
          Sigma[ntotspec+1,i1] = Sigma[i1,ntotspec+1];
        }
      }
      if (k1==ind4400) Sigma[ntotspec+1,ntotspec+1] = 2 * c_rho_sq[k1] * c_eta_sq[k1];

      #calculate the c GP
      if (k1 != ind4400){
        cfn[k1] ~ multi_normal(zerovec1, Sigma);
      } 
      else{
        {
          vector[dim] c;
          for (dum in 1:ntotspec) c[dum]=cfn[k1][dum];
          c[ntotspec+1]=0;
          c ~ multi_normal(zerovec2, Sigma);
        }
      }           
    }

    # probability for data
    f_model = D_ + cfn[k1] + d1_*dust1[k1] + d2_*dust2[k1];
    for (dum in 1:ntotspec) f_model[dum]=10.^(-f_model[dum]/2.5);
    flux[k1] ~ normal(f_model, sig[k1]);
    c_eta_sq[k1] ~ cauchy(0, 10);
    c_inv_rho_sq[k1] ~ cauchy(0, 400);
    c_sigma_sq[k1] ~ cauchy(0, 1);  
  }
   
}