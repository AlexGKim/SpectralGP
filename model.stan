#./gerard11 sample num_warmup=5000 num_samples=5000 data file=data.R init=init11.R output file=output11.csv refresh=1000

data {
  int nsne;
  int nbands;
  int maxphases;

  int nphases[nsne];

  matrix[maxphases, nbands] flux[nsne];
  matrix[maxphases, nbands] vars[nsne];
  vector[maxphases] phase[nsne];
}

transformed data {
  int dmatrix;
  vector[dmatrix] mu;
  for (i in 1:dmatrix) mu[i]=0;
  dmatrix=1;
  for (i in 1:nsne) dmatrix=dmatrix*nphases[i];
}

parameters {
  real<lower=0> c_eta_sq;
  real<lower=0> c_inv_rho_sq;
  real<lower=0> c_sigma_sq;

  real <lower=0> Delta_scale;
  simplex[nsne] Delta_unit;

  real t_max[nsne];
}

transformed parameters {
  vector[nsne] Delta;
  real<lower=0> c_rho_sq;
  Delta = Delta_scale*(Delta_unit-1./nsne);
  c_rho_sq = inv(c_inv_rho_sq);
}

model {

  matrix[dmatrix, dmatrix] Sigma;
  vector[dmatrix] y;
  int matind1;
  int matind2;

  for (k1 in 1:nbands){
    #treat each band indepedently
    matind1=1;
    for (i1 in 1:nsne){
      for (j1 in 1:nphases[i1]){
        y[matind1]=flux[i1,j1,k1];
        matind2=1;
        for (i2 in 1:nsne){
          for (j2 in 1:nphases[i2]){
            if (matind1 == matind2){
              Sigma[matind1,matind2]= c_eta_sq + c_sigma_sq+vars[i1,j1,k1];
            } else {
              Sigma[matind1,matind2]= c_eta_sq +c_eta_sq * exp(-c_rho_sq * pow(phase[i1,j1]-t_max[i1] - phase[i2,j2]+t_max[i2],2));
            }
            matind2=matind2+1;
          }
        }
        matind1=matind1+1;
      }
    }
    y ~ multi_normal(mu, Sigma);
  }
  c_eta_sq ~ cauchy(0, 5);
  c_inv_rho_sq ~ cauchy(0, 5);
  c_sigma_sq ~ cauchy(0, 5);     
}