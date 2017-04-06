#./gerard11 sample num_warmup=5000 num_samples=5000 data file=data.R init=init11.R output file=output11.csv refresh=1000

data {
  int nsne;
  int nbands;
  int maxphases;

  int nphases[nsne];
  int ind4400;

  matrix[maxphases, nbands] flux[nsne];
  matrix[maxphases, nbands] vars[nsne];
  vector[maxphases] phase[nsne];
}

transformed data {
  int ntotspec;
  ntotspec=0;
  for (i in 1:nsne) ntotspec=ntotspec+nphases[i];  
}

parameters {
  real<lower=0> c_eta_sq;
  real<lower=0> c_inv_rho_sq;
  real<lower=0> c_sigma_sq;

  real <lower=0> Delta_scale;
  simplex[nsne] Delta_unit;

  real t_max[nsne];

  matrix[maxphases, nbands] cfn[nsne];
}

transformed parameters {
  vector[nsne] Delta;
  real<lower=0> c_rho_sq;
  Delta = Delta_scale*(Delta_unit-1./nsne);
  c_rho_sq = inv(c_inv_rho_sq);
}

model {

  int dim;
  int matind1;
  int matind2;

  vector[ntotspec] y0;
  vector[ntotspec] sig0;

  #treat each band indepedently
  for (k1 in 1:nbands){
    #assign vectors of the fluxes and uncertainties for data only in this band   
    matind1=1;
    for (i1 in 1:nsne){
      for (j1 in 1:nphases[i1]){
        y0[matind1]=flux[i1,j1,k1];
        sig0[matind1]=sqrt(vars[i1,j1,k1]);
        matind1=matind1+1;
      }
    }

    if (k1==ind4400){
      dim = ntotspec+1;
    }
    else {
      dim = ntotspec;
    }

    {
      vector[dim] zerovec;
      matrix[dim, dim] Sigma;
      vector[dim] c;
      vector[ntotspec] y_trunc;
      vector[ntotspec] D_;
      for (i in 1:dim) zerovec[i]=0;

      matind1=1;
      for (i1 in 1:nsne){
        for (j1 in 1:nphases[i1]){
          c[matind1]=cfn[i1,j1,k1];
          D_[matind1]=Delta[i1];

          #fill in data
          matind2=1;
          for (i2 in 1:nsne){
            for (j2 in 1:nphases[i2]){
              if (matind1 <= matind2){
                if (matind1 == matind2){
                  Sigma[matind1,matind2]= c_eta_sq + c_sigma_sq;
                  } else {
                    Sigma[matind1,matind2]= c_eta_sq * exp(-c_rho_sq * pow(phase[i1,j1]-t_max[i1] - phase[i2,j2]+t_max[i2],2));
                    Sigma[matind2,matind1]=Sigma[matind1,matind2];
                  }
                }
                matind2=matind2+1;
              }
            }

            #if 4400A add the data -- slope covariance
            if (k1==ind4400){
              Sigma[matind1,matind2]= -2 * c_rho_sq * (phase[i1,j1]-t_max[i1]) * c_eta_sq * exp(-c_rho_sq * pow(phase[i1,j1]-t_max[i1],2));
              Sigma[matind2,matind1]=Sigma[matind1,matind2];
            }
            matind1=matind1+1;
          }
        }

        # if 4400A add the zero slope condition
        if (k1==ind4400){
          Sigma[matind1,matind1]= 2 * c_rho_sq * c_eta_sq;
          c[matind1]=0;
        }

        c ~ multi_normal(zerovec, Sigma);
        for (dum in 1:ntotspec) y_trunc[dum]=10.^(-(D_[dum]+c[dum])/2.5);
        y0 ~ normal(y_trunc,sig0);
      }
    }
    c_eta_sq ~ cauchy(0, 5);
    c_inv_rho_sq ~ cauchy(0, 100);
    c_sigma_sq ~ cauchy(0, 5);     
  }