#!/usr/bin/env python

import pickle
import numpy
import pystan
import sys
import run

def main(data):
   nchain=4
   nsne=data['nsne']
   nbands=data['nbands']
   ntotspec=data['ntotspec']
   R_simplex = ((-1.)**numpy.arange(nsne)*.25 + .5)*2./nsne
   R_simplex = R_simplex/R_simplex.sum()

   init = [{'c_eta_sq_raw':1/1e3+numpy.zeros(nbands),
      'c_inv_rho_sq_raw':200./1e4+numpy.zeros(nbands),
      'c_sigma_sq':0.1+numpy.zeros(nbands),
      'Delta_scale':2,
      'Delta_unit':R_simplex,
      't_max_raw':numpy.zeros(nsne),
      'cfn': numpy.zeros((nbands,ntotspec)),
      'd1_unit':R_simplex,
      'd2_unit':R_simplex,
      # 'dust1_0':10.,
      'dust2_0':-1.,
      'dust1':numpy.arange(1.,nbands+1)/nbands,
      'dust2_rest':numpy.zeros(nbands-1)
      } for _ in range(nchain)]

   sm = pystan.StanModel(file='model2.stan')
   control = {'stepsize':1}
   fit = sm.sampling(data=data, iter=2000, chains=nchain,control=control,init=init, thin=1)

   output = open('temp3.pkl','wb')
   pickle.dump((fit.extract(),fit.get_sampler_params()), output, protocol=2)
   output.close()
   print fit


if __name__ == "__main__":
   nsne=10
   data=run.makedata(nsne)
   main(data)