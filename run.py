#!/usr/bin/env python

import pickle
import numpy
import pystan
import sys

def mastertosnspec(m, cumsum):
   index = numpy.searchsorted(cumsum, m,side='right')-1
   return (index,m-cumsum[index])

def makedata(nsne):
   if nsne is None:
      nsne = len(allphases)
   pkl_file = open('data.pkl', 'r')
   (allphases, allfluxes, allvaris, ledge) = pickle.load(pkl_file)
   pkl_file.close()

   nphases = numpy.zeros(nsne,dtype='int')
   for i in xrange(nsne):
      nphases[i]=len(allphases[i])

   cumsum = numpy.cumsum(nphases)
   cumsum=numpy.concatenate(([0],cumsum))

   ntotspec = nphases.sum()

   nbands = allfluxes[0][0].shape[0]

   flux = numpy.zeros((nbands,ntotspec))
   sig = numpy.array(flux)
   phase = numpy.zeros(ntotspec)

   for i in xrange(nbands):
      index = 0
      for j in xrange(nsne):
         for k in xrange(nphases[j]):
            flux[i,index]=allfluxes[j][k][i]
            sig[i,index]=allvaris[j][k][i]
            phase[index]=allphases[j][k]
            index +=1

   sig = numpy.sqrt(sig)
   avflux = flux.mean()
   flux=flux/avflux
   sig=sig/avflux

   ind4400 = numpy.searchsorted(ledge, 4400.,side='right')  #note that the index used by STAN is 1-based
   data = {'nsne': nsne, 'nbands': nbands, 'nphases': nphases, 'ntotspec':ntotspec,'ind4400': ind4400, 'flux': flux, 'sig':sig, 'phase':phase}
   return data

def main(data):
   nchain=4
   nsne=data['nsne']
   nbands=data['nbands']
   ntotspec=data['ntotspec']
   R_simplex = ((-1.)**numpy.arange(nsne)*.25 + .5)*2./nsne
   R_simplex = R_simplex/R_simplex.sum()

   init = [{'c_eta_sq':1+numpy.zeros(nbands),
      'c_inv_rho_sq':200+numpy.zeros(nbands),
      'c_sigma_sq':0.01+numpy.zeros(nbands),
      'Delta_scale':5,
      'Delta_unit':R_simplex,
      't_max':numpy.zeros(nsne),
      'cfn': numpy.zeros((nbands,ntotspec))
      } for _ in range(nchain)]

   sm = pystan.StanModel(file='model.stan')
   control = {'stepsize':1}
   fit = sm.sampling(data=data, iter=1000, chains=nchain,control=control,init=init, thin=1)


   output = open('test1.pkl','wb')
   pickle.dump((fit.extract(),fit.get_sampler_params()), output, protocol=2)
   output.close()
   print fit

if __name__ == "__main__":
   nsne=10
   data=makedata(nsne)
   main(data)