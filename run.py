#!/usr/bin/env python

import pickle
import numpy
import pystan
import sys

# numpy.seterr(all='raise')

nchain=4

pkl_file = open('data.pkl', 'r')
(allphases, allfluxes, allvaris, ledge) = pickle.load(pkl_file)
pkl_file.close()

nsne = len(allphases)
nsne=4
nphases = numpy.zeros(nsne,dtype='int')
for i in xrange(nsne):
   nphases[i]=len(allphases[i]) 
maxphases = nphases.max()
nbands = allfluxes[0][0].shape[0]
# print nsne, maxphases, nbands

flux = numpy.zeros((nsne,maxphases,nbands))
var = numpy.zeros((nsne,maxphases,nbands))
phases = numpy.zeros((nsne,maxphases))

runnum=0.
for i in xrange(nsne):
   for j in xrange(len(allphases[i])):
      flux[i,j,:]=allfluxes[i][j]
      var[i,j,:]=allvaris[i][j]
      runnum = runnum+len(allfluxes[i][j])
   phases[i][:len(allphases[i])]=allphases[i]
av=flux.sum()/runnum
flux=flux/av
var = var/av**2

ind4400 = numpy.searchsorted(ledge, 4400.,side='right')  #note that the index used by STAN is 1-based
data = {'nsne': nsne, 'nbands': nbands, 'maxphases': maxphases, 'nphases': nphases, 'ind4400': ind4400, 'flux': flux, 'vars': var, 'phase':phases}

R_simplex = ((-1.)**numpy.arange(nsne)*.25 + .5)*2./nsne
R_simplex = R_simplex/R_simplex.sum()

init = [{'c_eta_sq':0.02,
   'c_inv_rho_sq':2*16.,
   'c_sigma_sq':1e-4,
   'Delta_scale':0.1,
   'Delta_unit':R_simplex,
   't_max':numpy.mean(phases,axis=1),
   'cfn': numpy.zeros((nsne,maxphases,nbands))
   } for _ in range(nchain)]

sm = pystan.StanModel(file='model.stan')
control = {'stepsize':1}
fit = sm.sampling(data=data, iter=2000, chains=nchain,control=control,init=init, thin=1)


output = open('temp.pkl','wb')
pickle.dump((fit.extract(),fit.get_sampler_params()), output, protocol=2)
output.close()
# print fit