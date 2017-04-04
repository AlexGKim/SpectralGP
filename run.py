#!/usr/bin/env python

import pickle
import numpy
import pystan


pkl_file = open('data.pkl', 'r')
(allphases, allfluxes, allvaris) = pickle.load(pkl_file)
pkl_file.close()

nsne = len(allphases)
nphases = numpy.zeros(nsne,dtype='int')
for i in xrange(nsne):
   nphases[i]=len(allphases[i]) 
maxphases = nphases.max()
nbands = allfluxes[0][0].shape[0]
print nsne, maxphases, nbands

flux = numpy.zeros((nsne,maxphases,nbands))
var = numpy.zeros((nsne,maxphases,nbands))
phases = numpy.zeros((nsne,maxphases))
for i in xrange(nsne):
   for j in xrange(len(allphases[i])):
      flux[i,j,:]=allfluxes[i][j]
      var[i,j,:]=allvaris[i][j]
   phases[i][:len(allphases[i])]=allphases[i]

data = {'nsne': nsne, 'nbands': nbands, 'nphases': nphases, 'maxphases': maxphases,'flux': flux, 'vars': var, 'phase':phases}

R_simplex = ((-1.)**numpy.arange(nsne)*.25 + .5)*2./nsne
R_simplex = R_simplex/R_simplex.sum()
init = [{'c_eta_sq':0.01,
   'c_inv_rho_sq':0.01,
   'c_sigma_sq':1,
   'Delta_scale':0.1,
   'Delta_unit':R_simplex,
   't_max':numpy.mean(phases,axis=1)
   } for _ in range(4)]

sm = pystan.StanModel(file='model.stan')
control = {'stepsize':1}
fit = sm.sampling(data=data, iter=2000, chains=4,control=control,init=init, thin=1)


output = open('temp.pkl'.format(seed),'wb')
pickle.dump((fit.extract(),fit.get_sampler_params()), output, protocol=2)
output.close()
print fit