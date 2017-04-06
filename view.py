#!/usr/bin/env python

import pickle
import matplotlib.pyplot as plt
import corner
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl
import numpy
mpl.rcParams['font.size'] = 18

f = open('temp.pkl','rb')
(fit,_) = pickle.load(f)

mega=numpy.array([fit['Delta_scale'],fit['c_eta_sq'],fit['c_inv_rho_sq'],fit['c_sigma_sq']])
mega=numpy.transpose(mega)
figure=corner.corner(mega) #,labels=[r'$\eta^2$',r'$w$',r'$\sigma$'])
plt.show()


figure=corner.corner(fit['t_max'])
plt.show()

wefwe

pkl_file = open('data.pkl', 'r')
(allphases, allfluxes, allvaris, ledge) = pickle.load(pkl_file)
pkl_file.close()

nsne = len(allphases)
nsne=4
nbands = allfluxes[0][0].shape[0]

for i in xrange(nbands):
    for j in xrange(nsne):
        for dum in xrange(len(fit['Delta_scale'])):
            plt.scatter(allphases[j]-fit['t_max'][dum,j],fit['cfn'][dum,j,0:len(allphases[j]),i]+nbands*0.5)
plt.show()