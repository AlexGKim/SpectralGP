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

mega=numpy.array([fit['Delta_scale'],fit['c_eta_sq'],numpy.sqrt(fit['c_inv_rho_sq']/2),numpy.sqrt(fit['c_sigma_sq'])])
mega=numpy.transpose(mega)
print mega.shape
figure=corner.corner(mega)
plt.show()

pkl_file = open('data.pkl', 'r')
(allphases, allfluxes, allvaris, ledge) = pickle.load(pkl_file)
pkl_file.close()

nsne = len(allphases)
nsne=4
nbands = allfluxes[0][0].shape[0]

for i in xrange(nbands):
    for j in xrange(nsne):
        for dum in xrange(len(fit['Delta_scale'])):
            plt.scatter(allphases[j]-fit['t_max'][j,dum],fit['cfn'][j,0:len(allphases[j]),dum])
    plt.show()