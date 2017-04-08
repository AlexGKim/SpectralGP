#!/usr/bin/env python

import pickle
import matplotlib.pyplot as plt
import corner
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl
import numpy
import sys
import run

# mpl.rcParams['font.size'] = 18

def main(argv, data):
    outdir='output{}/'.format(argv)

    ntotspec=data['ntotspec']
    nbands=data['nbands']

    f = open('temp{}.pkl'.format(argv),'rb')
    (fit,_) = pickle.load(f)

    figure=corner.corner(fit['t_max'])
    plt.savefig(outdir+'t_max_corner.pdf')

    for i in xrange(nbands):
        mega=numpy.array([fit['c_eta_sq'][:,i],fit['c_inv_rho_sq'][:,i],fit['c_sigma_sq'][:,i]])
        mega=numpy.transpose(mega)
        figure=corner.corner(mega) #,labels=[r'$\eta^2$',r'$w$',r'$\sigma$'])
        plt.savefig(outdir+'param_corner{}.pdf'.format(i))

    cumsum = numpy.cumsum(data['nphases'])
    cumsum=numpy.concatenate(([0],cumsum))
    snin, pin = run.mastertosnspec(xrange(ntotspec),cumsum)
    absphase = data['phase'][None,:]-fit['t_max'][:,snin]
    relphase = absphase-absphase[:,0][:,None]      
    plt.clf()
    for i in xrange(nbands):
        sub = numpy.random.randint(0,ntotspec*fit['t_max'].shape[0],2000)
        plt.scatter(relphase.flatten()[sub], fit['cfn'][:,i,:].flatten()[sub]+i, \
            label='{:}+{:}'.format(i,i*0.5), marker='.',s=2)
    plt.legend()
    plt.gca().invert_yaxis()
    plt.savefig(outdir+'lc.pdf')

if __name__ == "__main__":
    data = run.makedata()
   main(sys.argv[1],data)