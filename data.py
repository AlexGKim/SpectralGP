#!/usr/bin/env python

import pickle
import cPickle
import numpy
from astropy.io import fits
from specutils.io import read_fits


prange = [-10,15]
lmin = 3300.
lmax = 9400.
nlbin = 10
ledge  = numpy.arange(numpy.log(lmin),numpy.log(lmax)+0.001,numpy.log(lmax/lmin)/nlbin)
ledge = numpy.exp(ledge)

# two color parameter model
dic_meta=cPickle.load(open("SNF-0203-CABALLO2/META.pkl"))

allphases=[]
allfluxes=[]
allvaris = []

for sn in dic_meta.itervalues():
    if sn['idr.subset'] != 'auxiliary':
        # container for all spectra
        phases=[]
        fluxes=[]
        varis=[]
        for spec in sn['spectra'].itervalues():
            # salt2.phase is in the supernova restframe
            if (spec['salt2.phase'] >= prange[0] and spec['salt2.phase'] < prange[1]):
                filename = 'SNF-0203-CABALLO2/'+spec['idr.spec_restframe']

                #get wavelength solution
                myspec = read_fits.read_fits_spectrum1d(filename)

                if myspec.dispersion[0] < lmin or myspec.dispersion[1] > lmax:
                    sedge = (myspec.dispersion+numpy.roll(myspec.dispersion,1))/2
                    sedge = sedge[1:]

                    hdulist = fits.open(filename)

                    #containers for all bands in a spectrum
                    f_=[]
                    v_=[]
                    for fin in xrange(len(ledge)-1):

                        insertloc = [numpy.searchsorted(sedge, ledge[fin],side='right'),numpy.searchsorted(sedge, ledge[fin+1],side='left')]
                        # print ledge[fin], sedge[insertloc[0]:insertloc[1]],ledge[fin+1]
                        intedge = numpy.concatenate(([ledge[fin]], sedge[insertloc[0]:insertloc[1]],[ledge[fin+1]]))
                        width = numpy.roll(intedge,-1)-intedge
                        # print myspec.dispersion[insertloc[1]-1],myspec.dispersion[insertloc[1]],myspec.dispersion[insertloc[1]+1]
                        # print intedge
                        # print width
                        # print len(width),insertloc[1]-insertloc[0]
                        f_.append(numpy.sum(width[:-1]*hdulist[0].data[insertloc[0]:insertloc[1]+1]))
                        v_.append(numpy.sum(width[:-1]**2*hdulist[1].data[insertloc[0]:insertloc[1]+1]))

                        # print myspec.dispersion[insertloc],ledge[fin:fin+2]
                        # print insertloc
                    fluxes.append(numpy.array(f_))
                    varis.append(numpy.array(v_))
                    phases.append(spec['salt2.phase'])

        phases=numpy.array(phases)
        if (len(phases) >=5 and (phases<0).sum() >=2 and (phases>0).sum() >=2):
            allfluxes.append(fluxes)
            allphases.append(phases)
            allvaris.append(varis)

output = open('data.pkl','wb')
pickle.dump((allphases, allfluxes, allvaris, ledge),output)
output.close()