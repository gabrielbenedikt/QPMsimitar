#!/usr/bin/python3

import numpy

class Filters:
    def __init__(self):
        self.FilterList = ['None','Rectangular','Gaussian']

    def getFilterFunction(self, type):
        if type=='None':
            return self.nofilter
        elif type=='Rectangular':
            return self.rectangularfilter
        elif type=='Gaussian':
            return self.gaussianfilter

    def nofilter(self, wl,cwl,fwhm):
        return wl

    def rectangularfilter(self,wl,cwl,fwhm):
        if (numpy.absolute(cwl-wl)>fwhm/2):
            return 0
        else:
            return 1

    def gaussianfilter(self,wl,cwl,fwhm):
        s=fwhm/2*numpy.sqrt(2*numpy.log(2))
        return numpy.exp(-((wl-cwl)**2)/(2*(s**2)))
