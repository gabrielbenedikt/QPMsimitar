#!/usr/bin/python3

import numpy

class Filters:
    def __init__(self):
        self.FilterList = ['None','Rectangular','Gaussian']

    def getFilterFunction(self, ftype,cwl,fwhm):
        if ftype=='None':
            #return self.nofilteronlywl(1,1)
            return None
        elif ftype=='Rectangular':
            return self.rectangularfilteronlywl(cwl,fwhm)
        elif ftype=='Gaussian':
            return self.gaussianfilteronlywl(cwl,fwhm)


    def nofilter(self, wl,cwl,fwhm):
        return 1

    def nofilteronlywl(self, cwl, fwhm):
        def func(wl):
            return self.nofilter(wl, cwl, fwhm)
        return func

    def rectangularfilter(self, wl, cwl, fwhm):
        if (numpy.absolute(cwl-wl)>fwhm/2):
            return 0
        else:
            return 1

    def rectangularfilteronlywl(self, cwl, fwhm):
        def func(wl):
            return self.rectangularfilter(wl, cwl, fwhm)
        return func


    def gaussianfilter(self,wl,cwl,fwhm):
        s=fwhm/2*numpy.sqrt(2*numpy.log(2))
        return numpy.exp(-((wl-cwl)**2)/(2*(s**2)))

    def gaussianfilteronlywl(self, cwl, fwhm):
        def func(wl):
            return self.gaussianfilter(wl, cwl, fwhm)
        return func
