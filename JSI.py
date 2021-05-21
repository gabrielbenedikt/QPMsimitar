#!/usr/bin/python3

import concurrent.futures
import numpy
import scipy
import scipy.optimize
import scipy.interpolate
from Constants import Constants
#import matplotlib
#import matplotlib.pyplot as plt
# import array
# from matplotlib.offsetbox import AnchoredText
# import os


class JSI:
    def __init__(self):
        # Thermal expansion coefficients of KTP
        # (Emanueli 2003)
        self.TXCa = 6.7 * 10 ** (-6)
        self.TXCb = 11 * 10 ** (-9)
        self.TXrefT = 25

        # signal and idler wavelength filtering
        self.useFilter = False  # change filterfunction(ls,li) below
        self.filtertype = 'rectangular'  # type of filter used. TODO: Implement different filter functions
        self.useabs = False  # take absolute value of pea?

        # Pump and crystal properties
        self.pwl = 775  * 10 ** (-9)
        self.tau = 2    * 10 ** (-12)

        self.m   = -1
        self.PP  = 46.2 * 10 ** (-6)
        self.T   = 25
        self.L   = 30   * 10 ** (-3)

        # plot JSA or JSI (choose only 1!)
        self.calcJSA = False
        self.calcJSI = True
        
        self.usetaucf = False

    # thermal expansion factor
    def thermexpfactor(self, T):
        return (1 + self.TXCa * (T - self.TXrefT) + self.TXCb * (T - self.TXrefT) * (T - self.TXrefT))

    # phasematching conditions
    def econv(self, lp, ls, li):
        return 1 / lp - 1 / ls - 1 / li

    def pconv(self, lp, ls, li, T, PP):
        # wavelengths in nm, but refractive-index-functions take µm
        pp = self.ny(lp, T) / lp
        ps = self.ny(ls, T) / ls
        pi = self.nz(li, T) / li
        pc = self.m / PP
        return pp - ps - pi - pc

    # this wrapper returns a n arroy of functions(phasematching conditions)
    # that only takes signal- and idler-wavelength as arguments
    def epconvonlywl(self, lp, T, PP):
        def epconv(x):
            # x[0]:lambda_p
            # x[0]:lambda_s
            # x[1]:lambda_i
            return [self.econv(lp, x[0], x[1]), self.pconv(lp, x[0], x[1], T, PP)]

        return epconv

    # returns signal and idler wavelengths that satisfy phasematching
    # conditions for a given pumpwavelength, temperature and poling period
    def SIwls(self, x):
        # x[0]: lambda_pump
        # x[1]: Temperature
        # x[2]: Poling period
        return scipy.optimize.fsolve(self.epconvonlywl(x[0], x[1], x[2]), [2 * x[0], 2 * x[0]], xtol=1e-6)

    # returns difference between signal and idler wavelength
    # for a given pumpwavelength, temperature and poling period
    def wlgap(self, x):
        # x[0]: lambda_pump
        # x[1]: Temperature
        # x[2]: Poling period
        swl, iwl = self.SIwls([x[0], x[1], x[2]])
        return swl - iwl

    # returns a function that only depends on the poling period
    def wlgaponlyT(self, lp, PP):
        def wlgap2(T):
            swl, iwl = self.SIwls([lp, T, PP])
            return swl - iwl

        return wlgap2

    # returns a function that only depends on the poling period
    def wlgaponlyPP(self, lp, T):
        def wlgap2(PP):
            swl, iwl = self.SIwls([lp, T, PP * self.thermexpfactor(T)])
            return swl - iwl
        return wlgap2

    def deltak(self, lp, ls, li, T, PP):
        tmp = self.pconv(lp, ls, li, T, PP)
        return 2 * Constants().pi * tmp

    # calculates pumpwavelength out of signal and idler wavelength, using energy conservation
    def lambdap(self, ls, li):
        return 1 / (1 / ls + 1 / li)

    # phase matching amplitude
    def PMA(self, dk, cl):
        # note: numpy.sinc(x) evaluates to Sin(pi*x)/(pi*x)
        # this really produced some headache.
        return numpy.sinc(dk * cl / (2 * Constants().pi))  # *numpy.exp(numpy.complex(0,1)*dk*cl/2)

    def PMAgauss(self, dk, cl):
        return self.PMA(dk, cl)

    def PMAsech(self, dk, cl):
        return self.PMA(dk, cl)

    def PMAsinc(self, dk, cl):
        return self.PMA(dk, cl)

    # phase matching intensity
    def PMIsech(self, dk, cl):
        return numpy.square(self.PMA(dk, cl))

    def PMIgauss(self, dk, cl):
        return numpy.square(self.PMA(dk, cl))

    def PMIsinc(self, dk, cl):
        return numpy.square(self.PMA(dk, cl))
    ###################
    ###   gaussian  ###
    ###################
    # pump envelope (gauss)
    # amplitude
    def PEAgauss(self, lp, ls, li, sp):
        dl = 1 / ls + 1 / li - 1 / lp
        return (numpy.exp(- (Constants().pi * Constants().c * (dl) / (sp)) ** 2))  # note: 2*pi*c/(2*sp)

    # intensity
    def PEIgauss(self, lp, ls, li, sp):
        return self.PEAgauss(lp, ls, li, sp) ** 2

    # joint spectral amplitude for gaussian beam
    def JSAgauss(self, lp, ls, li, tauac, t, pp, cl):
        # tauac: autocorrelator measured pulsewidth
        if self.usetaucf:
            tau = tauac * Constants().taucfgauss
        else:
            tau = tauac
        # tau=tauac
        dnu = Constants().tbwpgauss / tau  # FWHM in frequency
        dw = 2 * Constants().pi * dnu  # FWHM in angular frequency
        sp = dw / (2 * numpy.sqrt(2 * numpy.log(2)))  # gaussian standard deviation from FWHM
        dk = self.deltak(self.lambdap(ls, li), ls, li, t, pp)
        jsa = self.PEAgauss(lp, ls, li, sp) * self.PMAgauss(dk, cl)
        if self.useabs:
            return numpy.absolute(jsa)
        else:
            return jsa

    # joint spectral intensity for gaussian beam
    def JSIgauss(self, lp, ls, li, tauac, t, pp, cl):
        return (self.JSAgauss(lp, ls, li, tauac, t, pp, cl)) ** 2

    def PEAnPMAnJSAgauss(self, lp, ls, li, tauac, t, pp, cl):
        # tauac: autocorrelator measured pulsewidth
        if self.usetaucf:
            tau = tauac * Constants().taucfgauss
        else:
            tau = tauac
        # tau=tauac
        dnu = Constants().tbwpgauss / tau  # FWHM in frequency
        dw = 2 * Constants().pi * dnu  # FWHM in angular frequency
        sp = dw / (2 * numpy.sqrt(2 * numpy.log(2)))  # gaussian standard deviation from FWHM
        dk = self.deltak(self.lambdap(ls, li), ls, li, t, pp)
        pea = self.PEAgauss(lp, ls, li, sp)
        pma = self.PMAgauss(dk, cl)
        jsa = pea * pma
        if self.useabs:
            return [numpy.absolute(pea), numpy.absolute(pma), numpy.absolute(jsa)]
        else:
            return [pea, pma, jsa]

    def PEInPMInJSIgauss(self, lp, ls, li, tauac, t, pp, cl):
        [pei, pea, jsa] = self.PEAnPMAnJSAgauss(lp, ls, li, tauac, t, pp, cl)
        return [pei ** 2, pea ** 2, jsa ** 2]

    ###################
    ###   sech^2    ###
    ###################
    # pump envelope amplitude for sech^2 beam
    def PEAsech(self, lp, ls, li, B):
        wfact = 2 * Constants().pi * Constants().c * (1 / ls + 1 / li - 1 / lp)
        argument = wfact * B
        return (1 / numpy.cosh(argument))

    # pump envelope intensity for gaussian beam
    def PEIsech(self, lp, ls, li, B):
        return (self.PEAsech(lp, ls, li, B)) ** 2

    # joint spectral amplitude for sech^2 beam
    def JSAsech(self, lp, ls, li, tauac, t, pp, cl):
        if self.usetaucf:
            tau = tauac * Constants().taucfsech
        else:
            tau = tauac
        # tau=tauac
        dw = 2 * Constants().pi * Constants().tbwpsech / (tau)
        B = 2 * numpy.arccosh(numpy.sqrt(2)) / dw
        dk = self.deltak(self.lambdap(ls, li), ls, li, t, pp)
        jsa = self.PEAsech(lp, ls, li, B) * self.PMAsech(dk, cl)
        if self.useabs:
            return numpy.absolute(jsa)
        else:
            return jsa

    # joint spectral intensity for sech^2 beam
    def JSIsech(self, lp, ls, li, tauac, t, pp, cl):
        return self.JSAsech(lp, ls, li, tauac, t, pp, cl) ** 2

    def PEAnPMAnJSAsech(self, lp, ls, li, tauac, t, pp, cl):
        if self.usetaucf:
            tau = tauac * Constants().taucfsech
        else:
            tau = tauac
        # tau=tauac
        dw = 2 * Constants().pi * Constants().tbwpsech / (tau)
        B = 2 * numpy.arccosh(numpy.sqrt(2)) / dw
        dk = self.deltak(self.lambdap(ls, li), ls, li, t, pp)
        pea = self.PEAsech(lp, ls, li, B)
        pma = self.PMAsech(dk, cl)
        jsa = pea * pma
        if self.useabs:
            return [numpy.absolute(pea), numpy.absolute(pma), numpy.absolute(jsa)]
        else:
            return [pea, pma, jsa]

    def PEInPMInJSIsech(self, lp, ls, li, tauac, t, pp, cl):
        [pea, pma, jsa] = self.PEAnPMAnJSAsech(lp, ls, li, tauac, t, pp, cl)
        return [pea ** 2, pma ** 2, jsa ** 2]
    
    ###################
    ###    sinc     ###
    ###################
    # pump envelope amplitude for sinc beam
    def PEAsinc(self, lp, ls, li, B):
        dl = 1 / ls + 1 / li - 1 / lp
        wfact = 2 * Constants().pi * Constants().c * (1 / ls + 1 / li - 1 / lp)
        #argument = wfact * B
        argument = wfact*B
        return (numpy.sinc(argument))
        
    # pump envelope intensity for sinc beam
    def PEIsinc(self, lp, ls, li, B):
        return (self.PEAsinc(lp, ls, li, B)) ** 2
        
    # joint spectral amplitude for sinc beam
    def JSAsinc(self, lp, ls, li, tauac, t, pp, cl):
        #if self.usetaucf:
            #tau = tauac * Constants().taucfsech
        #else:
            #tau = tauac
        tau=tauac
        dw = 2 * Constants().pi * Constants().tbwpsinc / (tau)
        B = 3.79099 / dw
        dk = self.deltak(self.lambdap(ls, li), ls, li, t, pp)
        jsa = self.PEAsinc(lp, ls, li, B) * self.PMAsinc(dk, cl)
        if self.useabs:
            return numpy.absolute(jsa)
        else:
            return jsa

    # joint spectral intensity for sinc beam
    def JSIsinc(self, lp, ls, li, tauac, t, pp, cl):
        return self.JSAsinc(lp, ls, li, tauac, t, pp, cl) ** 2
        
    def PEAnPMAnJSAsinc(self, lp, ls, li, tauac, t, pp, cl):
        #if self.usetaucf:
            #tau = tauac * Constants().taucfsech
        #else:
            #tau = tauac
        tau=tauac
        dw = 2 * Constants().pi * Constants().tbwpsinc / (tau)
        B = 3.79099 / dw
        dk = self.deltak(self.lambdap(ls, li), ls, li, t, pp)
        pea = self.PEAsinc(lp, ls, li, B)
        pma = self.PMAsinc(dk, cl)
        jsa = pea * pma
        if self.useabs:
            return [numpy.absolute(pea), numpy.absolute(pma), numpy.absolute(jsa)]
        else:
            return [pea, pma, jsa]

    def PEInPMInJSIsinc(self, lp, ls, li, tauac, t, pp, cl):
        [pea, pma, jsa] = self.PEAnPMAnJSAsinc(lp, ls, li, tauac, t, pp, cl)
        return [pea ** 2, pma ** 2, jsa ** 2]


    def filterfunction(self, ls, li):
        # rectangular filter with:
        # lc......central wavelength
        # bwh..half bandwidth (such that transmittivity=1 in the range (lc-bwh,lc+bwh) ).
        lc = 1546.2 * 10 ** -(9)
        bwh = 1.5 * 10 ** (-9)

        # rectangular filter on both signal and idler
        if 0:
            if (ls < (lc - bwh) or ls > (lc + bwh) or li < (lc - bwh) or li > (lc + bwh)):
                return (0)
            else:
                return (1)
        # rectangular filter only on signal
        if 1:
            if (ls < (lc - bwh) or ls > (lc + bwh)):
                return (0)
            else:
                return (1)
        # rectangular filter only on idler
        if 0:
            if (li < (lc - bwh) or li > (lc + bwh)):
                return (0)
            else:
                return (1)

    def GetEffectivePP(self, m, Tcp, PPguess, lp, refidxfunc):
        self.nx = refidxfunc[0]
        self.ny = refidxfunc[1]
        self.nz = refidxfunc[2]
        return scipy.optimize.newton_krylov(self.wlgaponlyPP(lp, Tcp), PPguess, f_tol=1e-14)

    def getplots(self,pumpwl,signalrange,idlerrange,tau,temp,polingp,crystallength, refidxfunc,qpmorder,filter,plotJSI,pumpshape):
        print('start calculating JSA or JSI')
        #
        # pumpwl: Pump wavelength
        # signalrange: [double,double]: Signal wavelength range
        # idlerrange: [double,double]: Idler wavelength range
        # tau: pump pulse duration
        # temp: Temperature
        # polingp: Crystal poling period
        # crystallength: Length of crystal
        # refidxfunc: [nx,ny,nz]: Functions for refractive indices of crystal
        # qpmorder: Quasi phase matching order
        # filter: [function,function]: [filterfunction for signal, filterfunction for idler]
        # plotJSI: bool: True for JSI, false for JSA
        # pumpshape: string: Shape of pump beam (gaussian, sech^2)
        #

        self.sigrange = signalrange
        self.idrange = idlerrange
        self.nx = refidxfunc[0]
        self.ny = refidxfunc[1]
        self.nz = refidxfunc[2]

        self.pwl = pumpwl
        self.tau = tau
        self.m = qpmorder
        self.T = temp
        self.PP = polingp * self.thermexpfactor(self.T)
        self.L = crystallength * self.thermexpfactor(self.T)

        self.pumpshape = pumpshape

        self.filtermatrix = []

        self.useFilter = True
        self.filtersignalfunction = filter[0]
        self.filteridlerfunction = filter[1]
        for i in range(0, len(self.sigrange)):
            filtervector = []
            for j in range(0, len(self.idrange)):
                filterval=self.filteridlerfunction(self.idrange[j])*self.filtersignalfunction(self.sigrange[i])
                filtervector.append(filterval)
            self.filtermatrix.append(filtervector)

        if plotJSI==False:
            self.calcJSA = True
            self.calcJSI = False
        else:
            self.calcJSA = False
            self.calcJSI = True
            
        self.calcGaussian = False
        self.calcSech = False
        self.calcSinc = False
        if self.pumpshape.casefold() =='gaussian':
            self.calcGaussian = True
        elif self.pumpshape.casefold() == 'sinc':
            self.calcSinc = True
        else:
            self.calcSech = True

        X, Y = numpy.meshgrid(self.sigrange, self.idrange)

        if self.calcJSA:
            if self.calcGaussian:
                [PE, PM, JS] = self.PEAnPMAnJSAgauss(self.pwl, X, Y, self.tau, self.T, self.PP, self.L)
            elif self.calcSech:
                [PE, PM, JS] = self.PEAnPMAnJSAsech(self.pwl, X, Y, self.tau, self.T, self.PP, self.L)
            elif self.calcSinc:
                [PE, PM, JS] = self.PEAnPMAnJSAsinc(self.pwl, X, Y, self.tau, self.T, self.PP, self.L)
        elif self.calcJSI:
            if self.calcGaussian:
                [PE, PM, JS] = self.PEInPMInJSIgauss(self.pwl, X, Y, self.tau, self.T, self.PP, self.L)
            elif self.calcSech:
                [PE, PM, JS] = self.PEInPMInJSIsech(self.pwl, X, Y, self.tau, self.T, self.PP, self.L)
            elif self.calcSinc:
                [PE, PM, JS] = self.PEInPMInJSIsinc(self.pwl, X, Y, self.tau, self.T, self.PP, self.L)
        else:
            print('Error: Calc neither JSA nor JSI.')
            return
        
        if self.useFilter:
            JS = JS * self.filtermatrix

        return [PE, PM, JS]

    def getpurity_vsTau(self,pumpwl,signalrange,idlerrange,taurange,temp,polingp,crystallength,refidxfunc,qpmorder,filter,pumpshape):
        #
        # pumpwl: Pump wavelength
        # signalrange: [double,double]: Signal wavelength range
        # idlerrange: [double,double]: Idler wavelength range
        # taurange: range of pump pulse durations
        # temp: Temperature
        # polingp: Crystal poling period
        # crystallength: Length of crystal
        # refidxfunc: [nx,ny,nz]: Functions for refractive indices of crystal
        # qpmorder: Quasi phase matching order
        # filter: [string,bool,bool]: [Type of filter to use, True: use filter for signal, True: use filter for idler]
        # pumpshape: string: Shape of pump beam (gaussian, sech^2)
        #

        print('start calculating JSA or JSI')

        self.sigrange = signalrange
        self.idrange = idlerrange
        self.nx = refidxfunc[0]
        self.ny = refidxfunc[1]
        self.nz = refidxfunc[2]

        self.pwl = pumpwl
        self.taurange = taurange
        self.m = qpmorder
        self.T = temp
        self.PP = polingp * self.thermexpfactor(self.T)
        self.L = crystallength * self.thermexpfactor(self.T)

        self.pumpshape = pumpshape

        self.filtermatrix = []

        self.useFilter = True
        self.filtersignalfunction = filter[0]
        self.filteridlerfunction = filter[1]
        for i in range(0, len(self.sigrange)):
            filtervector = []
            for j in range(0, len(self.idrange)):
                filterval=self.filteridlerfunction(self.sigrange[i])*self.filtersignalfunction(self.idrange[j])
                filtervector.append(filterval)
            self.filtermatrix.append(filtervector)
            
        self.calcGaussian = False
        self.calcSech = False
        self.calcSinc = False
        if self.pumpshape.casefold() =='gaussian':
            self.calcGaussian = True
        elif self.pumpshape.casefold() == 'sinc':
            self.calcSinc = True
        else:
            self.calcSech = True

        X, Y = numpy.meshgrid(self.sigrange, self.idrange)

        purity = []

        for i in range(0, len(self.taurange)):
            self.tau = self.taurange[i]
            # JSA
            if self.calcGaussian:
                JSA = self.JSAgauss(self.pwl, X, Y, self.tau, self.T, self.PP, self.L)
            elif self.calcSech:
                JSA = self.JSAsech(self.pwl, X, Y, self.tau, self.T, self.PP, self.L)
            elif self.calcSinc:
                JSA = self.JSAsinc(self.pwl, X, Y, self.tau, self.T, self.PP, self.L)
                
            if self.useFilter:
                JSA = JSA * self.filtermatrix

            # Purity
            # SDV
            sA = scipy.linalg.svd(JSA, overwrite_a=True, compute_uv=False)
            # singular values in s. need to normalize s to get the schmidt magnitudes
            snA = sA / scipy.linalg.norm(sA, 2)
            # sum of squares of schmidt magnitudes is purity
            purity.append(numpy.sum(snA ** 4))
            #infostring = "Pulsewidth: {0:.5f}ps\tpurity: {1:.5f}".format(self.tau * 10 ** (12), purity[i])
            #print(infostring)

        # interpolate purity curve
        increased_taurange = numpy.linspace(self.taurange[0], self.taurange[-1], 100000)
        interPur = scipy.interpolate.InterpolatedUnivariateSpline(self.taurange, purity)
        interPurvals = interPur(increased_taurange)
        max = numpy.max(interPurvals)
        maxidx = numpy.argmax(interPurvals)

        #print("maximum:\n\t\t\tpurity: {0:.4f}".format(max))

        return [purity,max,increased_taurange[maxidx]]

    def getpurity_vsL(self,pumpwl,signalrange,idlerrange,tau,temp,polingp,crystallengthrange,refidxfunc,qpmorder,filter,pumpshape):
        #
        # pumpwl: Pump wavelength
        # signalrange: [double,double]: Signal wavelength range
        # idlerrange: [double,double]: Idler wavelength range
        # taurange: range of pump pulse durations
        # temp: Temperature
        # polingp: Crystal poling period
        # crystallength: Length of crystal
        # refidxfunc: [nx,ny,nz]: Functions for refractive indices of crystal
        # qpmorder: Quasi phase matching order
        # filter: [string,bool,bool]: [Type of filter to use, True: use filter for signal, True: use filter for idler]
        # pumpshape: string: Shape of pump beam (gaussian, sech^2)
        #

        #print('start calculating JSA or JSI')

        self.sigrange = signalrange
        self.idrange = idlerrange
        self.nx = refidxfunc[0]
        self.ny = refidxfunc[1]
        self.nz = refidxfunc[2]

        self.pwl = pumpwl
        self.tau = tau
        self.m = qpmorder
        self.T = temp
        self.PP = polingp * self.thermexpfactor(self.T)
        self.Lrange = crystallengthrange * self.thermexpfactor(self.T)

        self.pumpshape = pumpshape

        self.filtermatrix = []

        self.useFilter = True
        self.filtersignalfunction = filter[0]
        self.filteridlerfunction = filter[1]
        for i in range(0, len(self.sigrange)):
            filtervector = []
            for j in range(0, len(self.idrange)):
                filterval=self.filteridlerfunction(self.sigrange[i])*self.filtersignalfunction(self.idrange[j])
                filtervector.append(filterval)
            self.filtermatrix.append(filtervector)

        self.calcGaussian = False
        self.calcSech = False
        self.calcSinc = False
        if self.pumpshape.casefold() == 'gaussian':
            self.calcGaussian = True
        elif self.pumpshape.casefold() == 'sech^2':
            self.calcSech = True
        elif self.pumpshape.casefold() == 'sinc':
            self.calcSinc = True

        X, Y = numpy.meshgrid(self.sigrange, self.idrange)

        purity = []

        for i in range(0, len(self.Lrange)):
            self.L = self.Lrange[i]
            # JSA
            if self.calcGaussian:
                JSA = self.JSAgauss(self.pwl, X, Y, self.tau, self.T, self.PP, self.L)
            elif self.calcSech:
                JSA = self.JSAsech(self.pwl, X, Y, self.tau, self.T, self.PP, self.L)
            elif self.calcSinc:
                JSA = self.JSAsinc(self.pwl, X, Y, self.tau, self.T, self.PP, self.L)
                
            if self.useFilter:
                JSA = JSA * self.filtermatrix

            # Purity
            # SDV
            sA = scipy.linalg.svd(JSA, overwrite_a=True, compute_uv=False)
            # singular values in s. need to normalize s to get the schmidt magnitudes
            snA = sA / scipy.linalg.norm(sA, 2)
            # sum of squares of schmidt magnitudes is purity
            purity.append(numpy.sum(snA ** 4))
            #infostring = "Crystal length: {0:.5f}ps\tpurity: {1:.5f}".format(self.L * 10 ** (3), purity[i])
            #print(infostring)

        # interpolate purity curve
        increased_Lrange = numpy.linspace(self.Lrange[0], self.Lrange[-1], 100000)
        interPur = scipy.interpolate.InterpolatedUnivariateSpline(self.Lrange, purity)
        interPurvals = interPur(increased_Lrange)
        max = numpy.max(interPurvals)
        maxidx = numpy.argmax(interPurvals)

        #print("maximum:\n\t\t\tpurity: {0:.4f}".format(max))

        return [purity,max,increased_Lrange[maxidx]]

    def getTcpVslp(self,pwlrange,temp,polingp,refidxfunc,qpmorder):
        self.PP=polingp
        self.T = temp
        [self.nx,self.ny,self.nz] = refidxfunc
        self.m=qpmorder

        plotrange = pwlrange
        Tcp = []
        initialTguess = self.T
        Tguess = initialTguess
        for j in range(0, len(plotrange)):
            if j:
                Tguess = Tcp[j - 1]
            tmp = scipy.optimize.fsolve(self.wlgaponlyT(plotrange[j], self.PP), Tguess)
            Tcp.append(tmp)

        return Tcp

    def getTcpVsPP(self,PPrange,temp,pwl,refidxfunc,qpmorder):
        self.pwl=pwl
        self.T = temp
        [self.nx,self.ny,self.nz] = refidxfunc
        self.m=qpmorder

        plotrange = PPrange
        Tcp = []
        initialTguess = self.T
        Tguess = initialTguess
        for j in range(0, len(plotrange)):
            if j:
                Tguess = Tcp[j - 1]
            tmp = scipy.optimize.fsolve(self.wlgaponlyT(pwl, plotrange[j]), Tguess)
            Tcp.append(tmp)

        return Tcp

    def getHOMinterference(self, pwl, temp, polingp, qpmorder, tau, cl, signalrange, idlerrange,JSIresolution, pumpshape, delayrange, refidxfunc, filter):
        [self.nx, self.ny, self.nz] = refidxfunc
        #
        # https://arxiv.org/pdf/1211.0120.pdf (On the Purity and Indistinguishability of Down-Converted Photons. Osorio, Sangouard, thew 2012)
        # Ansari, 2013 msc thesis
        #
        self.calcGaussian = False
        self.calcSech = False
        self.calcSinc = False
        if pumpshape.casefold() =='gaussian':
            self.calcGaussian = True
        elif pumpshape.casefold() =='sech^2':
            self.calcSech = True
        elif pumpshape.casefold() =='sinc':
            self.calcSinc = True

        self.filtermatrix = []

        self.useFilter = True
        self.filtersignalfunction = filter[0]
        self.filteridlerfunction = filter[1]
        for i in range(0, len(signalrange)):
            filtervector = []
            for j in range(0, len(idlerrange)):
                filterval = self.filteridlerfunction(signalrange[i]) * self.filtersignalfunction(idlerrange[j])
                filtervector.append(filterval)
            self.filtermatrix.append(filtervector)

        X, Y = numpy.meshgrid(signalrange, idlerrange)

        HOMI = []

        if self.calcSech:
            jsi = self.JSIsech(pwl, X, Y, tau, temp, polingp, cl)
            jsa = self.JSAsech(pwl, X, Y, tau, temp, polingp, cl)
            jsa_cc = numpy.conjugate(self.JSAsech(pwl, Y, X, tau, temp, polingp, cl))
        elif self.calcGaussian:
            jsi = self.JSIgauss(pwl, X, Y, tau, temp, polingp, cl)
            jsa = self.JSAgauss(pwl, X, Y, tau, temp, polingp, cl)
            jsa_cc = numpy.conjugate(self.JSAgauss(pwl, Y, X, tau, temp, polingp, cl))
        elif self.calcSinc:
            jsi = self.JSIsinc(pwl, X, Y, tau, temp, polingp, cl)
            jsa = self.JSAsinc(pwl, X, Y, tau, temp, polingp, cl)
            jsa_cc = numpy.conjugate(self.JSAsinc(pwl, Y, X, tau, temp, polingp, cl))
        else:
            print('ERROR: Unknown pump beamshape')
        if self.useFilter:
            jsa = jsa*self.filtermatrix
            jsi = jsi*self.filtermatrix
            jsa_cc = jsa_cc*self.filtermatrix
        
        if 0:
            for i in range(0,len(delayrange)):
                phase = numpy.exp(1j*2*numpy.pi*Constants().c*(1/X-1/Y)*delayrange[i])
                ProbMX = jsi-phase*jsa*jsa_cc
                
                #tmp=0
                #for i in range(0,len(X)):
                    #for j in range(0,len(Y)):
                        #tmp = tmp + ProbMX[i][j]
                HOMI.append(numpy.sum(ProbMX))
        else:
            def homf(i):
                    phase = numpy.exp(1j*2*numpy.pi*Constants().c*(1/X-1/Y)*delayrange[i])
                    ProbMX = jsi-phase*jsa*jsa_cc
                    return numpy.sum(ProbMX)
            with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
                futures = [ex.submit(homf, i) for i in range(0,len(delayrange))]
                HOMI = [f.result() for f in futures]

        # omit tiny imaginary parts
        HOMI = numpy.abs(HOMI)

        # norm to 1/2
        max = numpy.max(HOMI)
        HOMI = HOMI/max
        HOMI = 0.5*HOMI

        # determine visibility
        homimax=numpy.max(HOMI)
        homimin=numpy.min(HOMI)
        vis = numpy.abs((homimax-homimin)/(homimax))

        # calc FWHM
        #visinterpolf = scipy.interpolate.interp1d(delayrange, HOMI-0.25)
        #negrootstartest=delayrange[int(numpy.floor(len(delayrange)/3))]
        #posrootstartest=delayrange[int(numpy.floor(len(delayrange)*2/3))]
        #print('Estimate for negative root: ', negrootstartest)
        #print('Estimate for positive root: ', posrootstartest)
        #negroot = scipy.optimize.fsolve(visinterpolf, negrootstartest)
        #posroot = scipy.optimize.fsolve(visinterpolf, posrootstartest)
        
        #anorther try, maybe more stable=
        delayrangeneg = delayrange[:int(numpy.floor(len(delayrange)/2))]
        delayrangepos = delayrange[int(numpy.floor(len(delayrange)/2)):]
        HOMIneg = HOMI[:int(numpy.floor(len(HOMI)/2))]
        HOMIpos = HOMI[int(numpy.floor(len(HOMI)/2)):]
        visinterpolfneg = scipy.interpolate.interp1d(delayrangeneg, HOMIneg-0.25, fill_value='extrapolate')
        visinterpolfpos = scipy.interpolate.interp1d(delayrangepos, HOMIpos-0.25, fill_value='extrapolate')
        negrootstartest=delayrangeneg[int(numpy.floor(len(delayrangeneg)/2))]
        posrootstartest=delayrangepos[int(numpy.floor(len(delayrangepos)/2))]
        negroot=scipy.optimize.fsolve(visinterpolfneg, negrootstartest)
        posroot=scipy.optimize.fsolve(visinterpolfpos, posrootstartest)

        homfwhm=posroot[0]-negroot[0]

        return [HOMI,vis,homfwhm]

    def getFWHMvstau(self, pwl, signalrange, idlerrange, temp, polingp, qpmorder, cl, taurange, refidxfunc, filter, JSIresolution, pumpshape, decprec, usetaucf):
        [self.nx, self.ny, self.nz] = refidxfunc
        self.usetaucf = usetaucf
        
        self.calcGaussian = False
        self.calcSech = False
        self.calcSinc = False
        if pumpshape.casefold() =='gaussian':
            self.calcGaussian = True
        elif pumpshape.casefold() =='sech^2':
            self.calcSech = True
        elif pumpshape.casefold() =='sinc':
            self.calcSinc = True

        self.filtermatrix = []

        self.useFilter = True
        self.filtersignalfunction = filter[0]
        self.filteridlerfunction = filter[1]
        for i in range(0, len(signalrange)):
            filtervector = []
            for j in range(0, len(idlerrange)):
                filterval = self.filteridlerfunction(signalrange[i]) * self.filtersignalfunction(idlerrange[j])
                filtervector.append(filterval)
            self.filtermatrix.append(filtervector)
        
        fwhmsig=[]
        fwhmid=[]
        
        #TMP
        self.useFilter=False
        
        for h in range(0,len(taurange)):
            if self.calcGaussian:
                result =self.JSIgauss(pwl, signalrange[:,None], idlerrange[None,:], taurange[h], temp, polingp*self.thermexpfactor(temp), cl)-0.5
            elif self.calcSech: 
                result = self.JSIsech(pwl, signalrange[:,None], idlerrange[None,:], taurange[h], temp, polingp*self.thermexpfactor(temp), cl)-0.5
            elif self.calcSinc: 
                result = self.JSIsinc(pwl, signalrange[:,None], idlerrange[None,:], taurange[h], temp, polingp*self.thermexpfactor(temp), cl)-0.5
            else:
                print('ERROR: Unknown pump beamshape')
            if self.useFilter:
                result = result*self.filtermatrix
            hmptss=[]
            hmptsi=[]
            hmpts=[]
            
            #option 1
            result=result.round(decprec)
            idcs=numpy.where(result==0) 
            for i in range(0,len(idcs[0])):
                    hmptss.append(signalrange[idcs[0][i]]*10**9)
                    hmptsi.append(idlerrange[idcs[1][i]]*10**9)
                    
            #option 2
            #for i in range(0,len(signalrange)):
            #    for j in range(0,len(idlerrange)):
            #        if (numpy.abs(result[i][j])<10**(-decprec)):
            #            hmptss.append(signalrange[i]*10**9)
            #            hmptsi.append(idlerrange[j]*10**9)
            
            if len(hmptss) != 0:
                if len(hmptsi) !=0:
                    fwhmsig.append(max(hmptss)-min(hmptss))
                    fwhmid.append(max(hmptsi)-min(hmptsi))
            
        return [fwhmsig,fwhmid]
