#!/usr/bin/env python3

import cProfile

import concurrent.futures
import multiprocessing as mp
import numpy as np
import scipy
import scipy.optimize
import scipy.interpolate
import scipy.integrate
from datetime import datetime
from Constants import Constants

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
            return [self.econv(lp, x[0], x[1]), float(self.pconv(lp, x[0], x[1], T, PP))]

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
        return 2 * Constants().pi * self.pconv(lp, ls, li, T, PP)

    # calculates pumpwavelength out of signal and idler wavelength, using energy conservation
    def lambdap(self, ls, li):
        return 1 / (1 / ls + 1 / li)

    # phase matching amplitude
    def PMA(self, dk, cl):
        # note: np.sinc(x) evaluates to Sin(pi*x)/(pi*x)
        # this really produced some headache.
        
        arg = cl * dk / 2
        return cl * np.sinc(arg/np.pi ) * np.exp( 1j * arg )
        
        #return np.sinc(dk * cl / (2 * Constants().pi)) #* cl * np.exp(1j*dk*cl/2)

    def PMAgauss(self, dk, cl):
        return self.PMA(dk, cl)

    def PMAsech(self, dk, cl):
        return self.PMA(dk, cl)

    def PMAsinc(self, dk, cl):
        return self.PMA(dk, cl)

    # phase matching intensity
    def PMIsech(self, dk, cl):
        return np.square(self.PMA(dk, cl))

    def PMIgauss(self, dk, cl):
        return np.square(self.PMA(dk, cl))

    def PMIsinc(self, dk, cl):
        return np.square(self.PMA(dk, cl))
    ###################
    ###   gaussian  ###
    ###################
    # pump envelope (gauss)
    # amplitude
    def PEAgauss(self, lp, ls, li, sp):
        dl = 1 / ls + 1 / li - 1 / lp
        return (np.exp(- (Constants().pi * Constants().c * (dl) / (sp)) ** 2))  # note: 2*pi*c/(2*sp)

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
        sp = dw / (2 * np.sqrt(2 * np.log(2)))  # gaussian standard deviation from FWHM
        dk = self.deltak(self.lambdap(ls, li), ls, li, t, pp)
        jsa = self.PEAgauss(lp, ls, li, sp) * self.PMAgauss(dk, cl)
        if self.useabs:
            return np.absolute(jsa)
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
        sp = dw / (2 * np.sqrt(2 * np.log(2)))  # gaussian standard deviation from FWHM
        dk = self.deltak(self.lambdap(ls, li), ls, li, t, pp)
        pea = self.PEAgauss(lp, ls, li, sp)
        pma = self.PMAgauss(dk, cl)
        jsa = pea * pma
        if self.useabs:
            return [np.absolute(pea), np.absolute(pma), np.absolute(jsa)]
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
        return (1 / np.cosh(argument))

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
        B = 2 * np.arccosh(np.sqrt(2)) / dw
        
        dk = self.deltak(self.lambdap(ls, li), ls, li, t, pp)
        jsa = self.PEAsech(lp, ls, li, B) * self.PMAsech(dk, cl)
        
        if self.useabs:
            return np.absolute(jsa)
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
        B = 2 * np.arccosh(np.sqrt(2)) / dw
        dk = self.deltak(self.lambdap(ls, li), ls, li, t, pp)
        pea = self.PEAsech(lp, ls, li, B)
        pma = self.PMAsech(dk, cl)
        jsa = pea * pma
        if self.useabs:
            return [np.absolute(pea), np.absolute(pma), np.absolute(jsa)]
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
        return (np.sinc(argument))
        
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
            return np.absolute(jsa)
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
            return [np.absolute(pea), np.absolute(pma), np.absolute(jsa)]
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

    def getplots(self,pumpwl,signalrange,idlerrange,tau,temp,polingp,crystallength, refidxfunc,qpmorder,filterfuncs,plotJSI,pumpshape):
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
        # filterfuncs: [function,function]: [filterfunction for signal, filterfunction for idler]
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

        self.calcJSI = False
        self.calcJSA = False
        if plotJSI==False:
            self.calcJSA = True
        else:
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

        X, Y = np.meshgrid(self.sigrange, self.idrange)

        if self.calcJSA:
            if self.calcGaussian:
                [PE, PM, JS] = self.PEAnPMAnJSAgauss(self.pwl, X, Y, self.tau, self.T, self.PP, self.L)
            elif self.calcSech:
                [PE, PM, JS] = self.PEAnPMAnJSAsech(self.pwl, X, Y, self.tau, self.T, self.PP, self.L)
            elif self.calcSinc:
                [PE, PM, JS] = self.PEAnPMAnJSAsinc(self.pwl, X, Y, self.tau, self.T, self.PP, self.L)
        elif self.calcJSI:
            if self.calcGaussian:
                #arguments=[(self.pwl, swl[0], iwl[0], self.tau, self.T, self.PP, self.L) for swl, iwl in zip(self.sigrange, self.idrange)]
                #print(arguments[0])
                #with mp.Pool(processes=8) as pool:
                    #PE, PM, JS = pool.starmap(self.PEInPMInJSIgauss, arguments)
                ##PE = np.array(PE)
                ##PM = np.array(PM)
                ##JS = np.array(JS)
                [PE, PM, JS] = self.PEInPMInJSIgauss(self.pwl, X, Y, self.tau, self.T, self.PP, self.L)
            elif self.calcSech:
                [PE, PM, JS] = self.PEInPMInJSIsech(self.pwl, X, Y, self.tau, self.T, self.PP, self.L)
            elif self.calcSinc:
                [PE, PM, JS] = self.PEInPMInJSIsinc(self.pwl, X, Y, self.tau, self.T, self.PP, self.L)
        else:
            print('Error: Calc neither JSA nor JSI.')
            return
        if filterfuncs[0]==None and filterfuncs[1]==None:
            self.useFilter = False
        else:
            self.useFilter = True
        print('self.usefilter: ', self.useFilter)
        if self.useFilter:
            #self.useFilter = True
            self.filtermatrix = []
            self.filtersignalfunction = filterfuncs[0]
            self.filteridlerfunction = filterfuncs[1]
            for i in range(0, len(self.sigrange)):
                filtervector = []
                for j in range(0, len(self.idrange)):
                    filterval=self.filteridlerfunction(self.idrange[j])*self.filtersignalfunction(self.sigrange[i])
                    filtervector.append(filterval)
                self.filtermatrix.append(filtervector)
            JS = JS * self.filtermatrix

        return [PE, PM, JS]

    def getpurity_vsTau(self,pumpwl,signalrange,idlerrange,taurange,temp,polingp,crystallength,refidxfunc,qpmorder,filterfuncs,pumpshape):
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
            
        self.calcGaussian = False
        self.calcSech = False
        self.calcSinc = False
        if self.pumpshape.casefold() =='gaussian':
            self.calcGaussian = True
        elif self.pumpshape.casefold() == 'sinc':
            self.calcSinc = True
        else:
            self.calcSech = True

        X, Y = np.meshgrid(self.sigrange, self.idrange)

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
                self.filtersignalfunction = filterfuncs[0]
                self.filteridlerfunction = filterfuncs[1]
                if not (self.filtersignalfunction==None and self.filteridlerfunction==None):
                    self.filtermatrix = []
                    self.useFilter = True
                    for i in range(0, len(self.sigrange)):
                        filtervector = []
                        for j in range(0, len(self.idrange)):
                            filterval=self.filteridlerfunction(self.sigrange[i])*self.filtersignalfunction(self.idrange[j])
                            filtervector.append(filterval)
                        self.filtermatrix.append(filtervector)
                    JSA = JSA * self.filtermatrix

            # Purity
            # SDV
            sA = scipy.linalg.svd(JSA, overwrite_a=True, compute_uv=False)
            # singular values in s. need to normalize s to get the schmidt magnitudes
            snA = sA / scipy.linalg.norm(sA, 2)
            # sum of squares of schmidt magnitudes is purity
            purity.append(np.sum(snA ** 4))
            #infostring = "Pulsewidth: {0:.5f}ps\tpurity: {1:.5f}".format(self.tau * 10 ** (12), purity[i])
            #print(infostring)

        # interpolate purity curve
        increased_taurange = np.linspace(self.taurange[0], self.taurange[-1], 100000)
        interPur = scipy.interpolate.InterpolatedUnivariateSpline(self.taurange, purity)
        interPurvals = interPur(increased_taurange)
        max = np.max(interPurvals)
        maxidx = np.argmax(interPurvals)

        #print("maximum:\n\t\t\tpurity: {0:.4f}".format(max))

        return [purity,max,increased_taurange[maxidx]]

    def getpurity_vsL(self,pumpwl,signalrange,idlerrange,tau,temp,polingp,crystallengthrange,refidxfunc,qpmorder,filterfuncs,pumpshape):
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
        # filterfuncs: [string,bool,bool]: [Type of filter to use, True: use filter for signal, True: use filter for idler]
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

        self.calcGaussian = False
        self.calcSech = False
        self.calcSinc = False
        if self.pumpshape.casefold() == 'gaussian':
            self.calcGaussian = True
        elif self.pumpshape.casefold() == 'sech^2':
            self.calcSech = True
        elif self.pumpshape.casefold() == 'sinc':
            self.calcSinc = True

        X, Y = np.meshgrid(self.sigrange, self.idrange)

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
                self.filtersignalfunction = filterfuncs[0]
                self.filteridlerfunction = filterfuncs[1]
                if not (self.filtersignalfunction==None and self.filteridlerfunction==None):
                    self.filtermatrix = []
                    self.useFilter = True
                    for i in range(0, len(self.sigrange)):
                        filtervector = []
                        for j in range(0, len(self.idrange)):
                            filterval=self.filteridlerfunction(self.sigrange[i])*self.filtersignalfunction(self.idrange[j])
                            filtervector.append(filterval)
                    self.filtermatrix.append(filtervector)

                    JSA = JSA * self.filtermatrix

            # Purity
            # SDV
            sA = scipy.linalg.svd(JSA, overwrite_a=True, compute_uv=False)
            # singular values in s. need to normalize s to get the schmidt magnitudes
            snA = sA / scipy.linalg.norm(sA, 2)
            # sum of squares of schmidt magnitudes is purity
            purity.append(np.sum(snA ** 4))
            #infostring = "Crystal length: {0:.5f}ps\tpurity: {1:.5f}".format(self.L * 10 ** (3), purity[i])
            #print(infostring)

        # interpolate purity curve
        increased_Lrange = np.linspace(self.Lrange[0], self.Lrange[-1], 100000)
        interPur = scipy.interpolate.InterpolatedUnivariateSpline(self.Lrange, purity)
        interPurvals = interPur(increased_Lrange)
        max = np.max(interPurvals)
        maxidx = np.argmax(interPurvals)

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

#     #
#     # "original" method that didn's work as expected
#     #
#     def getHOMinterference_orig(self, pwl, temp, polingp, qpmorder, tau, cl, signalrange, idlerrange,JSIresolution, pumpshape, delayrange, homphase, refidxfunc, filterfuncs):
#         t0=datetime.now()
#         [self.nx, self.ny, self.nz] = refidxfunc
#         #
#         # https://arxiv.org/pdf/1211.0120.pdf (On the Purity and Indistinguishability of Down-Converted Photons. Osorio, Sangouard, thew 2012)
#         # Ansari, 2013 msc thesis
#         #
#         self.calcGaussian = False
#         self.calcSech = False
#         self.calcSinc = False
#         if pumpshape.casefold() =='gaussian':
#             self.calcGaussian = True
#         elif pumpshape.casefold() =='sech^2':
#             self.calcSech = True
#         elif pumpshape.casefold() =='sinc':
#             self.calcSinc = True
# 
#         X, Y = np.meshgrid(signalrange, idlerrange)
# 
#         HOMI = []
# 
#         if self.calcSech:
#             jsa = self.JSAsech(pwl, X, Y, tau, temp, polingp, cl)
#             jsa_t = self.JSAsech(pwl, Y, X, tau, temp, polingp, cl)
#         elif self.calcGaussian:
#             jsa = self.JSAgauss(pwl, X, Y, tau, temp, polingp, cl)
#             jsa_t = self.JSAgauss(pwl, Y, X, tau, temp, polingp, cl)
#         elif self.calcSinc:
#             jsa = self.JSAsinc(pwl, X, Y, tau, temp, polingp, cl)
#             jsa_t = self.JSAsinc(pwl, Y, X, tau, temp, polingp, cl)
#         else:
#             print('ERROR: Unknown pump beamshape')
#         #jsa_t = jsa_t / np.sum(jsa_t)
#         jsa_tcc = np.conjugate(jsa_t)
#         # jsa should integrate to 1]
#         np.save('jsa.npy', jsa)
#         np.save('jsa_t.npy', jsa_t)
#         np.save('X.npy', X)
#         np.save('Y.npy', Y)
#         jsa = jsa / np.sum(jsa)
#         jsa_tcc = jsa_tcc / np.sum(jsa_tcc)
#         
#         #filters
#         #self.useFilter = True ## TODO: why was this here? Oo
#         self.filtermatrix = []
#         self.filtersignalfunction = filterfuncs[0]
#         self.filteridlerfunction = filterfuncs[1]
#         if not (self.filtersignalfunction==None and self.filteridlerfunction == None):
#             for i in range(0, len(signalrange)):
#                 filtervector = []
#                 for j in range(0, len(idlerrange)):
#                     filterval = self.filteridlerfunction(signalrange[i]) * self.filtersignalfunction(idlerrange[j])
#                     filtervector.append(filterval)
#                 self.filtermatrix.append(filtervector)
#             jsa = jsa*self.filtermatrix#
#             jsa_tcc = jsa_tcc*np.transpose(self.filtermatrix)
#         
#         
#         def homf(i):
#             #phase = np.exp(1j*(2*np.pi*Constants().c*(1/X-1/Y)*delayrange[i]+homphase))
#             phase = np.exp(1j*(2*np.pi*Constants().c*(1/X-1/Y)*delayrange[i]))
#             jj = jsa*jsa_tcc*phase
#             #ProbMX = 0.5 * (1 - np.sum(jsa_tcc*jsa*phase))
#             ProbMX = 0.5 * (1 - np.sum(jj))
#             return ProbMX
#         
#         with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
#             futures = [ex.submit(homf, i) for i in range(0,len(delayrange))]
#             HOMI = [f.result() for f in futures]
# 
#         # omit tiny imaginary parts
#         HOMI = np.abs(HOMI)
#         
#         # determine visibility
#         vis = np.abs((np.max(HOMI)-np.min(HOMI))/(np.max(HOMI)))
# 
#         # calc FWHM
#         delayrangeneg = delayrange[:int(np.floor(len(delayrange)/2))]
#         delayrangepos = delayrange[int(np.floor(len(delayrange)/2)):]
#         HOMIneg = HOMI[:int(np.floor(len(HOMI)/2))]
#         HOMIpos = HOMI[int(np.floor(len(HOMI)/2)):]
#         visinterpolfneg = scipy.interpolate.interp1d(delayrangeneg, HOMIneg-0.25, fill_value='extrapolate')
#         visinterpolfpos = scipy.interpolate.interp1d(delayrangepos, HOMIpos-0.25, fill_value='extrapolate')
#         negrootstartest=delayrangeneg[int(np.floor(len(delayrangeneg)/2))]
#         posrootstartest=delayrangepos[int(np.floor(len(delayrangepos)/2))]
#         negroot=scipy.optimize.fsolve(visinterpolfneg, negrootstartest)
#         posroot=scipy.optimize.fsolve(visinterpolfpos, posrootstartest)
# 
#         homfwhm=posroot[0]-negroot[0]
# 
#         t1=datetime.now()
#         print('calculating HOM took', (t1-t0).total_seconds(), 's')
#         return [HOMI,vis,homfwhm]
    
#     #adaptation of SPDCalc.org
#     def getHOMinterference_new2(self, pwl, temp, polingp, qpmorder, tau, cl, signalrange, idlerrange,JSIresolution, pumpshape, delayrange, homphase, refidxfunc, filterfuncs):
#         t0=datetime.now()
#         [self.nx, self.ny, self.nz] = refidxfunc
#         #
#         # https://arxiv.org/pdf/1211.0120.pdf (On the Purity and Indistinguishability of Down-Converted Photons. Osorio, Sangouard, thew 2012)
#         # Ansari, 2013 msc thesis
#         #
#         self.calcSech = False
#         self.calcSech = True
# 
#         X, Y = np.meshgrid(signalrange, idlerrange)
# 
#         HOMI = []
#         
#         jsa1 = self.JSAsech(pwl, X, Y, tau, temp, polingp, cl)
#         jsa2 = self.JSAsech(pwl, Y, X, tau, temp, polingp, cl)
#  
#         jsi = self.JSIsech(pwl, X, Y, tau, temp, polingp, cl)
#         norm = np.sum(jsi)
# 
#         jsa1_r = np.real(jsa1)
#         jsa1_i = np.imag(jsa1)
#         jsa2_r = np.real(jsa2[::-1,::-1].T)
#         jsa2_i = np.imag(jsa2[::-1,::-1].T)
#         
#         def homf(i):
#             ProbMX = 0;
#             phase = 2*np.pi*Constants().c *(1/signalrange - 1/idlerrange.T)*delayrange[i]
#             Tosc_real = np.cos(phase)
#             Tosc_imag = np.sin(phase)
#             arg2_real = Tosc_real*jsa2_r - Tosc_imag*jsa2_i
#             arg2_imag = Tosc_real*jsa2_i + Tosc_imag*jsa2_r
#             PM_real = (jsa1_r - arg2_real)
#             PM_imag = (jsa1_i - arg2_imag)
#             val =  PM_real**2 + PM_imag**2
#             ProbMX = 0.25 * np.sum(val)
#             return ProbMX
#         
#         with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
#             futures = [ex.submit(homf, i) for i in range(0,len(delayrange))]
#             HOMI = [f.result() for f in futures]
# 
# 
#         HOMI = HOMI/norm
# 
#         # determine visibility
#         vis = np.abs((np.max(HOMI)-np.min(HOMI))/(np.max(HOMI)))
# 
#         # calc FWHM
#         delayrangeneg = delayrange[:int(np.floor(len(delayrange)/2))]
#         delayrangepos = delayrange[int(np.floor(len(delayrange)/2)):]
#         HOMIneg = HOMI[:int(np.floor(len(HOMI)/2))]
#         HOMIpos = HOMI[int(np.floor(len(HOMI)/2)):]
#         visinterpolfneg = scipy.interpolate.interp1d(delayrangeneg, HOMIneg-0.25, fill_value='extrapolate')
#         visinterpolfpos = scipy.interpolate.interp1d(delayrangepos, HOMIpos-0.25, fill_value='extrapolate')
#         negrootstartest=delayrangeneg[int(np.floor(len(delayrangeneg)/2))]
#         posrootstartest=delayrangepos[int(np.floor(len(delayrangepos)/2))]
#         negroot=scipy.optimize.fsolve(visinterpolfneg, negrootstartest)
#         posroot=scipy.optimize.fsolve(visinterpolfpos, posrootstartest)
# 
#         homfwhm=posroot[0]-negroot[0]
# 
#         t1=datetime.now()
#         print('calculating HOM took', (t1-t0).total_seconds(), 's')
#         return [HOMI,vis,homfwhm]
    
    #by numerical integration
    def getHOMinterference(self, pwl, temp, polingp, qpmorder, tau, cl, signalrange, idlerrange,JSIresolution, pumpshape, delayrange, homphase, refidxfunc, filterfuncs):
        t0=datetime.now()
        [self.nx, self.ny, self.nz] = refidxfunc
        #
        # https://arxiv.org/pdf/1211.0120.pdf (On the Purity and Indistinguishability of Down-Converted Photons. Osorio, Sangouard, thew 2012)
        # Ansari, 2013 msc thesis
        #
        X, Y = np.meshgrid(signalrange, idlerrange)
        
        jsa1 = self.JSAsech(pwl, X, Y, tau, temp, polingp, cl)
        jsa2 = self.JSAsech(pwl, Y, X, tau, temp, polingp, cl)
        jsa2 = np.conjugate(jsa2[::-1,::-1].T)
        #jsa2 = np.conjugate(jsa2)
        jsi = jsa1*np.conjugate(jsa1)
        norm = np.sum(jsi)
    
        def homf(i):
            return  np.sum(scipy.integrate.simpson(jsa1*jsa2*np.exp(1j*2*np.pi*Constants().c*(1/signalrange-1/idlerrange.T)*delayrange[i])))
        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
            futures = [ex.submit(homf, i) for i in range(0,len(delayrange))]
            HOMI = [f.result() for f in futures]

        HOMI = 0.5 - 0.5 * np.abs(HOMI/norm)
        
        # determine visibility
        vis = np.abs((np.max(HOMI)-np.min(HOMI))/(np.max(HOMI)))

        # calc FWHM
        delayrangeneg = delayrange[:int(np.floor(len(delayrange)/2))]
        delayrangepos = delayrange[int(np.floor(len(delayrange)/2)):]
        HOMIneg = HOMI[:int(np.floor(len(HOMI)/2))]
        HOMIpos = HOMI[int(np.floor(len(HOMI)/2)):]
        visinterpolfneg = scipy.interpolate.interp1d(delayrangeneg, HOMIneg-0.25, fill_value='extrapolate')
        visinterpolfpos = scipy.interpolate.interp1d(delayrangepos, HOMIpos-0.25, fill_value='extrapolate')
        negrootstartest=delayrangeneg[int(np.floor(len(delayrangeneg)/2))]
        posrootstartest=delayrangepos[int(np.floor(len(delayrangepos)/2))]
        negroot=scipy.optimize.fsolve(visinterpolfneg, negrootstartest)
        posroot=scipy.optimize.fsolve(visinterpolfpos, posrootstartest)

        homfwhm=posroot[0]-negroot[0]

        t1=datetime.now()
        print('calculating HOM took', (t1-t0).total_seconds(), 's')
        return [HOMI,vis,homfwhm]
    
    def getFWHMvstau(self, pwl, signalrange, idlerrange, temp, polingp, qpmorder, cl, taurange, refidxfunc, filterfuncs, JSIresolution, pumpshape, decprec, usetaucf):
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
                self.filtersignalfunction = filterfuncs[0]
                self.filteridlerfunction = filterfuncs[1]
                if not (self.filtersignalfunction==None and self.filteridlerfunction==None):
                    self.filtermatrix = []
                    for i in range(0, len(signalrange)):
                        filtervector = []
                        for j in range(0, len(idlerrange)):
                            filterval = self.filteridlerfunction(signalrange[i]) * self.filtersignalfunction(idlerrange[j])
                            filtervector.append(filterval)
                        self.filtermatrix.append(filtervector)
                    
                    result = result*self.filtermatrix
            hmptss=[]
            hmptsi=[]
            hmpts=[]
            
            result=result.round(decprec)
            idcs=np.where(result==0) 
            for i in range(0,len(idcs[0])):
                    hmptss.append(signalrange[idcs[0][i]]*10**9)
                    hmptsi.append(idlerrange[idcs[1][i]]*10**9)
                    
            if len(hmptss) != 0:
                if len(hmptsi) !=0:
                    fwhmsig.append(max(hmptss)-min(hmptss))
                    fwhmid.append(max(hmptsi)-min(hmptsi))
            
        return [fwhmsig,fwhmid]
