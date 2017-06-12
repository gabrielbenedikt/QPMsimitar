#!/usr/bin/python3

import numpy, scipy, scipy.optimize
import matplotlib
import matplotlib.pyplot as plt
from Constants import Constants
import array
from matplotlib.offsetbox import AnchoredText
import os

class JSI():
    def __init__(self):
        # Thermal expansion coefficients of KTP
        # (Emanueli 2003)
        self.TXCa = 6.7 * 10 ** (-6)
        self.TXCb = 11 * 10 ** (-9)
        self.TXrefT = 25

        #signal and idler wavelength filtering
        self.useFilter = False  # change filterfunction(ls,li) below
        self.filtertype = 'rectangular' #type of filter used. TODO: Implement different filter functions
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

        # resolution
        # number of subintervals the s/i wavelength range will be devided in
        self.numpts = 500
        # wavelengthrange to calculate [[a,b],[c,d]]
        # will calculate signal from 2*pwl-a to 2*pwl+b
        # will calculate idler from 2*pwl-c to 2*pwl+d
        self.wlrange = [[5, 5], [5, 5]]

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

    # this wrapper returns a n arroy of functions(phasematching conditions) that only takes signal- and idler-wavelength as arguments
    def epconvonlywl(self, lp, T, PP):
        def epconv(x):
            # x[0]:lambda_p
            # x[0]:lambda_s
            # x[1]:lambda_i
            return [self.econv(lp, x[0], x[1]), self.pconv(lp, x[0], x[1], T, PP)]

        return epconv

    # returns signal and idler wavelengths that satisfy phasematching conditions for a given pumpwavelength, temperature and poling period
    def SIwls(self, x):
        # x[0]: lambda_pump
        # x[1]: Temperature
        # x[2]: Poling period
        return scipy.optimize.fsolve(self.epconvonlywl(x[0], x[1], x[2]), [2 * x[0], 2 * x[0]], xtol=1e-6)

    # returns difference between signal and idler wavelength for a given pumpwavelength, temperature and poling period
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

    # phase matching intensity
    def PMIsech(self, dk, cl):
        return numpy.square(self.PMA(dk, cl))

    def PMIgauss(self, dk, cl):
        return numpy.square(self.PMA(dk, cl))

    ###################
    ###         gaussian         ###
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
        tau = tauac * Constants().taucfgauss
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
        tau = tauac * Constants().taucfgauss
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
    ###          sech^2          ###
    ###################
    # pump envelope amplitude for sech^2 beam
    def PEAsech(self, lp, ls, li, B):
        wfact = 2 * Constants().pi * Constants().c * (1 / ls + 1 / li - 1 / lp)
        argument = wfact * B
        return (1 / numpy.cosh(argument))

    # pump envelope intensity for gaussian beam
    def PEIsech(self, lp, ls, li, dw):
        return (self.PEAsech(lp, ls, li, dw)) ** 2

    # joint spectral amplitude for sech^2 beam
    def JSAsech(self, lp, ls, li, tauac, t, pp, cl):
        tau = tauac * Constants().taucfsech
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
        tau = tauac * Constants().taucfsech
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

    def getplots(self,pumpwl,signalrange,idlerrange,tau,temp,polingp,crystallength,refidxfunc,qpmorder,filter,plotJSI,resolution):
        #pumpwl: Pump wavelength
        #signalrange: [double,double]: Signal wavelength range
        #idlerrange: [double,double]: Idler wavelength range
        #tau: pump pulse duration
        #temp: Temperature
        #polingp: Crystal poling period
        #crystallength: Length of crystal
        #refidxfunc: [nx,ny,nz]: Functions for refractive indices of crystal
        #qpmorder: Quasi phase matching order
        #filter: [string,bool,bool]: [Type of filter to use, True to use filter for signal, True to use filter for idler]
        #plotJSI: bool: True for JSI, false for JSA
        #resolution: number of points signalrange and idlerrange will be devided into

        #get signalrange!
        #get idlerrange!
        self.nx = refidxfunc[0]
        self.ny = refidxfunc[1]
        self.nz = refidxfunc[2]

        self.pwl = pumpwl
        self.tau = tau
        self.m = qpmorder
        self.T = temp
        self.PP = polingp * self.thermexpfactor(self.T)
        self.L = crystallength * self.thermexpfactor(self.T)

        self.useFilter=True
        if (filter[0] == 'none'):
            self.useFilter=False
        else:
            if (filter[0] == 'rectangular'):
                self.filtertype = 'rectangular'
            elif (filter[0] == 'gaussian'):
                self.filtertype = 'gaussian'

            if filter[1]==True:
                self.filtersignal=True
            else:
                self.filtersignal=False

            if filter[2]==True:
                self.filteridler=True
            else:
                self.filteridler=False

            if (filter[1]==False):
                if (filter[2]==False):
                    print ('ERROR: filtertype not none but neither used for idler nor signal.')

        if plotJSI==False:
            self.calcJSA = True
            self.calcJSI = False
        else:
            self.calcJSA = False
            self.calcJSI = True

        self.numpts=resolution



    #################################
    ###           JSI             ###
    #################################
    signalrange = numpy.linspace(2 * pwl - wlrange[0][0] * 10 ** (-9), 2 * pwl + wlrange[0][1] * 10 ** (-9), numpts)
    sigpoints = len(signalrange)
    idlerrange = numpy.linspace(2 * pwl - wlrange[1][0] * 10 ** (-9), 2 * pwl + wlrange[1][1] * 10 ** (-9), numpts)
    idlpoints = len(idlerrange)

    X, Y = numpy.meshgrid(signalrange, idlerrange)
    filtermatrix = []
    for i in range(0, len(signalrange)):
        filtervector = []
        for j in range(0, len(idlerrange)):
            filtervector.append(filterfunction(signalrange[i], idlerrange[j]))
        filtermatrix.append(filtervector)

    if calcJSA:
        print("Calculating PEA, PMA and JSA for gaussian beam shape...")
        [Zpeig, Zpmig, Zjsig] = PEAnPMAnJSAgauss(pwl, X, Y, tau, T, PP, L)
        print("Calculating PEA, PMA and JSA for sech2 beam shape...")
        [Zpeis, Zpmis, Zjsis] = PEAnPMAnJSAsech(pwl, X, Y, tau, T, PP, L)
        if useFilter:
            Zjsag = Zjsag * filtermatrix
            Zjsas = Zjsas * filtermatrix
    elif calcJSI:
        print("Calculating PEI, PMI and JSI for gaussian beam shape...")
        [Zpeig, Zpmig, Zjsig] = PEInPMInJSIgauss(pwl, X, Y, tau, T, PP, L)
        print("Calculating PEI, PMI and JSI for sech2 beam shape...")
        [Zpeis, Zpmis, Zjsis] = PEInPMInJSIsech(pwl, X, Y, tau, T, PP, L)
        if useFilter:
            Zjsig = Zjsig * filtermatrix
            Zjsis = Zjsis * filtermatrix

    print("plotting..")
    f, ((peigplt, pmigplt, jsigplt), (peisplt, pmisplt, jsisplt)) = plt.subplots(2, 3, sharex='col', sharey='row',
                                                                                 figsize=(16, 11))
    colormap = matplotlib.cm.jet
    # axes range
    xmin = numpy.min(signalrange) * 10 ** 9
    xmax = numpy.max(signalrange) * 10 ** 9
    ymin = numpy.min(idlerrange) * 10 ** 9
    ymax = numpy.max(idlerrange) * 10 ** 9
    # prepare subplots
    ppeig = peigplt.imshow(Zpeig, cmap=colormap, vmin=abs(Zpeig).min(), vmax=abs(Zpeig).max(), aspect='auto',
                           origin='lower', interpolation='none', extent=[xmin, xmax, ymin, ymax])
    ppmig = pmigplt.imshow(Zpmig, cmap=colormap, vmin=abs(Zpmig).min(), vmax=abs(Zpmig).max(), aspect='auto',
                           origin='lower', interpolation='none', extent=[xmin, xmax, ymin, ymax])
    # pjsig = jsigplt.imshow(Zjsig, cmap=colormap, vmin=abs(Zjsig).min(), vmax=abs(Zjsig).max(),aspect='auto', origin='lower',interpolation='none', extent=[xmin,xmax,ymin,ymax])
    pjsig = jsigplt.imshow(Zjsig, cmap=colormap, vmin=abs(Zjsig).min(), vmax=abs(Zjsig).max(), aspect='auto',
                           origin='lower', interpolation='none', extent=[xmin, xmax, ymin, ymax])
    ppeis = peisplt.imshow(Zpeis, cmap=colormap, vmin=abs(Zpeis).min(), vmax=abs(Zpeis).max(), aspect='auto',
                           origin='lower', interpolation='none', extent=[xmin, xmax, ymin, ymax])
    peisplt.set_xlabel(r'$\lambda_s$ [nm]')
    peisplt.set_ylabel(r'$\lambda_i$ [nm]')
    ppmis = pmisplt.imshow(Zpmis, cmap=colormap, vmin=abs(Zpmis).min(), vmax=abs(Zpmis).max(), aspect='auto',
                           origin='lower', interpolation='none', extent=[xmin, xmax, ymin, ymax])
    # pjsis = jsisplt.imshow(Zjsis, cmap=colormap, vmin=abs(Zjsis).min(), vmax=abs(Zjsis).max(),aspect='auto', origin='lower',interpolation='none', extent=[xmin,xmax,ymin,ymax])
    pjsis = jsisplt.imshow(Zjsis, cmap=colormap, vmin=abs(Zjsis).min(), vmax=abs(Zjsis).max(), aspect='auto',
                           origin='lower', interpolation='none', extent=[xmin, xmax, ymin, ymax])
    # size of axes label
    for plot in [peigplt, pmigplt, jsigplt, peisplt, pmisplt, jsisplt]:
        plot.tick_params(axis='both', which='major', labelsize='medium')
        plot.tick_params(axis='both', which='minor', labelsize='medium')
    # label plot
    if calcJSA:
        peigplt.set_title('PEA', fontsize=20)
        pmigplt.set_title('PMA', fontsize=20)
        jsigplt.set_title('JSA', fontsize=20)
        # f.suptitle('PEA, PMA and JSA for gaussian and sech2 beamshape',fontsize=24)
    elif calcJSI:
        peigplt.set_title('PEI', fontsize=20)
        pmigplt.set_title('PMI', fontsize=20)
        jsigplt.set_title('JSI', fontsize=20)
        f.suptitle('PEI, PMI and JSI for gaussian and sech2 beamshape', fontsize=24)
    # create legend
    f.subplots_adjust(bottom=0.2)
    cbar_ax = f.add_axes([0.05, 0.05, 0.9, 0.025])
    cbar = f.colorbar(ppeig, cax=cbar_ax, orientation='horizontal')
    cbar.set_label('a.u.', fontsize='medium', labelpad=-1)
    # label rows with beamshape
    pad = 20
    jsigplt.annotate('gaussian', xy=(1, 0.5), xytext=(jsigplt.yaxis.labelpad + pad, 0),
                     xycoords='axes fraction', textcoords='offset points',
                     ha='right', va='center', rotation=270, fontsize=20)
    jsisplt.annotate('sech2', xy=(1, 0.5), xytext=(jsisplt.yaxis.labelpad + pad, 0),
                     xycoords='axes fraction', textcoords='offset points',
                     ha='right', va='center', rotation=270, fontsize=20)
    # state additional paramters on plot
    parameterstring = 'Pump wavelength: {0:.2f}nm\n'.format(pwl * 10 ** 9) + 'Crystal Length: {0:.2f}mm\n'.format(
        L * 10 ** 3) + 'Poling period: {0:.2f}µm\n'.format(PP * 10 ** 6) + 'Temperature: {0:.2f}°C\n'.format(
        T) + 'Pulse duration: {0:.2f}ps'.format(tau * 10 ** 12)
    plt.annotate(parameterstring, xy=(0.005, 0.885), xycoords='figure fraction', fontsize=10, color='r')
    plt.tight_layout(pad=8, w_pad=0.5, h_pad=2)
    if saveplot:
        print("saveing pdf...")
        plt.savefig(plotsavename, format="pdf", transparent=True, bbox_inches='tight', pad_inches=0)
    print("showing plot...")
    plt.show()