#!/usr/bin/python3

import numpy
import scipy
import scipy.optimize

class PMC:
        def __init__(self):
            # Thermal expansion coefficients of KTP
            # (Emanueli 2003)
            self.TXCa = 6.7 * 10 ** (-6)
            self.TXCb = 11 * 10 ** (-9)
            self.TXrefT = 25

            #default values
            self.lp=775 * 10 ** (-9)
            self.PP=46.2 * 10 ** (-6)
            self.m=1

        # thermal expansion factor
        def thermexpfactor(self, T):
            return (1 + self.TXCa * (T - self.TXrefT) + self.TXCb * (T - self.TXrefT) * (T - self.TXrefT))

        #phasematching conditions
        def econv(self,ls, li):
                return 1/self.lp - 1/ls - 1/li
        def pconv(self,ls, li, T, PP):
                pp=self.ny(self.lp,T)/self.lp
                ps=self.ny(ls,T)/ls
                pi=self.nz(li,T)/li
                pc=self.m/PP
                return pp-ps-pi-pc

        #this wrapper returns a list of functions(phasematching conditions) that only takes signal- and idler-wavelength as arguments
        def epconvonlywl(self,T,PP):
                def epconv(x):
                        #x[0]:lambda_s
                        #x[1]:lambda_i
                        return [self.econv(x[0],x[1]), self.pconv(x[0],x[1],T,PP)]
                return epconv

        #returns signal and idler wavelengths that satisfy phasematching conditions for a given pumpwavelength, temperature and poling period
        def SIwls(self,x):
                #x[0]: lambda_pump
                #x[1]: Temperature
                #x[2]: Poling period
                #return scipy.optimize.newton_krylov(self.epconvonlywl(x[1],x[2]),[2*x[0],2*x[0]],f_tol=1e-6) #slower!
                return scipy.optimize.fsolve(self.epconvonlywl(x[1], x[2]), [2 * x[0], 2 * x[0]], xtol=1e-6)

        #returns a function that only depends on the poling period
        def wlgaponlyT(self,PP,lp):
                def wlgap2(T):
                        swl,iwl=self.SIwls([lp,T,PP])
                        return swl-iwl
                return wlgap2

        # returns a function that only depends on the poling period
        def wlgaponlyPP(self, T, lp):
                def wlgap2(PP):
                        swl, iwl = self.SIwls([lp, T, PP])
                        return swl - iwl
                return wlgap2

        #calculate signal and idler wavelengths. Temperature is variated
        def getSI_wl_varT(self,pumpwl,polingp,Trange,refidxfunc,qpmorder):

                self.m=qpmorder
                self.lp=pumpwl
                self.PP = polingp

                self.nx = refidxfunc[0]
                self.ny = refidxfunc[1]
                self.nz = refidxfunc[2]

                sigwl=numpy.zeros(len(Trange))
                idwl=numpy.zeros(len(Trange))
                txf = self.thermexpfactor(Trange)
                for i in range(0,len(Trange)):
                        [sigwl[i],idwl[i]]=self.SIwls([pumpwl,Trange[i],polingp*txf[i]])

                #calculate the crossing point temperature
                Tcp=0
                Tcpguess=50
                Tcp=scipy.optimize.fsolve(self.wlgaponlyT(polingp,pumpwl),Tcpguess)

                #return:
                #signal wavelength, idler wavelength, crossing point temperature
                return [sigwl,idwl,Tcp[0]]

        # calculate signal and idler wavelengths. PP is variated
        def getSI_wl_varPP(self, pumpwl, PPrange, T, refidxfunc, qpmorder):

                self.m = qpmorder
                self.lp = pumpwl

                self.nx = refidxfunc[0]
                self.ny = refidxfunc[1]
                self.nz = refidxfunc[2]
                self.T = T

                sigwl = numpy.zeros(len(PPrange))
                idwl = numpy.zeros(len(PPrange))
                txf = self.thermexpfactor(PPrange)
                for i in range(0, len(PPrange)):
                        [sigwl[i], idwl[i]] = self.SIwls([pumpwl, T, PPrange[i] * txf[i]])

                # calculate the crossing point temperature
                PPcp = 0
                PPguess = ((PPrange[-1]-PPrange[0])/2)
                PPcp = scipy.optimize.fsolve(self.wlgaponlyPP(T, pumpwl), PPguess)

                # return:
                # signal wavelength, idler wavelength, crossing point temperature
                return [sigwl, idwl, PPcp[0]]
