#!/usr/bin/python3

import numpy
import scipy
import scipy.optimize

class PMC():
        def __init__(self): 
            # Thermal expansion coefficients of KTP
            # (Emanueli 2003)
            self.TXCa = 6.7 * 10 ** (-6)
            self.TXCb = 11 * 10 ** (-9)
            self.TXrefT = 25

        # thermal expansion factor
        def thermexpfactor(self, T):
            return (1 + self.TXCa * (T - self.TXrefT) + self.TXCb * (T - self.TXrefT) * (T - self.TXrefT))
        
        #phasematching conditions
        def econv(self,ls, li):
                return 1/lp - 1/ls - 1/li
        def pconv(self,ls, li, T, PP):
                pp=ny(lp,T)/lp
                ps=ny(ls,T)/ls
                pi=nz(li,T)/li
                pc=m/PP
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
                return scipy.optimize.newton_krylov(self.epconvonlywl(x[1],x[2]),[2*lp,2*lp],f_tol=1e-9)

        #returns a function that only depends on the poling period
        def wlgaponlyT(self,PP,lp):
                def wlgap2(T):
                        swl,iwl=self.SIwls([lp,T,PP])
                        return swl-iwl
                return wlgap2
        
        #calculate signal and idler wavelengths
        def getSI_wl(self,pumpwl,PP,Trange):
                sigwl=numpy.zeros(len(Trange))
                idwl=numpy.zeros(len(Trange))
                sigwltherm=numpy.zeros(len(Trange))
                idwltherm=numpy.zeros(len(Trange))
                txf=1
                for i in range(0,len(Trange)):
                    txf = self.thermexpfactor(Trange[i])
                    [sigwl[i],idwl[i]]=self.SIwls([pumpwl,Trange[i],PP*txf])

                #calculate the crossing point temperature
                Tcp=0
                Tcpguess=50
                Tcp=scipy.optimize.fsolve(wlgaponlyT(PP,lp),Tcguess)
                
                #return:
                #signal wavelength, idler wavelength, crossing point temperature
                return [sigwl,idwl,Tcp[0]]
