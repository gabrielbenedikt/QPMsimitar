#!/usr/bin/python3
#this script calculates the JSA for a sech^2 shaped pump, applies a SDV to it and calculates the purity
#it does so for different pump pulse widths and plots purity vs. pulse width

##################################
#                           TODO                                      ###
##################################
#
###
### Include thermal expansion leading to PP change
###
#
#

import numpy, scipy, scipy.optimize, scipy.linalg, scipy.interpolate
import matplotlib
import matplotlib.pyplot as plt
import array
from matplotlib.offsetbox import AnchoredText
import os
import pprint

#TESTING
from time import time
#/TESTING

useabs=False#take absolute value of pea?

#use spectral filter?
usefilter=False #change filterfunction(l) below
FilterString='Filter used: $\lambda_c$={0:.2f}nm, bandwidth={1:.2f}nm. (Only on idler)\n'.format(1546.2,3)

includegauss=True#do calc and plot for gaussian pulse too
includeintensity=False#do calc and plot for intensity too
plot=True
saveplot=True
#usetauAC: if true, treat pulsewidth as ACF-FWHM and convert to pulseduration FWHM using a deconvolution factor, depending on pulseshape
usetauAC=True

source=1

#savename of the plot
plotsavename="purity_source{0:}_nofilter.pdf".format(source)
#plotsavename="purity.pdf"
plottitle='Purity (source {0:})'.format(source)
#plottitle='Purity'

#Pump wavelength in m
pwl=773.1*10**(-9)
#Quasi phase matching order
m = -1
#Poling period in m
PPvec=[46.24657*10**(-6),46.24156*10**(-6)]
PP=PPvec[source-1]
#PP=46.03*10**(-6)
#PP=46.1*10**(-6)

#Temperature in °C
Tvec=[23.7,25.15]
T=Tvec[source-1]
#T=50
#T=32
#Crystal length in m
L=30*10**(-3)

#pump pulsewidth in s (measured by intensity autocorrelator)
#tauinter=[1.0*10**(-12),12.0*10**(-12)]
tauinter=[0.1*10**(-12),5.0*10**(-12)]
numtau=100
#factors to convert pulsewidth between autocorrelator measurement and real pulsewidth
if usetauAC:
    taucfgauss=(1/numpy.sqrt(2))     # gaussian shape
    taucfsech=0.647                         # sech2 shape
else:
    taucfgauss=1     # gaussian shape
    taucfsech=1                         # sech2 shape

#resolution
#number of subintervals the s/i wavelength range will be devided in
numpts=200
#wavelengthrange to calculate (interval centered around 2*pwl)
wlrange=8

#
#constants
#
#Pi
pi=numpy.pi
#time bandwidth product
tbwpgauss=0.441
tbwpsech=((2*numpy.log(1+numpy.sqrt(2)))/(pi))**2
print('TBWP (gauss): ', tbwpgauss)
print('TBWP (sech): ', tbwpsech)
#speed of light in m/s
c = 299792458
#refractive indices fit parameters
#
#
a10z = 9.9587*10**(-6)
a11z = 9.9228*10**(-6)
a12z = -8.9603*10**(-6)
a13z = 4.1010*10**(-6)
a20z = -1.1882*10**(-8)
a21z = 10.459*10**(-8)
a22z = -9.8136*10**(-8)
a23z = 3.1481*10**(-8)
a10y = 6.2897*10**(-6)
a11y = 6.3061*10**(-6)
a12y = -6.0629*10**(-6)
a13y = 2.6486*10**(-6)
a20y = -0.14445*10**(-8)
a21y = 2.2244*10**(-8)
a22y = -3.5770*10**(-8)
a23y = 1.3470*10**(-8)

def ny(lin, t):
        #formulas for µm...
        l=lin*10**(6)#formulas for µm...
        f1=a10y+(1/l)*(a11y+(1/l)*(a12y+(a13y/l) ) )
        f2=a20y+(1/l)*(a21y+(1/l)*(a22y+(a23y/l) ) )
        f3=numpy.sqrt(numpy.abs(-0.0138408*l*l+0.922683/(1-0.0467695/(l*l))+2.0993))
        return (t-25)*f1+(t-25)*(t-25)*f2+f3

def nz(lin, t):
        #formulas for µm...
        l=lin*10**(6)#formulas for µm...
        f1=a10z+(1/l)*(a11z+(1/l)*(a12z+(a13z/l) ) )
        f2=a20z+(1/l)*(a21z+(1/l)*(a22z+(a23z/l) ) )
        f3=numpy.sqrt(numpy.abs(-0.00968956*l*l+1.18431/(1-0.0514852/(l*l))+0.6603/(1-100.00507/(l*l))+2.12725))
        return (t-25)*f1+(t-25)*(t-25)*f2+f3

#
#Thermal expansion coefficients of KTP
#(Emanueli 2003)
#
TXCa=6.7*10**(-6)
TXCb=11*10**(-9)
TXrefT=25
#thermal expansion factor
def thermexpfactor(T):
        return (1+TXCa*(T-TXrefT)+TXCb*(T-TXrefT)*(T-TXrefT))

#phasematching conditions
def econv(lp, ls, li):
        return 1/lp - 1/ls - 1/li
def pconv(lp, ls, li, T, PP):
        #wavelengths in nm, but refractive-index-functions take µm
        pp=ny(lp,T)/lp
        ps=ny(ls,T)/ls
        pi=nz(li,T)/li
        pc=m/PP
        return pp-ps-pi-pc

#this wrapper returns a n arroy of functions(phasematching conditions) that only takes signal- and idler-wavelength as arguments
def epconvonlywl(lp,T,PP):
        def epconv(x):
                #x[0]:lambda_p
                #x[0]:lambda_s
                #x[1]:lambda_i
                return [econv(lp,x[0],x[1]), pconv(lp,x[0],x[1],T,PP)]
        return epconv

#returns signal and idler wavelengths that satisfy phasematching conditions for a given pumpwavelength, temperature and poling period
def SIwls(x):
        #x[0]: lambda_pump
        #x[1]: Temperature
        #x[2]: Poling period
        return scipy.optimize.fsolve(epconvonlywl(x[0],x[1],x[2]),[2*x[0],2*x[0]],xtol=1e-14)

#returns difference between signal and idler wavelength for a given pumpwavelength, temperature and poling period
def wlgap(x):
        #x[0]: lambda_pump
        #x[1]: Temperature
        #x[2]: Poling period
        swl,iwl=SIwls([x[0],x[1],x[2]])
        return swl-iwl

#returns a function that only depends on the poling period
def wlgaponlyT(lp, PP):
        def wlgap2(T):
                swl,iwl=SIwls([lp,T,PP])
                return swl-iwl
        return wlgap2

def deltak(lp, ls, li, T, PP):
        tmp=pconv(lp, ls, li, T, PP)
        return 2*pi*tmp

#calculates pumpwavelength out of signal and idler wavelength, using energy conservation
def lambdap(ls, li):
        return 1/(1/ls+1/li)


#phase matching amplitude
def PMA(dk,cl):
        #note: numpy.sinc(x) evaluates to Sin(pi*x)/(pi*x)
        #this really produced some headache.
        return numpy.sinc(dk*cl/(2*pi))#*numpy.exp(numpy.complex(0,1)*dk*cl/2)
def PMAgauss(dk, cl):
        return PMA(dk,cl)
def PMAsech(dk, cl):
        return PMA(dk,cl)
#phase matching intensity
def PMIsech(dk, cl):
        return numpy.square(PMA(dk,cl))
def PMIgauss(dk, cl):
        return numpy.square(PMA(dk,cl))

###################
###         gaussian         ###
###################
#pump envelope (gauss)
#amplitude
def PEAgauss(lp, ls, li, sp):
        dl = 1/ls + 1/li - 1/lp
        return (   numpy.exp(   - (   pi*c*( dl ) / (  sp )   )**2   )   ) #note: 2*pi*c/(2*sp)
#intensity
def PEIgauss(lp, ls, li, sp):
        return PEAgauss(lp, ls, li, sp)**2

#joint spectral amplitude for gaussian beam
def JSAgauss(lp, ls, li, tauac, t, pp, cl):
        #tauac: autocorrelator measured pulsewidth
        tau=tauac*taucfgauss
        #tau=tauac
        dnu=tbwpgauss/tau # FWHM in frequency
        dw=2*pi*dnu # FWHM in angular frequency
        sp=dw/(2*numpy.sqrt(2*numpy.log(2))) # gaussian standard deviation from FWHM
        dk=deltak(lambdap(ls,li), ls, li, t, pp)
        jsa=PEAgauss(lp,ls,li,sp)*PMAgauss(dk,cl)
        if useabs:
            return numpy.absolute(jsa)
        else:
            return jsa
#joint spectral intensity for gaussian beam
def JSIgauss(lp, ls, li, tauac, t, pp, cl):
        return (JSAgauss(lp, ls, li, tauac, t, pp, cl))**2

def PEAnPMAnJSAgauss(lp, ls, li, tauac, t, pp, cl):
        #tauac: autocorrelator measured pulsewidth
        tau=tauac*taucfgauss
        #tau=tauac
        dnu=tbwpgauss/tau # FWHM in frequency
        dw=2*pi*dnu # FWHM in angular frequency
        sp=dw/(2*numpy.sqrt(2*numpy.log(2))) # gaussian standard deviation from FWHM
        dk=deltak(lambdap(ls,li), ls, li, t, pp)
        pea=PEAgauss(lp,ls,li,sp)
        pma=PMAgauss(dk,cl)
        jsa=pea*pma
        if useabs:
            return [numpy.absolute(pea),numpy.absolute(pma),numpy.absolute(jsa)]
        else:
            return [pea,pma,jsa]
def PEInPMInJSIgauss(lp, ls, li, tauac, t, pp, cl):
        [pei,pea,jsa]=PEAnPMAnJSAgauss(lp,ls,li,tauac,t,pp,cl)
        return [pei**2,pea**2,jsa**2]


###################
###          sech^2          ###
###################
#pump envelope amplitude for sech^2 beam
def PEAsech(lp, ls, li, B):
        wfact=2*pi*c*(1/ls+1/li-1/lp)
        argument=wfact*B
        return (1/numpy.cosh(argument))
#pump envelope intensity for gaussian beam
def PEIsech(lp, ls, li, dw):
        return (PEAsech(lp, ls, li, dw))**2

#joint spectral amplitude for sech^2 beam
def JSAsech(lp, ls, li, tauac, t, pp, cl):
        tau = tauac*taucfsech
        #tau=tauac
        dw=2*pi*tbwpsech/(tau)
        B = 2*numpy.arccosh(numpy.sqrt(2))/dw
        dk=deltak(lambdap(ls,li), ls, li, t, pp)
        jsa=PEAsech(lp, ls, li, B)*PMAsech(dk, cl)
        if useabs:
            return numpy.absolute(jsa)
        else:
            return jsa
#joint spectral intensity for sech^2 beam
def JSIsech(lp, ls, li, tauac, t, pp, cl):
        return JSAsech(lp, ls, li, tauac, t, pp, cl)**2

def PEAnPMAnJSAsech(lp,ls,li,tauac,t,pp,cl):
        tau = tauac*taucfsech
        #tau=tauac
        dw=2*pi*tbwpsech/(tau)
        B = 2*numpy.arccosh(numpy.sqrt(2))/dw
        dk=deltak(lambdap(ls,li), ls, li, t, pp)
        pea=PEAsech(lp, ls, li, B)
        pma=PMAsech(dk, cl)
        jsa=pea*pma
        if useabs:
            return [numpy.absolute(pea),numpy.absolute(pma),numpy.absolute(jsa)]
        else:
            return [pea,pma,jsa]
def PEInPMInJSIsech(lp,ls,li,tauac,t,pp,cl):
        [pea,pma,jsa]=PEAnPMAnJSAsech(lp,ls,li,tauac,t,pp,cl)
        return [pea**2,pma**2,jsa**2]

def filterfunction(ls,li):
    #rectangular filter with:
    #lc......central wavelength
    #bwh..half bandwidth (such that transmittivity=1 in the range (lc-bwh,lc+bwh) ).
    lc=1546.2*10**-(9)
    bwh=1.5*10**(-9)
    
    #rectangular filter on both signal and idler
    if 0:
        if (ls<(lc-bwh) or ls>(lc+bwh) or li<(lc-bwh) or li>(lc+bwh)):
            return(0)
        else:
            return(1)
    #rectangular filter only on signal
    if 0:
        if (ls<(lc-bwh) or ls>(lc+bwh)):
            return(0)
        else:
            return(1)
    #rectangular filter only on idler
    if 1:
        if (li<(lc-bwh) or li>(lc+bwh)):
            return(0)
        else:
            return(1)


print('Calculating purity of source {0:}'.format(source))
print('Pump wavelength: {0:.2f}nm'.format(pwl*10**9))
print('Poling period: {0:.4f}µm'.format(PP*10**6))
print('Temperature: {0:.2f}°C'.format(T))
print('Crystal length: {0:.2f}mm'.format(L*10**3))
print('Pulsewidth range: {0:.2f}-{1:.2f} ps'.format(tauinter[0]*10**12,tauinter[1]*10**12))
print('Number of points in pulsewidth range: {0:.0f}'.format(numtau))
print('Wavelength range: {0:.2f}nm'.format(wlrange))
print('Number of points in wavelength range: {0:.0f}'.format(numpts))
print('Use spectral filtering: ',usefilter)
print('Do calculations for gaussian pulse shape too: ',includegauss)
print('Do calculations for intensity too: ',includeintensity)
print('Plot purity vs pulsewidth: ',plot)
print('Save plot: ',saveplot)
print('Use ACF-FWHM and convert to pulsewidth: ',usetauAC)
print('')

PP=PP*thermexpfactor(T)
L=L*thermexpfactor(T)
taurange=numpy.linspace(tauinter[0],tauinter[1],numtau)
signalrange=numpy.linspace(2*pwl-0.5*wlrange*10**(-9),2*pwl+0.5*wlrange*10**(-9),numpts)
sigpoints=len(signalrange)
idlerrange=numpy.linspace(2*pwl-0.5*wlrange*10**(-9),2*pwl+0.5*wlrange*10**(-9),numpts)
X,Y=numpy.meshgrid(signalrange,idlerrange)

filtermatrix=[]
if usefilter:
    for i in range(0,len(signalrange)):
        filtervector=[]
        for j in range(0,len(idlerrange)):
            filtervector.append(filterfunction(signalrange[i],idlerrange[j]))
        filtermatrix.append(filtervector)

idlpoints=len(idlerrange)
purityA=[]
purityI=[]
purityAg=[]
purityIg=[]
for i in range(0,len(taurange)):
        tau=taurange[i]
        if 1:
                #################################
                ###                          JSA                                     ###
                #################################
                #[ZpeasA,ZpmasA,ZjsasA]=PEAnPMAnJSAsech(pwl,X,Y,tau,T,PP,L)
                ZjsasA=JSAsech(pwl, X, Y, tau, T, PP, L)
                if includeintensity:
                    #[ZpeasI,ZpmasI,ZjsasI]=PEInPMInJSIsech(pwl,X,Y,tau,T,PP,L)
                    ZjsasI=JSIsech(pwl,X,Y,tau,T,PP,L)
                if usefilter:
                    ZjsasA=ZjsasA*filtermatrix
                    if includeintensity:
                        ZjsasI=ZjsasI*filtermatrix
                if includegauss:
                        #[ZpeagA,ZpmagA,ZjsagA]=PEAnPMAnJSAgauss(pwl,X,Y,tau,T,PP,L)
                        ZjsagA=JSAgauss(pwl,X,Y,tau,T,PP,L)
                        if includeintensity:
                            #[ZpeagI,ZpmagI,ZjsagI]=PEInPMInJSIgauss(pwl,X,Y,tau,T,PP,L)
                            ZjsagI=JSIgauss(pwl,X,Y,tau,T,PP,L)
                        if usefilter:
                            ZjsagA=ZjsagA*filtermatrix
                            if includeintensity:
                                ZjsagI=ZjsagI*filtermatrix

                #################################
                ###                       Purity                                    ###
                #################################
                #SDV
                sA=scipy.linalg.svd(ZjsasA,overwrite_a=True,compute_uv=False)
                #singular values in s. need to normalize s to get the schmidt magnitudes
                snA=sA/scipy.linalg.norm(sA,2)
                #sum of squares of schmidt magnitudes is purity
                purityA.append(numpy.sum(snA**4))
                infostring="Pulsewidth: {0:.5f}ps\tpurity(sech): {1:.5f}".format(tau*10**(12),purityA[i])
                if includeintensity:
                    sI=scipy.linalg.svd(ZjsasI,overwrite_a=True,compute_uv=False)
                    snI=sI/scipy.linalg.norm(sI,2)
                    purityI.append(numpy.sum(snI**4))
                    infostring=infostring+"\tpurity(intensity,sech): {0:.5f}".format(purityI[i])
                if includegauss:
                        sAg=scipy.linalg.svd(ZjsagA,overwrite_a=True,compute_uv=False)
                        snAg=sAg/scipy.linalg.norm(sAg,2)
                        purityAg.append(numpy.sum(snAg**4))
                        infostring=infostring+"\tpurity(gauss): {0:.5f}".format(purityAg[i])
                        if includeintensity:
                            sIg=scipy.linalg.svd(ZjsagI,overwrite_a=True,compute_uv=False)
                            snIg=sIg/scipy.linalg.norm(sIg,2)
                            purityIg.append(numpy.sum(snIg**4))
                            infostring=infostring+"\tpurity(intensity,gauss): {0:.5f}".format(purityIg[i])
                print(infostring)

#interpolate purity curve
increased_taurange=numpy.linspace(tauinter[0],tauinter[1],100000)
interPurA=scipy.interpolate.InterpolatedUnivariateSpline(taurange,purityA)
#interPurAsech=scipy.interpolate.InterpolatedUnivariateSpline(taurange*taucfsech,purityA)
interPurAsech=scipy.interpolate.InterpolatedUnivariateSpline(taurange,purityA)
interPurAvals=interPurA(increased_taurange)
interPurAsechvals=interPurAsech(increased_taurange)
maxA=numpy.max(interPurAvals)
maxAidx=numpy.argmax(interPurAvals)
if includeintensity:
    interPurI=scipy.interpolate.InterpolatedUnivariateSpline(taurange,purityI)
    #interPurIsech=scipy.interpolate.InterpolatedUnivariateSpline(taurange*taucfsech,purityI)
    interPurIsech=scipy.interpolate.InterpolatedUnivariateSpline(taurange,purityI)
    interPurIvals=interPurI(increased_taurange)
    interPurIsechvals=interPurIsech(increased_taurange)
    maxI=numpy.max(interPurIvals)
    maxIidx=numpy.argmax(interPurIvals)

if includegauss:
            interPurAg=scipy.interpolate.InterpolatedUnivariateSpline(taurange,purityAg)
            #interPurAgauss=scipy.interpolate.InterpolatedUnivariateSpline(taurange*taucfsech,purityAg)
            interPurAgauss=scipy.interpolate.InterpolatedUnivariateSpline(taurange,purityAg)
            interPurAgvals=interPurAg(increased_taurange)
            interPurAgaussvals=interPurAgauss(increased_taurange)
            maxAg=numpy.max(interPurAgvals)
            maxAgidx=numpy.argmax(interPurAgvals)
            if includeintensity:
                interPurIg=scipy.interpolate.InterpolatedUnivariateSpline(taurange,purityIg)
                #interPurIgauss=scipy.interpolate.InterpolatedUnivariateSpline(taurange*taucfsech,purityIg)
                interPurIgauss=scipy.interpolate.InterpolatedUnivariateSpline(taurange,purityIg)
                interPurIgvals=interPurIgauss(increased_taurange)
                interPurIgaussvals=interPurIgauss(increased_taurange)
                maxIg=numpy.max(interPurIgvals)
                maxIgidx=numpy.argmax(interPurIgvals)

if includeintensity:
    print("maximum:\n\t\t\tpurity: {0:.4f}\tpurity(intensity): {1:.4f}".format(maxA,maxI))
else:
    print("maximum:\n\t\t\tpurity: {0:.4f}".format(maxA))

if plot:
        print('plotting...')
        fig = plt.figure(figsize=(8.75*1.2, 5*1.2))
        xax1=fig.add_subplot(111)
        xax2=xax1.twiny()
        AnnotateString=''
        #plt.rc('text', usetex=True)
        AnnotateString=AnnotateString+r'Maximum purity (sech): {0:.3} at $\tau(sech)={1:.3}$ps, $\tau(AC)={2:.3}$ps'.format(maxA,increased_taurange[maxAidx]*taucfsech*10**12,increased_taurange[maxAidx]*10**12)+'\n'
        if usefilter:
            AnnotateString=AnnotateString+FilterString
        #AnnotateString=AnnotateString+'Maximum purity(intensity): {0:.3} at {1:.3}ps\n'.format(maxI,increased_taurange[maxIidx]*taucfsech*10**12)
        #AnnotateString=AnnotateString+'Maximum purity: {0:.3} at {1:.3}ps\n'.format(maxA,increased_taurange[maxAidx]*taucfsech*10**12)
        #AnnotateString=AnnotateString+'Maximum purity: {0:.3} at {1:.3}ps\n'.format(maxA,increased_taurange[maxAidx]*10**12)
        #AnnotateString=AnnotateString+'Maximum purity(intensity): {0:.3} at {1:.3}ps\n'.format(maxI,increased_taurange[maxIidx]*taucfsech*10**12)
        xax1.set_xlabel(r'Pulsewidth (sech$^2$ pulse)[ps]')
        xax2.set_xlabel(r'Pulsewidth (gaussian pulse)[ps]')
        #xax2.set_xlabel(r'Arbitrary scale')
        xax1.set_ylabel('Purity')
        xax1.plot(taurange*taucfsech*10**12,purityA,label='Purity(sech)',color='orange')
        if includeintensity:
            xax1.plot(taurange*taucfsech*10**12,purityI,label='Purity(intensity,sech)')
        if includegauss:
                xax2.plot(taurange*taucfgauss*10**12,purityAg,label='Purity(gauss)',color='blue')
                AnnotateString=AnnotateString+r'Maximum purity (gauss): {0:.3} at $\tau(gauss)={1:.3}$ps, $\tau(AC)={2:.3}$ps'.format(maxAg,increased_taurange[maxAgidx]*taucfgauss*10**12,increased_taurange[maxAgidx]*10**12)+'\n'
                if includeintensity:
                    xax1.plot(taurange*taucfgauss*10**12,purityIg,label='Purity(intensity,gauss)')
        
        AnnotateString = AnnotateString + 'Pump wavelength: {0:.2f}nm\nPoling period: {1:.4f}µm\nTemperature: {2:.1f}°C\nCrystal length: {3:.1f}mm'.format(pwl*10**9,PP*10**6,T,L*10**3)
        fig.suptitle(plottitle,fontsize=16)
        xax1.grid()
        #xax1.legend()
        handles1, labels1 = xax1.get_legend_handles_labels() 
        handles2, labels2 = xax2.get_legend_handles_labels() 
        xax2.legend(handles1+handles2, labels1+labels2) 
        plt.annotate(AnnotateString,xy=(0.2,0.01),xycoords='axes fraction')
        if saveplot:
                print('saving plot...')
                save_figure_fname = os.path.join(os.getcwd(), plotsavename)
                plt.savefig(save_figure_fname, format="pdf",transparent=True, bbox_inches='tight', pad_inches=0)
        plt.show()
