#!/usr/bin/env python3

import cProfile
import numba
import os

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

        # Spatial parameters for tight focusing and fiber coupling
        self.fibre_coupling_enable = False
        self.focusing_enable = False
        
        # Collimated beam diameters
        self.Beamdiameter_pump = 1.0 * 10 ** (-3)
        self.Beamdiameter_signal = 1.0 * 10 ** (-3)
        self.Beamdiameter_idler = 1.0 * 10 ** (-3)
        
        # Focal lengths
        self.focallength_pump = 10.0 * 10 ** (-3)
        self.focallength_signal = 10.0 * 10 ** (-3)
        self.focallength_idler = 10.0 * 10 ** (-3)

        # Spatial walk-off
        self.walkoff_enable = True

        # Number of z-samples used by batched spatial overlap integration.
        # Use an odd value for Simpson integration.
        self.spatial_z_points = 129

    def calculate_focused_waists(self, lp, ls, li):
        """
        Calculate the focused beam waists at the crystal.
        w_0 = \lambda f / (\pi w_collimated)
        spot size = spot diameter = 2x waist
        """
        w_col_p = self.Beamdiameter_pump / 2.0
        w_col_s = self.Beamdiameter_signal / 2.0
        w_col_i = self.Beamdiameter_idler / 2.0

        w0_p = (lp * self.focallength_pump) / (np.pi * w_col_p)
        w0_s = (ls * self.focallength_signal) / (np.pi * w_col_s)
        w0_i = (li * self.focallength_idler) / (np.pi * w_col_i)

        return w0_p, w0_s, w0_i

    def walkoff_angle(self, wavelength, T):
        """
        Calculate spatial walk-off angle for extraordinary (z-polarized) ray
        in KTP propagating along x-axis.

        For a biaxial crystal with propagation along x and z-polarization,
        the walk-off angle is approximately:
        rho ≈ (1/2) * (1/n_y^2 - 1/n_z^2) * n_z

        Returns angle in radians.
        """
        ny = self.ny(wavelength, T)
        nz = self.nz(wavelength, T)
        # Walk-off for e-ray in biaxial crystal (propagation along x, z-polarized)
        rho = 0.5 * (1.0 / ny**2 - 1.0 / nz**2) * nz
        return rho

    
    #def spatial_overlap(self, dk, cl, lp, ls, li):
    #    return self.spatial_overlap_numba(dk, cl, lp, ls, li, self.ny, self.nz, self.T, self.walkoff_enable, Beamdiameter_pump=self.Beamdiameter_pump, Beamdiameter_signal=self.Beamdiameter_signal, Beamdiameter_idler=self.Beamdiameter_idler, focallength_pump=self.focallength_pump, focallength_signal=self.focallength_signal, focallength_idler=self.focallength_idler)
    #@staticmethod
    #@numba.njit
    #def spatial_overlap_numba(dk, cl, lp, ls, li, ny, nz, T, walkoff_enable, Beamdiameter_pump, Beamdiameter_signal, Beamdiameter_idler, focallength_pump, focallength_signal, focallength_idler):
    #    """
    #    Spatial overlap integral for the biphoton amplitude including
    #    tight focusing, fiber coupling, and spatial walk-off.
    #    Performs rigorous numerical integration over the crystal length.
    #    """
    #    w_col_p = Beamdiameter_pump / 2.0
    #    w_col_s = Beamdiameter_signal / 2.0
    #    w_col_i = Beamdiameter_idler / 2.0#

  #      w0_p = (lp * focallength_pump) / (np.pi * w_col_p)
    #    w0_s = (ls * focallength_signal) / (np.pi * w_col_s)
    #    w0_i = (li * focallength_idler) / (np.pi * w_col_i)
  #
  #      # Rayleigh ranges including refractive index
  #      zr_p = ny(lp, T) * np.pi * w0_p**2 / lp
  #      zr_s = ny(ls, T) * np.pi * w0_s**2 / ls
  #      zr_i = nz(li, T) * np.pi * w0_i**2 / li
#
#        # Walk-off angle for idler (extraordinary beam)
#        if walkoff_enable:
#            rho_i = 0.5 * (1.0 / ny(li, T)**2 - 1.0 / nz(li, T)**2) * nz(li, T)
#        else:
#            rho_i = 0.0
#
#        def _spatial_integrand( z, dk, w0_p, w0_s, w0_i, zr_p, zr_s, zr_i, rho_i=0.0):
#
#            # Gaussian complex beam parameters
#            q_p = 1.0 + 1j * z / zr_p
#            q_s_conj = 1.0 - 1j * z / zr_s
#            q_i_conj = 1.0 - 1j * z / zr_i
#
#            # Inverse beam size parameters
#            a_p = 1.0 / (w0_p**2 * q_p)
#            a_s = 1.0 / (w0_s**2 * q_s_conj)
#            a_i = 1.0 / (w0_i**2 * q_i_conj)
#
#            # Sum of non-displaced beams (pump + signal)
#            A_prime = a_p + a_s
#            # Total including idler
#            A = A_prime + a_i
#
#            # Transverse displacement of idler beam center due to walk-off
#            d = rho_i * z
#
#            # Overlap integral with displaced Gaussian
#            # For three Gaussians where idler is displaced by distance d:
#            # int exp(-r^2*A_prime - |r-d|^2*a_i) d^2r = (pi/A) * exp(-d^2 * a_i * A_prime / A)
#            displacement_factor = np.exp(-d**2 * a_i * A_prime / A)
#
#            # Gouy phase, normalization, displacement, and longitudinal phase matching
#            amplitude = (1.0 / (q_p * q_s_conj * q_i_conj)) * (np.pi / A) * displacement_factor * np.exp(1j * dk * z)
#
#            return amplitude
#
#        def integrand_real(z):
#            return _spatial_integrand(z, dk, w0_p, w0_s, w0_i, zr_p, zr_s, zr_i, rho_i).real
#
#        def integrand_imag(z):
#            return _spatial_integrand(z, dk, w0_p, w0_s, w0_i, zr_p, zr_s, zr_i, rho_i).imag
#
#        # Integrate over crystal length from -cl/2 to cl/2
#        integral_real, _ = scipy.integrate.quad(integrand_real, -cl / 2, cl / 2)
#        integral_imag, _ = scipy.integrate.quad(integrand_imag, -cl / 2, cl / 2)
#
#        return integral_real + 1j * integral_imag
    
    def spatial_overlap(self, dk, cl, lp, ls, li, temp=None):
        """
        Spatial overlap integral for the biphoton amplitude including
        tight focusing, fiber coupling, and spatial walk-off.
        Performs rigorous numerical integration over the crystal length.
        """
        w0_p, w0_s, w0_i = self.calculate_focused_waists(lp, ls, li)#

        # Rayleigh ranges including refractive index
        t_eval = self.T if temp is None else temp
        zr_p = self.ny(lp, t_eval) * np.pi * w0_p**2 / lp
        zr_s = self.ny(ls, t_eval) * np.pi * w0_s**2 / ls
        zr_i = self.nz(li, t_eval) * np.pi * w0_i**2 / li

        # Walk-off angle for idler (extraordinary beam)
        rho_i = self.walkoff_angle(li, t_eval) if self.walkoff_enable else 0.0

        def integrand_real(z):
            return self._spatial_integrand(z, dk, w0_p, w0_s, w0_i, zr_p, zr_s, zr_i, rho_i).real

        def integrand_imag(z):
            return self._spatial_integrand(z, dk, w0_p, w0_s, w0_i, zr_p, zr_s, zr_i, rho_i).imag

        # Integrate over crystal length from -cl/2 to cl/2
        integral_real, _ = scipy.integrate.quad(integrand_real, -cl / 2, cl / 2)
        integral_imag, _ = scipy.integrate.quad(integrand_imag, -cl / 2, cl / 2)

        return integral_real + 1j * integral_imag

    def _spatial_overlap_grid(self, dk, cl, lp, ls, li, temp=None):
        """
        Batched spatial overlap integration on a shared z-grid.

        This avoids per-point scipy.quad calls and is substantially faster for
        evaluating large wavelength grids.
        """
        dk_arr, lp_arr, ls_arr, li_arr = np.broadcast_arrays(
            np.asarray(dk),
            np.asarray(lp),
            np.asarray(ls),
            np.asarray(li),
        )
        cl_val = float(np.asarray(cl).reshape(-1)[0])

        w0_p, w0_s, w0_i = self.calculate_focused_waists(lp_arr, ls_arr, li_arr)

        t_eval = self.T if temp is None else temp
        zr_p = self.ny(lp_arr, t_eval) * np.pi * w0_p ** 2 / lp_arr
        zr_s = self.ny(ls_arr, t_eval) * np.pi * w0_s ** 2 / ls_arr
        zr_i = self.nz(li_arr, t_eval) * np.pi * w0_i ** 2 / li_arr

        rho_i = self.walkoff_angle(li_arr, t_eval) if self.walkoff_enable else 0.0

        npts = int(self.spatial_z_points)
        if npts < 3:
            npts = 3
        if npts % 2 == 0:
            npts += 1

        z = np.linspace(-cl_val / 2, cl_val / 2, npts)
        dz = z[1] - z[0]

        # Memory-efficient Simpson integration: avoid storing a full (..., z) cube.
        weights = np.ones(npts)
        weights[1:-1:2] = 4.0
        weights[2:-1:2] = 2.0

        integral = np.zeros_like(dk_arr, dtype=np.complex128)

        for k, zk in enumerate(z):
            q_p = 1.0 + 1j * zk / zr_p
            q_s_conj = 1.0 - 1j * zk / zr_s
            q_i_conj = 1.0 - 1j * zk / zr_i

            a_p = 1.0 / (w0_p ** 2 * q_p)
            a_s = 1.0 / (w0_s ** 2 * q_s_conj)
            a_i = 1.0 / (w0_i ** 2 * q_i_conj)

            A_prime = a_p + a_s
            A = A_prime + a_i

            d = rho_i * zk
            displacement_factor = np.exp(-d ** 2 * a_i * A_prime / A)

            amplitude = (
                (1.0 / (q_p * q_s_conj * q_i_conj))
                * (np.pi / A)
                * displacement_factor
                * np.exp(1j * dk_arr * zk)
            )
            integral += weights[k] * amplitude

        return integral * (dz / 3.0)

    def _spatial_integrand(self, z, dk, w0_p, w0_s, w0_i, zr_p, zr_s, zr_i, rho_i=0.0):
        """
        Integrand for spatial overlap calculation with walk-off.

        Args:
            z: Longitudinal position in crystal (m)
            dk: Phase mismatch (1/m)
            w0_p, w0_s, w0_i: Beam waists (m)
            zr_p, zr_s, zr_i: Rayleigh ranges (m)
            rho_i: Walk-off angle of idler beam (radians)

        Returns:
            Complex amplitude contribution at position z
        """
        # Gaussian complex beam parameters
        q_p = 1.0 + 1j * z / zr_p
        q_s_conj = 1.0 - 1j * z / zr_s
        q_i_conj = 1.0 - 1j * z / zr_i

        # Inverse beam size parameters
        a_p = 1.0 / (w0_p**2 * q_p)
        a_s = 1.0 / (w0_s**2 * q_s_conj)
        a_i = 1.0 / (w0_i**2 * q_i_conj)

        # Sum of non-displaced beams (pump + signal)
        A_prime = a_p + a_s
        # Total including idler
        A = A_prime + a_i

        # Transverse displacement of idler beam center due to walk-off
        d = rho_i * z

        # Overlap integral with displaced Gaussian
        # For three Gaussians where idler is displaced by distance d:
        # int exp(-r^2*A_prime - |r-d|^2*a_i) d^2r = (pi/A) * exp(-d^2 * a_i * A_prime / A)
        displacement_factor = np.exp(-d**2 * a_i * A_prime / A)

        # Gouy phase, normalization, displacement, and longitudinal phase matching
        amplitude = (1.0 / (q_p * q_s_conj * q_i_conj)) * (np.pi / A) * displacement_factor * np.exp(1j * dk * z)

        return amplitude

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
        return self.lamdap_numba(ls, li)
    @staticmethod
    @numba.njit
    def lamdap_numba(ls, li):
        return 1 / (1 / ls + 1 / li)

    # phase matching amplitude
    def PMA(self, dk, cl, lp=None, ls=None, li=None, temp=None):
        if (self.focusing_enable or self.fibre_coupling_enable) and lp is not None and ls is not None and li is not None:
            return self._spatial_overlap_grid(dk, cl, lp, ls, li, temp=temp)
        # note: np.sinc(x) evaluates to Sin(pi*x)/(pi*x)
        # this really produced some headache.
        arg = cl * dk / 2
        return cl * np.sinc(arg/np.pi ) * np.exp( 1j * arg )
        #return np.sinc(dk * cl / (2 * Constants().pi)) #* cl * np.exp(1j*dk*cl/2)

    def PMAgauss(self, dk, cl, lp=None, ls=None, li=None, temp=None):
        return self.PMA(dk, cl, lp, ls, li, temp=temp)

    def PMAcwgauss(self, dk, cl, lp=None, ls=None, li=None, temp=None):
        return self.PMA(dk, cl, lp, ls, li, temp=temp)

    def PMAsech(self, dk, cl, lp=None, ls=None, li=None, temp=None):
        return self.PMA(dk, cl, lp, ls, li, temp=temp)

    def PMAsinc(self, dk, cl, lp=None, ls=None, li=None, temp=None):
        return self.PMA(dk, cl, lp, ls, li, temp=temp)

    # phase matching intensity
    def PMIsech(self, dk, cl, lp=None, ls=None, li=None):
        return np.square(np.abs(self.PMA(dk, cl, lp, ls, li)))

    def PMIgauss(self, dk, cl, lp=None, ls=None, li=None):
        return np.square(np.abs(self.PMA(dk, cl, lp, ls, li)))

    def PMIcwgauss(self, dk, cl, lp=None, ls=None, li=None):
        return np.square(np.abs(self.PMA(dk, cl, lp, ls, li)))

    def PMIsinc(self, dk, cl, lp=None, ls=None, li=None):
        return np.square(np.abs(self.PMA(dk, cl, lp, ls, li)))

    ###################
    ###  cw(gauss)  ###
    ###################
    def PEAcwgauss(self, lp, ls, li, bw):
        dnu = (1 / ls + 1 / li - 1 / lp) # /c
        snu = bw/(lp**2) # /c
        return (np.exp(- ( (dnu) / (2*snu)) ** 2))  # note: 2*pi*c/(2*sp)
        # return (np.exp(- (2*Constants().pi * Constants().c * (dl) / (2*2*Constants().pi * Constants().c*bw)) ** 2))  # note: 2*pi*c/(2*sp)

    # joint spectral amplitude for gaussian beam
    def JSAcwgauss(self, lp, ls, li, bw, t, pp, cl):
        sp = bw / (2 * np.sqrt(2 * np.log(2)))  # gaussian standard deviation from FWHM
        dk = self.deltak(self.lambdap(ls, li), ls, li, t, pp)
        jsa = self.PEAcwgauss(lp, ls, li, sp) * self.PMAcwgauss(dk, cl, lp, ls, li, temp=t)
        if self.useabs:
            return np.absolute(jsa)
        else:
            return jsa

    def PEIcwgauss(self, lp, ls, li, sp):
        return np.abs(self.PEAcwgauss(lp, ls, li, sp)) ** 2

    # joint spectral intensity for gaussian beam
    def JSIcwgauss(self, lp, ls, li, bw, t, pp, cl):
        return np.abs(self.JSAcwgauss(lp, ls, li, bw, t, pp, cl)) ** 2

    def PEAnPMAnJSAcwgauss(self, lp, ls, li, bw, t, pp, cl):
        sp = bw / (2 * np.sqrt(2 * np.log(2)))  # gaussian standard deviation from FWHM
        dk = self.deltak(self.lambdap(ls, li), ls, li, t, pp)
        pea = self.PEAcwgauss(lp, ls, li, sp)
        pma = self.PMAcwgauss(dk, cl, lp, ls, li, temp=t)
        jsa = pea * pma
        if self.useabs:
            return [np.absolute(pea), np.absolute(pma), np.absolute(jsa)]
        else:
            return [pea, pma, jsa]

    def PEInPMInJSIcwgauss(self, lp, ls, li, bw, t, pp, cl):
        [pea, pma, jsa] = self.PEAnPMAnJSAcwgauss(lp, ls, li, bw, t, pp, cl)
        return [np.abs(pea) ** 2, np.abs(pma) ** 2, np.abs(jsa) ** 2]
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
        return np.abs(self.PEAgauss(lp, ls, li, sp)) ** 2

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
        jsa = self.PEAgauss(lp, ls, li, sp) * self.PMAgauss(dk, cl, lp, ls, li, temp=t)
        if self.useabs:
            return np.absolute(jsa)
        else:
            return jsa

    # joint spectral intensity for gaussian beam
    def JSIgauss(self, lp, ls, li, tauac, t, pp, cl):
        return np.abs(self.JSAgauss(lp, ls, li, tauac, t, pp, cl)) ** 2

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
        pma = self.PMAgauss(dk, cl, lp, ls, li, temp=t)
        jsa = pea * pma
        if self.useabs:
            return [np.absolute(pea), np.absolute(pma), np.absolute(jsa)]
        else:
            return [pea, pma, jsa]

    def PEInPMInJSIgauss(self, lp, ls, li, tauac, t, pp, cl):
        [pea, pma, jsa] = self.PEAnPMAnJSAgauss(lp, ls, li, tauac, t, pp, cl)
        return [np.abs(pea) ** 2, np.abs(pma) ** 2, np.abs(jsa) ** 2]

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
        return np.abs(self.PEAsech(lp, ls, li, B)) ** 2

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
        jsa = self.PEAsech(lp, ls, li, B) * self.PMAsech(dk, cl, lp, ls, li, temp=t)
        
        if self.useabs:
            return np.absolute(jsa)
        else:
            return jsa

    # joint spectral intensity for sech^2 beam
    def JSIsech(self, lp, ls, li, tauac, t, pp, cl):
        return np.abs(self.JSAsech(lp, ls, li, tauac, t, pp, cl)) ** 2

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
        pma = self.PMAsech(dk, cl, lp, ls, li, temp=t)
        jsa = pea * pma
        if self.useabs:
            return [np.absolute(pea), np.absolute(pma), np.absolute(jsa)]
        else:
            return [pea, pma, jsa]

    def PEInPMInJSIsech(self, lp, ls, li, tauac, t, pp, cl):
        [pea, pma, jsa] = self.PEAnPMAnJSAsech(lp, ls, li, tauac, t, pp, cl)
        return [np.abs(pea) ** 2, np.abs(pma) ** 2, np.abs(jsa) ** 2]
    
    ###################
    ###    sinc     ###
    ###################
    # pump envelope amplitude for sinc beam
    def PEAsinc(self, lp, ls, li, B):
        dl = 1 / ls + 1 / li - 1 / lp
        wfact = 2 * Constants().pi * Constants().c * (1 / ls + 1 / li - 1 / lp)
        #argument = wfact * B
        argument = wfact*B / Constants().pi #/pi to counteract np.sincs spurious pi: np.sinc(x) = sin(pi*x)/(pi*x)
        return (np.sinc(argument))
        
    # pump envelope intensity for sinc beam
    def PEIsinc(self, lp, ls, li, B):
        return np.abs(self.PEAsinc(lp, ls, li, B)) ** 2
        
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
        jsa = self.PEAsinc(lp, ls, li, B) * self.PMAsinc(dk, cl, lp, ls, li, temp=t)
        if self.useabs:
            return np.absolute(jsa)
        else:
            return jsa

    # joint spectral intensity for sinc beam
    def JSIsinc(self, lp, ls, li, tauac, t, pp, cl):
        return np.abs(self.JSAsinc(lp, ls, li, tauac, t, pp, cl)) ** 2
        
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
        pma = self.PMAsinc(dk, cl, lp, ls, li, temp=t)
        jsa = pea * pma
        if self.useabs:
            return [np.absolute(pea), np.absolute(pma), np.absolute(jsa)]
        else:
            return [pea, pma, jsa]

    def PEInPMInJSIsinc(self, lp, ls, li, tauac, t, pp, cl):
        [pea, pma, jsa] = self.PEAnPMAnJSAsinc(lp, ls, li, tauac, t, pp, cl)
        return [np.abs(pea) ** 2, np.abs(pma) ** 2, np.abs(jsa) ** 2]

    def GetEffectivePP(self, m, Tcp, PPguess, lp, refidxfunc):
        self.nx = refidxfunc[0]
        self.ny = refidxfunc[1]
        self.nz = refidxfunc[2]
        return scipy.optimize.newton_krylov(self.wlgaponlyPP(lp, Tcp), PPguess, f_tol=1e-14)

    def getplots(self,pumpwl,signalrange,idlerrange,tau,temp,polingp,crystallength, refidxfunc,qpmorder,filterfuncs,plotJSI,pumpshape,pumpcwbw, focusing_enable, fibre_coupling_enable, focallength_pump, focallength_signal, focallength_idler, beamdiameter_pump, beamdiameter_signal, beamdiameter_idler):
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
        self.nx, self.ny, self.nz = refidxfunc

        self.pwl = pumpwl
        self.tau = tau
        self.pumpcwbw = pumpcwbw
        self.m = qpmorder
        self.T = temp
        self.PP = polingp * self.thermexpfactor(self.T)
        self.L = crystallength * self.thermexpfactor(self.T)

        self.pumpshape = pumpshape

        self.focusing_enable = focusing_enable
        self.fibre_coupling_enable = fibre_coupling_enable
        self.focallength_pump = focallength_pump
        self.focallength_signal = focallength_signal
        self.focallength_idler = focallength_idler
        self.Beamdiameter_pump = beamdiameter_pump
        self.Beamdiameter_signal = beamdiameter_signal
        self.Beamdiameter_idler = beamdiameter_idler

        self.calcJSI = plotJSI
        self.calcJSA = not plotJSI
            
        self.calcGaussian = False
        self.calcSech = False
        self.calcSinc = False
        self.calcCWGauss = False
        if self.pumpshape.casefold() =='gaussian':
            self.calcGaussian = True
            trifunc = self.PEAnPMAnJSAgauss if self.calcJSA else self.PEInPMInJSIgauss
            tau = self.tau
        elif self.pumpshape.casefold() == 'sinc':
            self.calcSinc = True
            trifunc = self.PEAnPMAnJSAsinc if self.calcJSA else self.PEInPMInJSIsinc
            tau = self.tau
        elif self.pumpshape.casefold() == 'sech^2':
            self.calcSech = True
            trifunc = self.PEAnPMAnJSAsech if self.calcJSA else self.PEInPMInJSIsech
            tau = self.tau
        elif self.pumpshape.casefold() == 'cw':
            self.calcCWGauss = True
            trifunc = self.PEAnPMAnJSAcwgauss if self.calcJSA else self.PEInPMInJSIcwgauss
            tau = self.pumpcwbw
        else:
            print('Error: No valid pump shape specified. Please choose between "gaussian", "sech^2", "sinc" and "cw".')

        X, Y = np.meshgrid(self.sigrange, self.idrange)

        [PE, PM, JS] = trifunc(self.pwl, X, Y, tau, self.T, self.PP, self.L)
        
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

    def getpurity_vsL(self,pumpwl,signalrange,idlerrange,tau,temp,polingp,crystallengthrange,refidxfunc,qpmorder,filterfuncs,pumpshape,pumpcwbw):
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
        self.pumpcwbw = pumpcwbw
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
            elif self.calcGaussian:
                JSA = self.JSAcwgauss(self.pwl, X, Y, self.pumpcwbw, self.T, self.PP, self.L)
                
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

    #by numerical integration
    def getHOMinterference(self, pwl, temp, polingp, qpmorder, tau, cl, signalrange, idlerrange,JSIresolution, pumpshape, delayrange, homphase, refidxfunc, filterfuncs, pumpcwbw, focusing_enable, fibre_coupling_enable, focallength_pump, focallength_signal, focallength_idler, beamdiameter_pump, beamdiameter_signal, beamdiameter_idler):
        t0=datetime.now()
        self.m = qpmorder
        [self.nx, self.ny, self.nz] = refidxfunc
        #
        # https://arxiv.org/pdf/1211.0120.pdf (On the Purity and Indistinguishability of Down-Converted Photons. Osorio, Sangouard, thew 2012)
        # Ansari, 2013 msc thesis
        #

        self.focusing_enable = focusing_enable
        self.fibre_coupling_enable = fibre_coupling_enable
        self.focallength_pump = focallength_pump
        self.focallength_signal = focallength_signal
        self.focallength_idler = focallength_idler
        self.Beamdiameter_pump = beamdiameter_pump
        self.Beamdiameter_signal = beamdiameter_signal
        self.Beamdiameter_idler = beamdiameter_idler

        

        X, Y = np.meshgrid(signalrange, idlerrange)

        self.calcGaussian, self.calcSech, self.calcSinc, self.calcCW = False, False, False, False
        if pumpshape.casefold() =='gaussian':
            self.calcGaussian = True
            jsafunc = self.JSAgauss
            peafunc = self.PEAgauss
            pmafunc = self.PMAgauss
            tau=tau
        elif pumpshape.casefold() =='sech^2':
            self.calcSech = True
            jsafunc = self.JSAsech
            peafunc = self.PEAsech
            pmafunc = self.PMAsech
            tau=tau
        elif pumpshape.casefold() =='sinc':
            self.calcSinc = True
            jsafunc = self.JSAsinc
            peafunc = self.PEAsinc
            pmafunc = self.PMAsinc
            tau=tau
        elif pumpshape.casefold() == 'cw':
            self.calcCW = True
            jsafunc = self.JSAcwgauss
            peafunc = self.PEAcwgauss
            pmafunc = self.PMAcwgauss
            tau=pumpcwbw

        cl = self.thermexpfactor(temp)*cl
        polingp = self.thermexpfactor(temp)*polingp
        self.PP = polingp * self.thermexpfactor(temp)
        self.L = cl * self.thermexpfactor(temp)


        jsa1 = jsafunc(pwl, X, Y, tau, temp, polingp, cl)
        jsa2 = jsafunc(pwl, Y, X, tau, temp, polingp, cl)
        jsa2c = np.conjugate(jsa2)

        norm = np.sum(scipy.integrate.simpson( ( np.abs(jsa1)**2)))
        def homf(i):
            return  scipy.integrate.simpson(scipy.integrate.simpson( np.abs(jsa1)**2 \
                   - jsa1 * jsa2c * np.exp(-1j*2*np.pi*Constants().c*(1/X-1/Y)*delayrange[i])   ))
        with concurrent.futures.ThreadPoolExecutor(max_workers=len(os.sched_getaffinity(0))) as ex:
            futures = [ex.submit(homf, i) for i in range(0,len(delayrange))]
            HOMI = np.array([f.result() for f in futures])


        HOMI = np.real(HOMI/norm)

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

    def getHOMinterferenceT(self, pwl, polingp, qpmorder, tau, cl, signalrange, idlerrange, JSIresolution, pumpshape, temprange, homphase, refidxfunc, filterfuncs, pumpcwbw, focusing_enable, fibre_coupling_enable, focallength_pump, focallength_signal, focallength_idler, beamdiameter_pump, beamdiameter_signal, beamdiameter_idler):
        t0=datetime.now()
        self.focusing_enable = focusing_enable
        self.fibre_coupling_enable = fibre_coupling_enable
        self.focallength_pump = focallength_pump
        self.focallength_signal = focallength_signal
        self.focallength_idler = focallength_idler
        self.Beamdiameter_pump = beamdiameter_pump
        self.Beamdiameter_signal = beamdiameter_signal
        self.Beamdiameter_idler = beamdiameter_idler
        
        self.m = qpmorder
        [self.nx, self.ny, self.nz] = refidxfunc
        #
        # https://arxiv.org/pdf/1211.0120.pdf (On the Purity and Indistinguishability of Down-Converted Photons. Osorio, Sangouard, thew 2012)
        # Ansari, 2013 msc thesis
        #
        X, Y = np.meshgrid(signalrange, idlerrange)

        self.calcGaussian, self.calcSech, self.calcSinc, self.calcCW = False, False, False, False
        if pumpshape.casefold() =='gaussian':
            self.calcGaussian = True
            jsafunc = self.JSAgauss
            peafunc = self.PEAgauss
            pmafunc = self.PMAgauss
            tau=tau
        elif pumpshape.casefold() =='sech^2':
            self.calcSech = True
            jsafunc = self.JSAsech
            peafunc = self.PEAsech
            pmafunc = self.PMAsech
            tau=tau
        elif pumpshape.casefold() =='sinc':
            self.calcSinc = True
            jsafunc = self.JSAsinc
            peafunc = self.PEAsinc
            pmafunc = self.PMAsinc
            tau=tau
        elif pumpshape.casefold() == 'cw':
            self.calcCW = True
            jsafunc = self.JSAcwgauss
            peafunc = self.PEAcwgauss
            pmafunc = self.PMAcwgauss
            tau=pumpcwbw

        if (pwl>400*10**(-9)) and (pwl<410*10**(-9)):
            delay=(1.47+0.28)*10**(-12) #404.87 source
            print("warning: custom delay applied")
        elif (pwl>770*10**(-9)) and (pwl<780*10**(-9)):
            delay=(4.1+0.32)*10**(-12) #773 source
            print("warning: custom delay applied")
        else:
            print("warning: custom delay NOT applied")

        def homf(i):
            clt = self.thermexpfactor(temprange[i])*cl
            jsa1 = jsafunc(pwl, X, Y, tau, temprange[i], polingp, clt)
            jsa2c = np.conjugate(jsafunc(pwl, Y, X, tau, temprange[i], polingp, clt))
            norm = scipy.integrate.simpson(scipy.integrate.simpson(  np.abs(jsa1)**2  )  )

            return  scipy.integrate.simpson(scipy.integrate.simpson( np.abs(jsa1)**2 \
                                    - jsa1 * jsa2c * np.exp(-1j*2*np.pi*Constants().c*(1/X-1/Y)*delay)   )) / norm
        max_workers = len(os.sched_getaffinity(0))
        if self.focusing_enable or self.fibre_coupling_enable:
            max_workers = max_workers//2

        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
            futures = [ex.submit(homf, i) for i in range(0,len(temprange))]
            HOMI = np.array([f.result() for f in futures])

        # calc FWHM
        if 0:
            temprangeneg = temprange[:int(np.floor(len(temprange)/2))]
            temprangepos = temprange[int(np.floor(len(temprange)/2)):]
            HOMIneg = HOMI[:int(np.floor(len(HOMI)/2))]
            HOMIpos = HOMI[int(np.floor(len(HOMI)/2)):]
            visinterpolfneg = scipy.interpolate.interp1d(temprangeneg, HOMIneg-0.25, fill_value='extrapolate')
            visinterpolfpos = scipy.interpolate.interp1d(temprangepos, HOMIpos-0.25, fill_value='extrapolate')
            negrootstartest=temprangeneg[int(np.floor(len(temprangeneg)/2))]
            posrootstartest=temprangepos[int(np.floor(len(temprangepos)/2))]
            negroot=scipy.optimize.fsolve(visinterpolfneg, negrootstartest)
            posroot=scipy.optimize.fsolve(visinterpolfpos, posrootstartest)

            homfwhm=posroot[0]-negroot[0]

            t1=datetime.now()
            print('calculating HOM took', (t1-t0).total_seconds(), 's')
        else:
            homfwhm=0
            vis=0
        return [HOMI,vis,homfwhm]

    def getFWHMvstau(self, pwl, signalrange, idlerrange, temp, polingp, qpmorder, cl, taurange, refidxfunc, filterfuncs, JSIresolution, pumpshape, decprec, usetaucf):
        [self.nx, self.ny, self.nz] = refidxfunc
        self.usetaucf = usetaucf
        
        self.calcGaussian = False
        self.calcSech = False
        self.calcSinc = False
        self.calcCWGauss = False
        if pumpshape.casefold() =='gaussian':
            self.calcGaussian = True
            jsifunc = self.JSIgauss
        elif pumpshape.casefold() =='sech^2':
            self.calcSech = True
            jsifunc = self.JSIsech
        elif pumpshape.casefold() =='sinc':
            self.calcSinc = True
            jsifunc = self.JSIsinc
        elif pumpshape.casefold() =='cw':
            self.calcCWGauss = True
            jsifunc = self.JSIcwgauss
        else:
            print('ERROR: Unknown pump beamshape')

        fwhmsig=[]
        fwhmid=[]

        #TMP
        self.useFilter=False

        for h in range(0,len(taurange)):
            result =jsifunc(pwl, signalrange[:,None], idlerrange[None,:], taurange[h], temp, polingp*self.thermexpfactor(temp), cl)
            maximum = np.max(result)
            result = result/maximum-0.5

            print(f"{result=}")
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
            print(f"{idcs=}")
            print(f"{result=}")
            for i in range(0,len(idcs[0])):
                    hmptss.append(signalrange[idcs[0][i]]*10**9)
                    hmptsi.append(idlerrange[idcs[1][i]]*10**9)
            print(f"{hmptss=}")
            if len(hmptss) != 0:
                if len(hmptsi) !=0:
                    fwhmsig.append(max(hmptss)-min(hmptss))
                    fwhmid.append(max(hmptsi)-min(hmptsi))
        return [fwhmsig,fwhmid]
