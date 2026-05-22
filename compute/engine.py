"""
ComputeEngine — pure-numerics module with no GUI dependencies.

Every computation the GUI can trigger becomes a method here.
Each method takes a plain params dict and returns a plain result dict.
Refractive indices are resolved by name (material, polarisation, paper).
"""

import numpy as np
import scipy
import scipy.linalg
import scipy.interpolate
import scipy.optimize

from RefractiveIndex import RefractiveIndex
from PMC import PMC
from JSI import JSI
from Filters import Filters
from Constants import Constants


class ComputeEngine:
    """Stateless compute engine. Each method receives all parameters it needs."""

    # ------------------------------------------------------------------ #
    #  Helpers                                                             #
    # ------------------------------------------------------------------ #

    def _resolve_refidx(self, params: dict):
        """Resolve refractive-index functions from name parameters."""
        material = params['material']
        ri = RefractiveIndex()
        nxfunc = ri.getSingleIDX(material, "X", params['nx_paper'])
        nyfunc = ri.getSingleIDX(material, "Y", params['ny_paper'])
        nzfunc = ri.getSingleIDX(material, "Z", params['nz_paper'])
        return [nxfunc, nyfunc, nzfunc]

    def _resolve_filters(self, params: dict):
        """Resolve filter functions from parameters."""
        f = Filters()
        ffs = f.getFilterFunction(
            params.get('signal_filter_type', 'None'),
            params.get('signal_filter_center_wl', 1550e-9),
            params.get('signal_filter_fwhm', 3e-9))
        ffi = f.getFilterFunction(
            params.get('idler_filter_type', 'None'),
            params.get('idler_filter_center_wl', 1550e-9),
            params.get('idler_filter_fwhm', 3e-9))
        return [ffs, ffi]

    def _get_si_wavelengths(self, pwl, PP, T, refidxfunc, m):
        """Get signal/idler centre wavelengths at a single temperature."""
        Tvec = np.arange(T, T + 1, 2)
        [ls, li, _] = PMC().getSI_wl_varT(pwl, PP, Tvec, refidxfunc, m)
        return ls, li

    def _configure_jsi(self, jsi, params):
        """Apply focusing/fibre-coupling settings to a JSI instance."""
        jsi.focusing_enable = params.get('focusing_enable', False)
        jsi.fibre_coupling_enable = params.get('fibre_coupling_enable', False)
        jsi.focallength_pump = params.get('focallength_pump', 10e-3)
        jsi.focallength_signal = params.get('focallength_signal', 10e-3)
        jsi.focallength_idler = params.get('focallength_idler', 10e-3)
        jsi.Beamdiameter_pump = params.get('beamdiameter_pump', 1e-3)
        jsi.Beamdiameter_signal = params.get('beamdiameter_signal', 1e-3)
        jsi.Beamdiameter_idler = params.get('beamdiameter_idler', 1e-3)

    # ------------------------------------------------------------------ #
    #  Metadata                                                            #
    # ------------------------------------------------------------------ #

    def get_available_materials(self) -> list:
        return RefractiveIndex().materialList

    def get_available_refractive_indices(self, material: str) -> list:
        return RefractiveIndex().getAvailableRefractiveIndices(material)

    # ------------------------------------------------------------------ #
    #  Compute methods                                                     #
    # ------------------------------------------------------------------ #

    def compute_refractive_index_vs_T(self, params: dict) -> dict:
        T_min = params['T_min']
        T_max = params['T_max']
        wl = params['wavelength']
        selected = params['selected_indices']  # list of (material, pol, paper)

        plotrange = np.linspace(T_min, T_max, 1000)
        results = {}
        for material, pol, paper in selected:
            func = RefractiveIndex().getSingleIDX(material, pol, paper)
            key = f"{material}:{pol}:{paper}"
            results[key] = func(wl, plotrange)

        return {'T_range': plotrange, 'n_values': results}

    def compute_refractive_index_vs_wl(self, params: dict) -> dict:
        wl_min = params['wl_min']
        wl_max = params['wl_max']
        T = params['temperature']
        selected = params['selected_indices']

        plotrange = np.linspace(wl_min, wl_max, 1000)
        results = {}
        for material, pol, paper in selected:
            func = RefractiveIndex().getSingleIDX(material, pol, paper)
            key = f"{material}:{pol}:{paper}"
            results[key] = func(plotrange, T)

        return {'wl_range': plotrange, 'n_values': results}

    def compute_pmc_vs_T(self, params: dict) -> dict:
        pwl = params['pump_wl']
        PP = params['poling_period']
        T_min = params['T_min']
        T_max = params['T_max']
        m = params['qpm_order']
        refidxfunc = self._resolve_refidx(params)

        plotrange = np.arange(T_min, T_max, (T_max - T_min) / 250)
        [siwl, idwl, Tcp] = PMC().getSI_wl_varT(pwl, PP, plotrange, refidxfunc, m)

        return {
            'T_range': plotrange,
            'signal_wl': siwl,
            'idler_wl': idwl,
            'Tcp': Tcp,
        }

    def compute_pmc_vs_PP(self, params: dict) -> dict:
        pwl = params['pump_wl']
        PP_min = params['PP_min']
        PP_max = params['PP_max']
        T = params['temperature']
        m = params['qpm_order']
        refidxfunc = self._resolve_refidx(params)

        plotrange = np.arange(PP_min, PP_max, (PP_max - PP_min) / 250)
        [siwl, idwl, PPcp] = PMC().getSI_wl_varPP(pwl, plotrange, T, refidxfunc, m)

        return {
            'PP_range': plotrange,
            'signal_wl': siwl,
            'idler_wl': idwl,
            'PPcp': PPcp,
        }

    def compute_jsi(self, params: dict) -> dict:
        pwl = params['pump_wl']
        PP = params['poling_period']
        L = params['crystal_length']
        T = params['temperature']
        m = params['qpm_order']
        tau = params['pulsewidth']
        pumpcwbw = params.get('pump_cw_bw', 0)
        wlrange = params['wl_range']
        numpts = params['resolution']
        pumpshape = params['pump_shape']
        plot_jsi = params.get('plot_jsi', True)
        refidxfunc = self._resolve_refidx(params)
        spectralfilters = self._resolve_filters(params)

        ls, li = self._get_si_wavelengths(pwl, PP, T, refidxfunc, m)
        signalrange = np.linspace(ls - wlrange / 2, ls + wlrange / 2, numpts)
        idlerrange = np.linspace(li - wlrange / 2, li + wlrange / 2, numpts)

        jsi = JSI()
        self._configure_jsi(jsi, params)

        [PE, PM, JS] = jsi.getplots(
            pwl, signalrange, idlerrange, tau, T, PP, L, refidxfunc,
            m, spectralfilters, plot_jsi, pumpshape, pumpcwbw,
            params.get('focusing_enable', False),
            params.get('fibre_coupling_enable', False),
            params.get('focallength_pump', 10e-3),
            params.get('focallength_signal', 10e-3),
            params.get('focallength_idler', 10e-3),
            params.get('beamdiameter_pump', 1e-3),
            params.get('beamdiameter_signal', 1e-3),
            params.get('beamdiameter_idler', 1e-3))

        return {
            'PE': PE, 'PM': PM, 'JS': JS,
            'signal_range': signalrange,
            'idler_range': idlerrange,
            'pump_shape': pumpshape,
            'plot_jsi': plot_jsi,
        }

    def compute_purity_vs_tau(self, params: dict) -> dict:
        pwl = params['pump_wl']
        PP = params['poling_period']
        L = params['crystal_length']
        T = params['temperature']
        m = params['qpm_order']
        tau_min = params['tau_min']
        tau_max = params['tau_max']
        wlrange = params['wl_range']
        wlpts = params['wl_resolution']
        taupts = params['tau_resolution']
        pumpshape = params['pump_shape']
        refidxfunc = self._resolve_refidx(params)
        spectralfilters = self._resolve_filters(params)

        taurange = np.arange(tau_min, tau_max, (tau_max - tau_min) / taupts)
        ls, li = self._get_si_wavelengths(pwl, PP, T, refidxfunc, m)
        signalrange = np.linspace(ls - wlrange / 2, ls + wlrange / 2, wlpts)
        idlerrange = np.linspace(li - wlrange / 2, li + wlrange / 2, wlpts)

        [purity, max_pur, max_tau] = JSI().getpurity_vsTau(
            pwl, signalrange, idlerrange, taurange, T,
            PP, L, refidxfunc, m, spectralfilters, pumpshape)

        return {
            'tau_range': taurange,
            'purity': purity,
            'max_purity': max_pur,
            'max_tau': max_tau,
        }

    def compute_purity_vs_L(self, params: dict) -> dict:
        pwl = params['pump_wl']
        PP = params['poling_period']
        L_min = params['L_min']
        L_max = params['L_max']
        T = params['temperature']
        m = params['qpm_order']
        tau = params['pulsewidth']
        pumpcwbw = params.get('pump_cw_bw', 0)
        wlrange = params['wl_range']
        wlpts = params['wl_resolution']
        pts = params['L_resolution']
        pumpshape = params['pump_shape']
        refidxfunc = self._resolve_refidx(params)
        spectralfilters = self._resolve_filters(params)

        Lrange = np.arange(L_min, L_max, (L_max - L_min) / pts)
        ls, li = self._get_si_wavelengths(pwl, PP, T, refidxfunc, m)
        signalrange = np.linspace(ls - wlrange / 2, ls + wlrange / 2, wlpts)
        idlerrange = np.linspace(li - wlrange / 2, li + wlrange / 2, wlpts)

        [purity, max_pur, max_L] = JSI().getpurity_vsL(
            pwl, signalrange, idlerrange, tau, T,
            PP, Lrange, refidxfunc, m, spectralfilters, pumpshape, pumpcwbw)

        return {
            'L_range': Lrange,
            'purity': purity,
            'max_purity': max_pur,
            'max_L': max_L,
        }

    def compute_purity_vs_L_and_tau(self, params: dict) -> dict:
        pwl = params['pump_wl']
        PP = params['poling_period']
        L_min = params['L_min']
        L_max = params['L_max']
        T = params['temperature']
        m = params['qpm_order']
        tau_min = params['tau_min']
        tau_max = params['tau_max']
        wlrange = params['wl_range']
        wlpts = params['wl_resolution']
        pts = params['resolution']
        pumpshape = params['pump_shape']
        refidxfunc = self._resolve_refidx(params)
        spectralfilters = self._resolve_filters(params)

        Lrange = np.arange(L_min, L_max, (L_max - L_min) / pts)
        Taurange = np.arange(tau_min, tau_max, (tau_max - tau_min) / pts)
        ls, li = self._get_si_wavelengths(pwl, PP, T, refidxfunc, m)
        signalrange = np.linspace(ls - wlrange / 2, ls + wlrange / 2, wlpts)
        idlerrange = np.linspace(li - wlrange / 2, li + wlrange / 2, wlpts)

        purity = JSI().getpurity_vsLandTau(
            pwl, signalrange, idlerrange, Taurange, T,
            PP, Lrange, refidxfunc, m, spectralfilters, pumpshape)

        return {
            'L_range': Lrange,
            'tau_range': Taurange,
            'purity': purity,
        }

    def compute_effective_PP(self, params: dict) -> dict:
        m = params['qpm_order']
        Tcp = params['temperature']
        PP_guess = params['PP_guess']
        pwl = params['pump_wl']
        refidxfunc = self._resolve_refidx(params)

        PP = JSI().GetEffectivePP(m, Tcp, PP_guess, pwl, refidxfunc)
        return {'effective_PP': PP}

    def compute_tcp_vs_PP(self, params: dict) -> dict:
        PP_min = params['PP_min']
        PP_max = params['PP_max']
        T = params['temperature']
        pwl = params['pump_wl']
        m = params['qpm_order']
        refidxfunc = self._resolve_refidx(params)

        PPrange = np.linspace(PP_min, PP_max, 100)
        Tcp = JSI().getTcpVsPP(PPrange, T, pwl, refidxfunc, m)

        return {'PP_range': PPrange, 'Tcp': Tcp}

    def compute_tcp_vs_lp(self, params: dict) -> dict:
        pwl_min = params['pwl_min']
        pwl_max = params['pwl_max']
        T = params['temperature']
        PP = params['poling_period']
        m = params['qpm_order']
        refidxfunc = self._resolve_refidx(params)

        pwlrange = np.linspace(pwl_min, pwl_max, 100)
        Tcp = JSI().getTcpVslp(pwlrange, T, PP, refidxfunc, m)

        return {'pwl_range': pwlrange, 'Tcp': Tcp}

    def compute_hom_interference(self, params: dict) -> dict:
        pwl = params['pump_wl']
        T = params['temperature']
        PP = params['poling_period']
        m = params['qpm_order']
        tau = params['pulsewidth']
        pumpcwbw = params.get('pump_cw_bw', 0)
        cl = params['crystal_length']
        pumpshape = params['pump_shape']
        delay_range_width = params['delay_range']
        resolution = params['hom_resolution']
        homphase = params.get('hom_phase', 0)
        jsi_resolution = params['jsi_resolution']
        jsi_wlrange = params['jsi_wl_range']
        refidxfunc = self._resolve_refidx(params)
        spectralfilters = self._resolve_filters(params)

        delayrange = np.linspace(-delay_range_width / 2,
                                  delay_range_width / 2, resolution)

        ls, li = self._get_si_wavelengths(pwl, PP, T, refidxfunc, m)
        signalrange = np.linspace(ls - jsi_wlrange / 2,
                                  ls + jsi_wlrange / 2, jsi_resolution)
        idlerrange = np.linspace(li - jsi_wlrange / 2,
                                 li + jsi_wlrange / 2, jsi_resolution)

        jsi = JSI()
        [CoincProb, vis, fwhm] = jsi.getHOMinterference(
            pwl, T, PP, m, tau, cl, signalrange, idlerrange,
            jsi_resolution, pumpshape, delayrange, homphase,
            refidxfunc, spectralfilters, pumpcwbw,
            params.get('focusing_enable', False),
            params.get('fibre_coupling_enable', False),
            params.get('focallength_pump', 10e-3),
            params.get('focallength_signal', 10e-3),
            params.get('focallength_idler', 10e-3),
            params.get('beamdiameter_pump', 1e-3),
            params.get('beamdiameter_signal', 1e-3),
            params.get('beamdiameter_idler', 1e-3))

        return {
            'delay_range': delayrange,
            'coinc_prob': CoincProb,
            'visibility': vis,
            'fwhm': fwhm,
        }

    def compute_hom_interference_T(self, params: dict) -> dict:
        pwl = params['pump_wl']
        T = params['temperature']
        PP = params['poling_period']
        m = params['qpm_order']
        tau = params['pulsewidth']
        pumpcwbw = params.get('pump_cw_bw', 0)
        cl = params['crystal_length']
        pumpshape = params['pump_shape']
        T_min = params['T_min']
        T_max = params['T_max']
        resolution = params['hom_resolution']
        homphase = params.get('hom_phase', 0)
        jsi_resolution = params['jsi_resolution']
        jsi_wlrange = params['jsi_wl_range']
        hom_temprange = params['hom_temp_range']
        refidxfunc = self._resolve_refidx(params)
        spectralfilters = self._resolve_filters(params)

        # Calculate crossing point temperature
        plotrange_tcp = np.arange(T_min, T_max, (T_max - T_min) / 250)
        [_, _, Tcp] = PMC().getSI_wl_varT(pwl, PP, plotrange_tcp, refidxfunc, m)

        temprange = np.linspace(Tcp - hom_temprange / 2,
                                Tcp + hom_temprange / 2, resolution)

        ls, li = self._get_si_wavelengths(pwl, PP, T, refidxfunc, m)
        signalrange = np.linspace(ls - jsi_wlrange / 2,
                                  ls + jsi_wlrange / 2, jsi_resolution)
        idlerrange = np.linspace(li - jsi_wlrange / 2,
                                 li + jsi_wlrange / 2, jsi_resolution)

        jsi = JSI()
        [CoincProb, vis, fwhm] = jsi.getHOMinterferenceT(
            pwl, PP, m, tau, cl, signalrange, idlerrange,
            jsi_resolution, pumpshape, temprange, homphase,
            refidxfunc, spectralfilters, pumpcwbw,
            params.get('focusing_enable', False),
            params.get('fibre_coupling_enable', False),
            params.get('focallength_pump', 10e-3),
            params.get('focallength_signal', 10e-3),
            params.get('focallength_idler', 10e-3),
            params.get('beamdiameter_pump', 1e-3),
            params.get('beamdiameter_signal', 1e-3),
            params.get('beamdiameter_idler', 1e-3))

        return {
            'temp_range': temprange,
            'coinc_prob': CoincProb,
            'visibility': vis,
            'fwhm': fwhm,
        }

    def compute_fwhm_vs_tau(self, params: dict) -> dict:
        pwl = params['pump_wl']
        T = params['temperature']
        PP = params['poling_period']
        m = params['qpm_order']
        cl = params['crystal_length']
        pumpshape = params['pump_shape']
        jsi_resolution = params['jsi_resolution']
        jsi_wlrange = params['jsi_wl_range']
        tau_min = params['tau_min']
        tau_max = params['tau_max']
        fwhmres = params['fwhm_resolution']
        decprec = params['fwhm_precision']
        usetaucf = params.get('use_taucf', False)
        refidxfunc = self._resolve_refidx(params)
        spectralfilters = self._resolve_filters(params)

        ls, li = self._get_si_wavelengths(pwl, PP, T, refidxfunc, m)
        signalrange = np.linspace(ls - jsi_wlrange / 2,
                                  ls + jsi_wlrange / 2, jsi_resolution)
        idlerrange = np.linspace(li - jsi_wlrange / 2,
                                 li + jsi_wlrange / 2, jsi_resolution)
        taurange = np.linspace(tau_min, tau_max, fwhmres)

        [sigfwhm, idfwhm] = JSI().getFWHMvstau(
            pwl, signalrange, idlerrange, T, PP, m, cl, taurange,
            refidxfunc, spectralfilters, jsi_resolution, pumpshape,
            decprec, usetaucf)

        return {
            'tau_range': taurange,
            'signal_fwhm': sigfwhm,
            'idler_fwhm': idfwhm,
        }

    def estimate_filter_losses(self, params: dict) -> dict:
        pwl = params['pump_wl']
        PP = params['poling_period']
        L = params['crystal_length']
        T = params['temperature']
        m = params['qpm_order']
        tau = params['pulsewidth']
        pumpcwbw = params.get('pump_cw_bw', 0)
        wlrange = params['wl_range']
        numpts = params['resolution']
        pumpshape = params['pump_shape']
        refidxfunc = self._resolve_refidx(params)
        spectralfilters = self._resolve_filters(params)

        # Reference: no filters
        ffiref = Filters().getFilterFunction('None', 1, 1)
        ffsref = Filters().getFilterFunction('None', 1, 1)
        spectralfiltersref = [ffsref, ffiref]

        ls, li = self._get_si_wavelengths(pwl, PP, T, refidxfunc, m)
        signalrange = np.linspace(ls - wlrange / 2, ls + wlrange / 2, numpts)
        idlerrange = np.linspace(li - wlrange / 2, li + wlrange / 2, numpts)

        jsi = JSI()
        [PE, PM, JS] = jsi.getplots(
            pwl, signalrange, idlerrange, tau, T, PP, L, refidxfunc,
            m, spectralfilters, True, pumpshape, pumpcwbw,
            False, False, 10e-3, 10e-3, 10e-3, 1e-3, 1e-3, 1e-3)
        [PEwoSL, PMwoSL, JSwoSL] = jsi.getplots(
            pwl, signalrange, idlerrange, tau, T, PP, L, refidxfunc,
            m, spectralfilters, True, pumpshape, pumpcwbw,
            False, False, 10e-3, 10e-3, 10e-3, 1e-3, 1e-3, 1e-3)
        [PEref, PMref, JSref] = jsi.getplots(
            pwl, signalrange, idlerrange, tau, T, PP, L, refidxfunc,
            m, spectralfiltersref, True, pumpshape, pumpcwbw,
            False, False, 10e-3, 10e-3, 10e-3, 1e-3, 1e-3, 1e-3)
        [PEwoSLref, PMwoSLref, JSwoSLref] = jsi.getplots(
            pwl, signalrange, idlerrange, tau, T, PP, L, refidxfunc,
            m, spectralfiltersref, True, pumpshape, pumpcwbw,
            False, False, 10e-3, 10e-3, 10e-3, 1e-3, 1e-3, 1e-3)

        # Find sidelobe boundaries (same algorithm as GUI.py)
        wlrangelen = len(signalrange)
        JSdiagonal = []
        for i in range(1, wlrangelen):
            JSdiagonal.append(JSwoSLref[i, wlrangelen - i])
        JSdiaginterp = scipy.interpolate.interp1d(
            np.linspace(signalrange[1], signalrange[-1], numpts - 1).flatten(),
            JSdiagonal, kind='cubic', bounds_error=False)

        wlrangehalf = int(len(signalrange) / 2)
        luvalold = JSdiaginterp(signalrange[wlrangehalf - 1])
        rlvalold = JSdiaginterp(signalrange[wlrangehalf + 1])
        foundlumin = False
        foundrlmin = False
        luminidx = wlrangehalf
        rlminidx = wlrangehalf
        for i in range(2, wlrangehalf):
            luval = JSdiaginterp(signalrange[wlrangehalf - i])
            rlval = JSdiaginterp(signalrange[wlrangehalf + i])
            if not foundlumin:
                if luval < luvalold:
                    luvalold = luval
                    luminidx = wlrangehalf - i
                else:
                    foundlumin = True
            if not foundrlmin:
                if rlval < rlvalold:
                    rlvalold = rlval
                    rlminidx = wlrangehalf + i
                else:
                    foundrlmin = True

        for i in range(0, len(signalrange)):
            for j in range(0, len(signalrange)):
                if i < (2 * luminidx - wlrangelen + j):
                    JSwoSLref[i, j] = 0
                    JSwoSL[i, j] = 0
                elif i > (2 * rlminidx - wlrangelen + j):
                    JSwoSLref[i, j] = 0
                    JSwoSL[i, j] = 0

        mag_wf = np.sum(JS)
        mag_wof = np.sum(JSref)
        mag_wf_wosl = np.sum(JSwoSL)
        mag_wof_wosl = np.sum(JSwoSLref)

        filterlosses = (1 - mag_wf / mag_wof) if mag_wof != 0 else 0
        nonsidelobefilterlosses = (1 - mag_wf_wosl / mag_wof_wosl) if mag_wof_wosl != 0 else 0
        sidelobelosses = (1 - mag_wof_wosl / mag_wof) if mag_wof != 0 else 0
        filteredsidelobelosses = (1 - mag_wf_wosl / mag_wf) if mag_wf != 0 else 0

        return {
            'signal_range': signalrange,
            'idler_range': idlerrange,
            'JS': JS, 'JSref': JSref,
            'JSwoSL': JSwoSL, 'JSwoSLref': JSwoSLref,
            'filter_losses': filterlosses,
            'nonsidelobe_filter_losses': nonsidelobefilterlosses,
            'sidelobe_losses': sidelobelosses,
            'filtered_sidelobe_losses': filteredsidelobelosses,
        }
