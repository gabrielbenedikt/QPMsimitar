import pytest
import numpy as np
from JSI import JSI
from RefractiveIndex import RefractiveIndex
from Filters import Filters

def test_jsi_initialization():
    jsi = JSI()
    assert np.isclose(jsi.pwl, 775e-9)
    assert np.isclose(jsi.tau, 2e-12)
    assert jsi.m == -1

def test_jsi_getplots_basic():
    jsi = JSI()
    ri = RefractiveIndex()
    material = "KTP"
    refidxfunc = [
        ri.getSingleIDX(material, "X", "kato"),
        ri.getSingleIDX(material, "Y", "kato"),
        ri.getSingleIDX(material, "Z", "kato")
    ]
    
    # We use very small resolution (e.g. 10x10) for testing speed
    pumpwl = 775e-9
    sigrange = np.linspace(1500e-9, 1600e-9, 10)
    idrange = np.linspace(1500e-9, 1600e-9, 10)
    tau = 2e-12
    T = 25.0
    PP = 46.2e-6
    L = 30e-3
    qpmorder = 1
    
    f = Filters()
    spectralfilters = [
        f.getFilterFunction("None", 1, 1),
        f.getFilterFunction("None", 1, 1)
    ]
    
    # PE, PM, JS
    PE, PM, JS = jsi.getplots(pumpwl, sigrange, idrange, tau, T, PP, L, refidxfunc, qpmorder, spectralfilters, True, "Gaussian", 0, False, False, 10e-3, 10e-3, 10e-3, 1e-3, 1e-3, 1e-3)
    
    assert PE.shape == (10, 10)
    assert PM.shape == (10, 10)
    assert JS.shape == (10, 10)
    
    # JSI values should be non-negative
    assert np.all(JS >= 0)
