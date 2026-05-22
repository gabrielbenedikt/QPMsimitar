import pytest
import numpy as np
from PMC import PMC
from RefractiveIndex import RefractiveIndex

def test_pmc_initialization():
    pmc = PMC()
    assert np.isclose(pmc.lp, 775e-9)
    assert np.isclose(pmc.PP, 46.2e-6)
    assert pmc.m == 1

def test_getSI_wl_varT():
    pmc = PMC()
    ri = RefractiveIndex()
    material = "KTP"
    nxfunc = ri.getSingleIDX(material, "X", "kato")
    nyfunc = ri.getSingleIDX(material, "Y", "kato")
    nzfunc = ri.getSingleIDX(material, "Z", "kato")
    refidxfunc = [nxfunc, nyfunc, nzfunc]
    
    pumpwl = 775e-9
    polingp = 46.2e-6
    Trange = np.array([25.0, 30.0, 35.0])
    qpmorder = 1
    
    sigwl, idwl, Tcp = pmc.getSI_wl_varT(pumpwl, polingp, Trange, refidxfunc, qpmorder)
    
    assert len(sigwl) == 3
    assert len(idwl) == 3
    assert Tcp > 0
    # The wavelengths must satisfy energy conservation
    assert np.allclose(1/pumpwl, 1/sigwl + 1/idwl)

def test_getSI_wl_varPP():
    pmc = PMC()
    ri = RefractiveIndex()
    material = "KTP"
    refidxfunc = [
        ri.getSingleIDX(material, "X", "kato"),
        ri.getSingleIDX(material, "Y", "kato"),
        ri.getSingleIDX(material, "Z", "kato")
    ]
    
    pumpwl = 775e-9
    PPrange = np.array([45e-6, 46.2e-6, 47e-6])
    T = 25.0
    qpmorder = 1
    
    sigwl, idwl, PPcp = pmc.getSI_wl_varPP(pumpwl, PPrange, T, refidxfunc, qpmorder)
    
    assert len(sigwl) == 3
    assert len(idwl) == 3
    assert PPcp > 0
