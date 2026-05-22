import pytest
import numpy as np
from compute.engine import ComputeEngine

def test_engine_available_materials():
    engine = ComputeEngine()
    mats = engine.get_available_materials()
    assert isinstance(mats, list)
    assert "KTP" in mats
    assert "PPKTP" in mats

def test_compute_pmc_vs_T():
    engine = ComputeEngine()
    params = {
        'material': 'KTP',
        'nx_paper': 'kato',
        'ny_paper': 'kato',
        'nz_paper': 'kato',
        'pump_wl': 775e-9,
        'poling_period': 46.2e-6,
        'T_min': 20.0,
        'T_max': 40.0,
        'qpm_order': 1
    }
    
    result = engine.compute_pmc_vs_T(params)
    
    assert 'T_range' in result
    assert 'signal_wl' in result
    assert 'idler_wl' in result
    assert 'Tcp' in result
    
    assert len(result['T_range']) > 0
    assert len(result['signal_wl']) == len(result['T_range'])
    assert result['Tcp'] > 0
    
    # Energy conservation check for a single random point
    idx = len(result['signal_wl']) // 2
    swl = result['signal_wl'][idx]
    iwl = result['idler_wl'][idx]
    assert np.isclose(1/params['pump_wl'], 1/swl + 1/iwl)

def test_compute_refractive_index_vs_T():
    engine = ComputeEngine()
    params = {
        'material': 'KTP',
        'T_min': 20.0,
        'T_max': 40.0,
        'wavelength': 775e-9,
        'selected_indices': [('KTP', 'X', 'kato')]
    }
    
    result = engine.compute_refractive_index_vs_T(params)
    
    assert 'T_range' in result
    assert 'n_values' in result
    assert 'KTP:X:kato' in result['n_values']
    assert len(result['n_values']['KTP:X:kato']) == len(result['T_range'])
