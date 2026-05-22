from compute.schemas import BaseComputeParams

raw_params = {
    'material': 'KTP',
    'nx_paper': 'kato',
    'ny_paper': 'kato',
    'nz_paper': 'kato',
    'pump_wl': 775e-9,
    'poling_period': 46.2e-6,
    'temperature': 20.0,
    'qpm_order': 1,
    'pulsewidth': 1e-12,
    'pump_shape': 'Gauss',
    'pump_cw_bw': 0.0,
    'signal_filter_type': 'None',
    'signal_filter_center_wl': 1550e-9,
    'signal_filter_fwhm': 3e-9,
    'idler_filter_type': 'None',
    'idler_filter_center_wl': 1550e-9,
    'idler_filter_fwhm': 3e-9,
    'focusing_enable': False,
    'fibre_coupling_enable': False,
    'focallength_pump': 10e-3,
    'focallength_signal': 10e-3,
    'focallength_idler': 10e-3,
    'beamdiameter_pump': 1e-3,
    'beamdiameter_signal': 1e-3,
    'beamdiameter_idler': 1e-3,
    'delay_range': 2e-12,
    'hom_resolution': 250.0, # Wait! is resolution a float?
    'hom_phase': 0.0,
    'jsi_resolution': 100.0,
    'jsi_wl_range': 100e-9
}

try:
    BaseComputeParams(**raw_params)
    print("Success")
except Exception as e:
    print(e)
