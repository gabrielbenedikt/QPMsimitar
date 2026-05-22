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
    'pump_cw_bw': 0,
    'delay_range': 2e-12,
    'hom_resolution': 250,
    'hom_phase': 0,
    'jsi_resolution': 100,
    'jsi_wl_range': 100e-9
}

try:
    BaseComputeParams(**raw_params)
    print("Success")
except Exception as e:
    print(e)
