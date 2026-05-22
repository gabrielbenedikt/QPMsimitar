import asyncio
import httpx
from compute.serialization import serialize

async def run():
    params = {
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
        'jsi_wl_range': 100e-9,
        'focusing_enable': False,
        'fibre_coupling_enable': False,
        'signal_filter_type': 'None',
        'idler_filter_type': 'None',
    }
    payload = serialize(params)
    async with httpx.AsyncClient(base_url="https://10.42.1.15:8443", verify=False) as client:
        response = await client.post("/compute/compute_hom_interference", content=payload, headers={"X-API-Token": "test-token"})
        print(response.status_code)
        print(response.text)

asyncio.run(run())
