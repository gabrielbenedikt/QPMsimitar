import pytest
import httpx
import os
import struct
import msgpack
from compute.server import app

from compute.serialization import serialize, deserialize

@pytest.fixture(autouse=True)
def setup_env():
    # Force a known token for tests
    os.environ["QPMSIMITAR_API_TOKEN"] = "test-token"
    yield
    os.environ.pop("QPMSIMITAR_API_TOKEN", None)

@pytest.mark.asyncio
async def test_server_materials_unauthorized():
    async with httpx.AsyncClient(app=app, base_url="http://test") as client:
        response = await client.get("/materials")
        assert response.status_code == 403

@pytest.mark.asyncio
async def test_server_materials_authorized():
    async with httpx.AsyncClient(app=app, base_url="http://test") as client:
        response = await client.get("/materials", headers={"X-API-Token": "test-token"})
        assert response.status_code == 200
        mats = response.json()
        assert isinstance(mats, list)
        assert "KTP" in mats

@pytest.mark.asyncio
async def test_compute_endpoint():
    # We send compressed msgpack bytes and receive streaming compressed msgpack
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
    payload = serialize(params)
    
    async with httpx.AsyncClient(app=app, base_url="http://test") as client:
        async with client.stream("POST", "/compute/compute_pmc_vs_T", 
                                 content=payload, 
                                 headers={"X-API-Token": "test-token"}) as response:
            assert response.status_code == 200
            
            body_bytes = await response.aread()
            assert len(body_bytes) >= 4
            length = struct.unpack("!I", body_bytes[:4])[0]
            msg = deserialize(body_bytes[4:4+length])
            if msg['type'] == 'result':
                assert 'T_range' in msg['value']
                assert len(msg['value']['T_range']) > 0
            elif msg['type'] == 'error':
                pytest.fail(f"Compute error: {msg['value']}")
