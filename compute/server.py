import asyncio
import struct
import inspect
from fastapi import FastAPI, Request, HTTPException, Depends
from fastapi.responses import StreamingResponse

from compute.engine import ComputeEngine
from compute.serialization import serialize, deserialize
from compute.security import verify_token

app = FastAPI(title="QPMsimitar Compute Server", dependencies=[Depends(verify_token)])
engine = ComputeEngine()

@app.get("/materials")
def get_materials():
    return engine.get_available_materials()

@app.get("/refractive_indices/{material}")
def get_refractive_indices(material: str):
    return engine.get_available_refractive_indices(material)

@app.post("/compute/{method_name}")
async def compute_endpoint(method_name: str, request: Request):
    """
    Executes a method on the ComputeEngine.
    Expects a msgpack serialized dictionary of parameters in the request body.
    Returns a stream of length-prefixed msgpack messages.
    """
    if not hasattr(engine, method_name):
        raise HTTPException(status_code=404, detail=f"Method {method_name} not found in ComputeEngine")
        
    body = await request.body()
    try:
        raw_params = deserialize(body)
        from compute.schemas import BaseComputeParams
        # Validate structure (raises ValidationError on failure)
        BaseComputeParams(**raw_params)
        params = raw_params
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Invalid payload or schema: {e}")

    loop = asyncio.get_running_loop()
    q = asyncio.Queue()

    def progress_cb(percent: int):
        loop.call_soon_threadsafe(q.put_nowait, {'type': 'progress', 'value': percent})

    def worker():
        try:
            method = getattr(engine, method_name)
            # Inject progress callback into the parameters dictionary
            params['progress_callback'] = progress_cb
                
            result = method(params)
            loop.call_soon_threadsafe(q.put_nowait, {'type': 'result', 'value': result})
        except Exception as e:
            loop.call_soon_threadsafe(q.put_nowait, {'type': 'error', 'value': str(e)})

    # Start the computation in a background thread
    task = loop.run_in_executor(None, worker)

    async def event_stream():
        while True:
            msg = await q.get()
            
            # Serialize the message and prefix with a 4-byte network-order length
            payload = serialize(msg)
            header = struct.pack("!I", len(payload))
            yield header + payload
            
            if msg['type'] in ('result', 'error'):
                break

    return StreamingResponse(event_stream(), media_type="application/x-msgpack")
