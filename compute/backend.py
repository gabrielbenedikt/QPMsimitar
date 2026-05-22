"""
Backend abstraction for QPMsimitar compute.

ComputeBackend ABC defines the interface. LocalBackend runs in-process.
RemoteBackend (Phase 2) will send requests over the network.
"""

from abc import ABC, abstractmethod
from typing import Any, Callable, Optional


class ProgressCallback:
    """Simple progress reporter. GUI can poll or connect signals."""
    def __init__(self):
        self.progress = 0.0      # 0.0 to 1.0
        self.message = ""
        self._cancelled = False
        self._on_progress: Optional[Callable] = None

    def set_callback(self, fn: Callable):
        """Set a function to call on each progress update: fn(progress, message)."""
        self._on_progress = fn

    def update(self, progress: float, message: str = ""):
        self.progress = progress
        self.message = message
        if self._on_progress:
            self._on_progress(progress, message)

    def cancel(self):
        self._cancelled = True

    @property
    def is_cancelled(self):
        return self._cancelled


class ComputeBackend(ABC):
    """Abstract interface between GUI and numerics."""

    @abstractmethod
    def call(self, method: str, params: dict,
             progress: Optional[ProgressCallback] = None) -> dict:
        """Execute a named computation with given parameters.

        Args:
            method: Name of the compute method (e.g. 'compute_pmc_vs_T')
            params: Dict of parameters (scalars, strings, numpy arrays)
            progress: Optional progress callback for long-running tasks

        Returns:
            Result dict with numpy arrays and scalars.
        """
        ...

    def get_available_materials(self) -> list:
        """Return list of available crystal materials."""
        ...

    def get_available_refractive_indices(self, material: str) -> list:
        """Return available refractive indices for a material."""
        ...


class LocalBackend(ComputeBackend):
    """Runs ComputeEngine in the same process."""

    def __init__(self):
        from compute.engine import ComputeEngine
        self._engine = ComputeEngine()

    def call(self, method: str, params: dict,
             progress: Optional[ProgressCallback] = None) -> dict:
        fn = getattr(self._engine, method, None)
        if fn is None:
            raise ValueError(f"Unknown compute method: {method}")
        if progress is not None:
            params['_progress'] = progress
        return fn(params)

    def get_available_materials(self) -> list:
        return self._engine.get_available_materials()

    def get_available_refractive_indices(self, material: str) -> list:
        return self._engine.get_available_refractive_indices(material)

class RemoteBackend(ComputeBackend):
    """Executes ComputeEngine methods on a remote FastAPI server."""

    def __init__(self, server_url: str, api_token: str = "", client_cert: str = "",
                 client_key: str = "", ca_cert: str = "", verify_ssl: bool = True):
        self.server_url = server_url.rstrip('/')
        import httpx
        import ssl
        
        headers = {}
        if api_token:
            headers["X-API-Token"] = api_token

        # Configure mTLS or cert verification
        verify = verify_ssl
        if verify_ssl:
            if ca_cert and client_cert and client_key:
                context = ssl.create_default_context(cafile=ca_cert)
                context.load_cert_chain(certfile=client_cert, keyfile=client_key)
                verify = context
            elif ca_cert:
                verify = ca_cert
        else:
            # If verify_ssl is False, we pass False to httpx which disables verification.
            # However, if client certs are still provided, we still need to load them to do client auth.
            if client_cert and client_key:
                # We create a context, but disable hostname/cert checking on it.
                context = ssl.create_default_context()
                context.check_hostname = False
                context.verify_mode = ssl.CERT_NONE
                context.load_cert_chain(certfile=client_cert, keyfile=client_key)
                verify = context

        # Configure httpx client with some sane timeouts for long-running compute tasks
        self.client = httpx.Client(timeout=httpx.Timeout(None), headers=headers, verify=verify)

    def call(self, method: str, params: dict,
             progress: Optional[ProgressCallback] = None) -> dict:
        from compute.serialization import serialize, deserialize
        import struct

        url = f"{self.server_url}/compute/{method}"
        payload = serialize(params)
        
        # Read streaming response
        with self.client.stream("POST", url, content=payload, headers={"Content-Type": "application/x-msgpack"}) as response:
            if response.status_code == 400:
                error_text = response.read().decode('utf-8', errors='ignore')
                raise RuntimeError(f"Server Validation Error: {error_text}")
            response.raise_for_status()
            
            buf = bytearray()
            for chunk in response.iter_bytes():
                buf.extend(chunk)
                while len(buf) >= 4:
                    msg_len = struct.unpack("!I", buf[:4])[0]
                    if len(buf) < 4 + msg_len:
                        break # wait for more data
                        
                    msg_bytes = bytes(buf[4:4+msg_len])
                    buf = buf[4+msg_len:]
                    
                    msg = deserialize(msg_bytes)
                    if msg['type'] == 'progress':
                        if progress:
                            progress.update(msg['value'] / 100.0) # assuming engine sends 0-100
                    elif msg['type'] == 'result':
                        if progress:
                            progress.update(1.0)
                        return msg['value']
                    elif msg['type'] == 'error':
                        raise RuntimeError(f"Remote compute error: {msg['value']}")
                        
        raise RuntimeError("Server closed connection without returning a result.")

    def get_available_materials(self) -> list:
        response = self.client.get(f"{self.server_url}/materials")
        response.raise_for_status()
        return response.json()

    def get_available_refractive_indices(self, material: str) -> list:
        response = self.client.get(f"{self.server_url}/refractive_indices/{material}")
        response.raise_for_status()
        return response.json()
